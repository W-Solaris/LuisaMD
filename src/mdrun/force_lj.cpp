#include "force_lj.h"
#include "math.h"
#include "stdio.h"

ForceLJ::ForceLJ(Stream &stream, Device &device, In &in, int ntypes_) {
  eng_vdwl = device.create_buffer<float>(1);
  virial = device.create_buffer<float>(1);

  style = FORCELJ;
  ntypes = ntypes_;
  cutforce = in.force_cut;

  cutforcesq = device.create_buffer<float>(ntypes * ntypes);
  std::vector<float> h_cut(ntypes * ntypes);

  epsilon = device.create_buffer<float>(ntypes * ntypes);
  std::vector<float> h_epsilon(ntypes * ntypes);

  sigma6 = device.create_buffer<float>(ntypes * ntypes);
  std::vector<float> h_sigma6(ntypes * ntypes);

  sigma = device.create_buffer<float>(ntypes * ntypes);
  std::vector<float> h_sigma(ntypes * ntypes);

  epsilon_scalar = in.epsilon;
  sigma_scalar = in.sigma;

  for (int i = 0; i < ntypes * ntypes; i++) {
    h_cut[i] = cutforce * cutforce;
    // if (i < MAX_STACK_TYPES * MAX_STACK_TYPES)
    //   cutforcesq_s[i] = cutforce * cutforce;
    h_epsilon[i] = in.epsilon;
    h_sigma[i] = in.sigma;
    h_sigma6[i] =
        in.sigma * in.sigma * in.sigma * in.sigma * in.sigma * in.sigma;

    // // TODO: epsilon_s should be GPU buffer!
    // if (i < MAX_STACK_TYPES * MAX_STACK_TYPES) {
    //   epsilon_s[i] = h_epsilon[i];
    //   sigma6_s[i] = h_sigma6[i];
    // }
  }

  // move to GPU
  stream << cutforcesq.copy_from(h_cut.data());
  stream << epsilon.copy_from(h_epsilon.data());
  stream << sigma.copy_from(h_sigma.data());
  stream << sigma6.copy_from(h_sigma6.data());
}

// void ForceLJ::setup(Stream &stream, Device &device) {
void ForceLJ::setup() {
  // pass
}

void ForceLJ::compute(Stream &stream, Device &device, Atom &atom,
                      Neighbor &neighbor) {
  nlocal = atom.nlocal;
  nall = atom.nlocal + atom.nghost;

  Kernel1D compute_lj_kernel = [&]() noexcept {
    eng_vdwl->write(0, 0.f);
    virial->write(0, 0.f);

    auto i = dispatch_x();
    auto jnum = neighbor.numneigh->read(i);
    Float xtmp = atom.x->read(i)[0];
    Float ytmp = atom.x->read(i)[1];
    Float ztmp = atom.x->read(i)[2];
    Int type_i = atom.type->read(i);

    $for(k, jnum) {
      Int j = neighbor.neighbors->read(make_uint2(i, k.cast<uint>()));
      Float delx = xtmp - atom.x->read(j)[0];
      Float dely = ytmp - atom.x->read(j)[1];
      Float delz = ztmp - atom.x->read(j)[2];
      Int type_j = atom.type->read(j);

      Float rsq = delx * delx + dely * dely + delz * delz;

      Int type_ij = type_i * ntypes + type_j;

      $if(rsq < cutforcesq->read(type_ij)) {
        Float sr2 = 1.0f / rsq;
        Float sr6 = sr2 * sr2 * sr2 * sigma6->read(type_ij);
        Float force = 48.0f * sr6 * (sr6 - 0.5f) * sr2 * epsilon->read(type_ij);
        Float3 tmp = f->read(i);
        tmp[0] += delx * force;
        tmp[1] += dely * force;
        tmp[2] += delz * force;
        f->write(i, tmp);
        tmp = f->read(j);
        tmp[0] -= delx * force;
        tmp[1] -= dely * force;
        tmp[2] -= delz * force;
        f->write(j, tmp);

        if (evflag) {
          eng_vdwl->write(0,
                          (4.0f * sr6 * (sr6 - 1.0f)) * epsilon->read(type_ij));
          virial->write(0, (delx * delx + dely * dely + delz * delz) * force);
        }
      };
    };
  };
  auto initial_integrate = device.compile(compute_lj_kernel);
  stream << initial_integrate().dispatch(nlocal) << synchronize();
}
