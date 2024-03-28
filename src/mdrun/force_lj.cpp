#include "force_lj.h"
#include "luisa/dsl/builtin.h"
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

ForceLJ::~ForceLJ() {}

void ForceLJ::setup(Stream &stream, Device &device, Atom &atom) {
  f = device.create_buffer<float3>(atom.nlocal);
  f_delta = device.create_buffer<float3>(atom.nlocal);
  eng_vdwl_buf = device.create_buffer<float>(atom.nlocal);
  virial_buf = device.create_buffer<float>(atom.nlocal);
}

void ForceLJ::setup_shader(Device &device, Atom &atom, Neighbor &neighbor) {
  nlocal = atom.nlocal;
  Kernel1D compute_lj_kernel = [&]() noexcept {
    auto i = dispatch_x();
    Float3 f_ = make_float3(0.f);

    auto jnum = neighbor.numneigh->read(i);
    Float3 tmpi = atom.x->read(i);
    Float xtmp = tmpi[0];
    Float ytmp = tmpi[1];
    Float ztmp = tmpi[2];
    Int type_i = atom.type->read(i);

    Float test = 0.f;
    Float eng_vdwl_ = 0.f;
    Float virial_ = 0.f;
    Int target = 2000;
    // Int target = 24;
    Int coord = 4690;
    Int min = 0;
    Float min_r = 100.f;
    $if(i == target) {
      device_log("x: {} ", xtmp);
      device_log("y: {} ", ytmp);
      device_log("z: {} ", ztmp);
    };
    $for(k, jnum) {
      Int j = neighbor.neighbors->read(make_uint2(i, k.cast<uint>()))[0];
      Float3 tmpj = atom.x->read(j);
      Float delx = xtmp - tmpj[0];
      Float dely = ytmp - tmpj[1];
      Float delz = ztmp - tmpj[2];
      Int type_j = atom.type->read(j);
      Float rsq = delx * delx + dely * dely + delz * delz;
      // $if(rsq < 0.8f) { device_log("too small rsq of {}", rsq); };

      Int type_ij = type_i * ntypes + type_j;

      // $if(i == target & j == coord) {
      //   device_log("x: {} and {}", xtmp, tmpj[0]);
      //   device_log("y: {} and {}", ytmp, tmpj[1]);
      //   device_log("z: {} and {}", ztmp, tmpj[2]);
      //   device_log("rsq {} : ", rsq);
      // };

      // $if(i == target) {
      //   device_log("to {}, rsq: {}", j, rsq);
      //   //   device_log("cut: {}", cutforcesq->read(type_ij));
      // };

      $if(rsq < cutforcesq->read(type_ij)) {
        $if(rsq < min_r) {
          min_r = rsq;
          min = j;
        };
        Float sr2 = 1.0f / rsq;
        Float sr6 = sr2 * sr2 * sr2 * sigma6->read(type_ij);
        Float force =
            48.0f * pow(1.f / (delx * delx + dely * dely + delz * delz), 3.f) *
            (pow(1.f / (delx * delx + dely * dely + delz * delz), 3.f) - 0.5f) *
            1.f / (delx * delx + dely * dely + delz * delz) *
            epsilon->read(type_ij);

        f_[0] += delx * force;
        f_[1] += dely * force;
        f_[2] += delz * force;
        // $if(i == target) {
        // $if(i == target & j == coord) {
        //   device_log("force with {} : ({},{},{})", j, delx * force,
        //              dely * force, delz * force);
        // };

        if (evflag) {
          eng_vdwl_ += (4.0f * sr6 * (sr6 - 1.0f)) * epsilon->read(type_ij);
          virial_ += (delx * delx + dely * dely + delz * delz) * force;
        }
      };
    };
    eng_vdwl_buf->write(i, eng_vdwl_);
    virial_buf->write(i, virial_);

    $if(i == target) {
      //   device_log("min: {}", min);
      //   device_log("min_rsq: {}", min_r);
      device_log("force: ({},{},{})", f_[0], f_[1], f_[2]);
    };
    f->write(i, f_);
  };

  Kernel1D gather_ev_kernel = [&]() noexcept {
    Float eng_vdwl_ = 0.f;
    Float virial_ = 0.f;
    $for(i, atom.nlocal) {
      eng_vdwl_ += eng_vdwl_buf->read(i);
      virial_ += virial_buf->read(i);
    };
    eng_vdwl->write(0, eng_vdwl_);
    virial->write(0, virial_);
  };

  auto o = luisa::compute::ShaderOption{.enable_fast_math = false};
  compute_shader = device.compile(compute_lj_kernel, o);
  gather_ev_shader = device.compile(gather_ev_kernel, o);
}
void ForceLJ::compute(Stream &stream) {
  stream << compute_shader().dispatch(nlocal) << synchronize();
  // if (evflag) {
  stream << gather_ev_shader().dispatch(1) << synchronize();
  // }
}
