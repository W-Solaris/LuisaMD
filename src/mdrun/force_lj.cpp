#include "force_lj.h"
#include "math.h"
#include "stdio.h"

#ifndef VECTORLENGTH
#define VECTORLENGTH 4
#endif

ForceLJ::ForceLJ(int ntypes_) {
  cutforce = 0.0;
  use_oldcompute = 0;
  reneigh = 1;
  style = FORCELJ;
  ntypes = ntypes_;

  BufferFloat d_cut = device.create_buffer<float>(ntypes * ntypes);
  std::vector<float> h_cut(ntypes * ntypes);
  cutforcesq = d_cut;

  float_1d_view_type d_epsilon("ForceLJ::epsilon", ntypes * ntypes);
  std::vector<float> h_epsilon(ntypes * ntypes);
  epsilon = d_epsilon;

  float_1d_view_type d_sigma6("ForceLJ::sigma6", ntypes * ntypes);
  std::vector<float> h_sigma6(ntypes * ntypes);
  sigma6 = d_sigma6;

  float_1d_view_type d_sigma("ForceLJ::sigma", ntypes * ntypes);
  std::vector<float> h_sigma(ntypes * ntypes);
  sigma = d_sigma;

  for (int i = 0; i < ntypes * ntypes; i++) {
    h_cut[i] = 0.0;
    h_epsilon[i] = 1.0;
    h_sigma6[i] = 1.0;
    h_sigma[i] = 1.0;
    if (i < MAX_STACK_TYPES * MAX_STACK_TYPES) {
      epsilon_s[i] = 1.0;
      sigma6_s[i] = 1.0;
    }
  }

  Kokkos::deep_copy(d_cut, h_cut);
  Kokkos::deep_copy(d_epsilon, h_epsilon);
  Kokkos::deep_copy(d_sigma6, h_sigma6);
  Kokkos::deep_copy(d_sigma, h_sigma);
}

ForceLJ::~ForceLJ() {}

void ForceLJ::setup() {
  float_1d_view_type d_cut("ForceLJ::cutforcesq", ntypes * ntypes);
  std::vector<float> h_cut(ntypes * ntypes);
  cutforcesq = d_cut;

  for (int i = 0; i < ntypes * ntypes; i++) {
    h_cut[i] = cutforce * cutforce;
    if (i < MAX_STACK_TYPES * MAX_STACK_TYPES)
      cutforcesq_s[i] = cutforce * cutforce;
  }

  Kokkos::deep_copy(d_cut, h_cut);
}

void ForceLJ::compute(Atom &atom, Neighbor &neighbor) {
  eng_vdwl = 0;
  virial = 0;

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABEL_ROCM)
  const int host_device = 0;
#else
  const int host_device = 1;
#endif

  nlocal = atom.nlocal;
  nall = atom.nlocal + atom.nghost;

  x = atom.x;
  f_a = atom.f;
  f = atom.f;
  type = atom.type;

  neighbors = neighbor.neighbors;
  numneigh = neighbor.numneigh;

  // clear force on own and ghost atoms

  Kokkos::deep_copy(f, 0.0);

  /* switch to correct compute */
  eng_vdwl = 0;
  virial = 0;

  // loop over all neighbors of my atoms
  // store force on both atoms i and j

  for (int i = 0; i < nlocal; i++) {
    const int jnum = numneigh[i];
    const md_float xtmp = x(i, 0);
    const md_float ytmp = x(i, 1);
    const md_float ztmp = x(i, 2);
    const int type_i = type[i];

    for (int k = 0; k < jnum; k++) {
      const int j = neighbors(i, k);
      const md_float delx = xtmp - x(j, 0);
      const md_float dely = ytmp - x(j, 1);
      const md_float delz = ztmp - x(j, 2);
      int type_j = type[j];
      const md_float rsq = delx * delx + dely * dely + delz * delz;

      const int type_ij = type_i * ntypes + type_j;

      if (rsq < cutforcesq(type_ij)) {
        const md_float sr2 = 1.0 / rsq;
        const md_float sr6 = sr2 * sr2 * sr2 * sigma6(type_ij);
        const md_float force =
            48.0 * sr6 * (sr6 - 0.5) * sr2 * epsilon(type_ij);
        f(i, 0) += delx * force;
        f(i, 1) += dely * force;
        f(i, 2) += delz * force;
        f(j, 0) -= delx * force;
        f(j, 1) -= dely * force;
        f(j, 2) -= delz * force;

        if (evflag) {
          eng_vdwl += (4.0 * sr6 * (sr6 - 1.0)) * epsilon(type_ij);
          virial += (delx * delx + dely * dely + delz * delz) * force;
        }
      }
    }
  }
  .
}

、