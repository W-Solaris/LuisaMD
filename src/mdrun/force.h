#ifndef FORCE_H_
#define FORCE_H_

#include <luisa/luisa-compute.h>

#include "atom.h"
#include "lj.h"
#include "neighbor.h"
using namespace luisa;
using namespace luisa::compute;

class Force {
public:
  md_float cutforce;
  md_float eng_vdwl;
  md_float mass;
  int evflag;
  md_float virial;
  int ntypes;
  Force(){};
  virtual ~Force(){};
  virtual void setup(){};
  virtual void finalise(){};
  virtual void compute(Atom &, Neighbor &, int){};

  int use_sse;
  int use_oldcompute;
  int nthreads;
  int reneigh;

  Float epsilon, sigma6, sigma; // Parameters for LJ only
  md_float epsilon_scalar, sigma_scalar;

  ForceStyle style;

  md_float cutforcesq_s[MAX_STACK_TYPES * MAX_STACK_TYPES];
  md_float epsilon_s[MAX_STACK_TYPES * MAX_STACK_TYPES];
  md_float sigma6_s[MAX_STACK_TYPES * MAX_STACK_TYPES];

protected:
  int nlocal;
  int nall;

  Int numneigh; // # of neighbors for each atom

  Buffer<float3> f;
  Buffer<float3> f_a;
  Buffer<int> type;

  int me;
};
#endif