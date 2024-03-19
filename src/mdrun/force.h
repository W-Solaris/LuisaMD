#ifndef FORCE_H_
#define FORCE_H_

#include <luisa/luisa-compute.h>

#include "atom.h"
#include "lj.h"
#include "luisa/runtime/bindless_array.h"
#include "neighbor.h"
using namespace luisa;
using namespace luisa::compute;

class Force {
 public:
  float cutforce;
  Buffer<float> cutforcesq;
  Buffer<float> eng_vdwl;
  Buffer<float> virial;
  float mass;
  int evflag;
  int ntypes;

  Force(){};
  virtual ~Force(){};
  // virtual void setup(){}; // no need
  virtual void finalise(){};
  virtual void compute(Atom &, Neighbor &, int){};

  Buffer<float> epsilon, sigma6, sigma;  // Parameters for LJ only
  float epsilon_scalar, sigma_scalar;

  ForceStyle style;

  float cutforcesq_s[MAX_STACK_TYPES * MAX_STACK_TYPES];
  float epsilon_s[MAX_STACK_TYPES * MAX_STACK_TYPES];
  float sigma6_s[MAX_STACK_TYPES * MAX_STACK_TYPES];

 protected:
  int nlocal;
  int nall;

  // Buffer<int> numneigh;  // # of neighbors for each atom
  // Image<int> neighbors;  // array of neighbors of each atom

  Buffer<float3> xf;
  Buffer<float3> f;
  Buffer<float3> f_a;
  Buffer<int> type;

  int me;
};
#endif