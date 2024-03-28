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
  Buffer<float> eng_vdwl_buf;
  Buffer<float> virial_buf;
  float mass;
  int evflag;
  int ntypes;

  Force(){};
  virtual ~Force(){};
  virtual void setup(Stream &, Device &, Atom &){};
  virtual void finalise(){};
  virtual void setup_shader(Device &, Atom &, Neighbor &){};
  virtual void compute(Stream &){};

  int use_sse;
  Buffer<float> epsilon, sigma6, sigma;  // Parameters for LJ only
  float epsilon_scalar, sigma_scalar;

  ForceStyle style;

  // float cutforcesq_s[MAX_STACK_TYPES * MAX_STACK_TYPES];
  // float epsilon_s[MAX_STACK_TYPES * MAX_STACK_TYPES];
  // float sigma6_s[MAX_STACK_TYPES * MAX_STACK_TYPES];

 public:
  int nlocal;
  int nall;

  // Buffer<float3> x;
  Buffer<float3> f;
  Buffer<float3> f_delta;
  Buffer<int> type;

 protected:
  Shader<1> compute_shader;
  Shader<1> gather_ev_shader;
};
#endif