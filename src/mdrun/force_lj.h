#ifndef FORCELJ_H
#define FORCELJ_H

#include "atom.h"
#include "force.h"
#include "neighbor.h"
#include "types.h"

class ForceLJ : Force {
public:
  ForceLJ(Stream &stream, Device &device, In &in, int ntypes_);
  virtual ~ForceLJ();
  void setup(Stream &, Device &, Atom &);
  // void setup(Stream &stream, Device &device);
  void setup_shader(Device &, Atom &, Neighbor &);
  void compute(Stream &);
};

#endif
