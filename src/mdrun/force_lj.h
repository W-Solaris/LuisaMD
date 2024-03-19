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
  void compute(Atom &, Neighbor &);
};

#endif
