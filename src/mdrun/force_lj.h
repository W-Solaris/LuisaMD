#ifndef FORCELJ_H
#define FORCELJ_H

#include "atom.h"
#include "force.h"
#include "neighbor.h"
#include "types.h"

class ForceLJ : Force {
 public:
  template <int EVFLAG, int GHOST_NEWTON, int STACK_PARAMS>
  struct TagComputeHalfNeighThread {
    enum { evflag = EVFLAG };
    enum { ghost_newton = GHOST_NEWTON };
    enum { stack_params = STACK_PARAMS };
  };
  template <int EVFLAG, int STACK_PARAMS>
  struct TagComputeFullNeigh {
    enum { evflag = EVFLAG };
    enum { stack_params = STACK_PARAMS };
  };

  ForceLJ(int ntypes_);
  virtual ~ForceLJ();
  void setup();
  void compute(Atom &, Neighbor &, int);

 protected:
  template <int EVFLAG>
  void compute_original(Atom &, Neighbor &, int);

 public:
};

#endif
