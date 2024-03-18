#ifndef ATOM_H
#define ATOM_H

#include <luisa/luisa-compute.h>

#include "string.h"
#include "types.h"
using namespace luisa;
using namespace luisa::compute;

class Neighbor;
struct Box {
  md_float xprd, yprd, zprd;
  md_float xlo, xhi;
  md_float ylo, yhi;
  md_float zlo, zhi;
};

class Atom {
public:
  typedef int value_type;
  int natoms;
  int nlocal, nghost;
  int nmax;

  Buffer<float3> x;
  Buffer<float3> v;
  Buffer<float3> f;
  std::vector<float3> h_x;
  std::vector<float3> h_v;
  std::vector<float3> h_f;

  int ntypes;
  Buffer<int> type, new_type, old_type;
  std::vector<int> h_type;

  Buffer<float3> xold, new_x, new_v, old_x, old_v;

  md_float virial, mass;

  int comm_size, reverse_size, border_size;

  Box box;

  Atom(int ntypes_);
  ~Atom();

  void addatom(md_float, md_float, md_float, md_float, md_float, md_float);

  // private:
  //   int_1d_view_type binpos;
  //   int_2d_view_type bins;
  //   x_view_type x_copy;
  //   x_view_type v_copy;
  //   int_1d_view_type type_copy;
  int copy_size;

  //   float_1d_view_type buf;
  //   int_1d_view_type list;
  //   int pbc_flags[4];
  //   int first;
};

// void Atom::copy(int i, int j) const {
//   x(j, 0) = x(i, 0);
//   x(j, 1) = x(i, 1);
//   x(j, 2) = x(i, 2);
//   v(j, 0) = v(i, 0);
//   v(j, 1) = v(i, 1);
//   v(j, 2) = v(i, 2);
//   type[j] = type[i];
// }

// void Atom::operator()(TagAtomPBC, const int &i) const {
//   if (x(i, 0) < 0.0)
//     x(i, 0) += box.xprd;

//   if (x(i, 0) >= box.xprd)
//     x(i, 0) -= box.xprd;

//   if (x(i, 1) < 0.0)
//     x(i, 1) += box.yprd;

//   if (x(i, 1) >= box.yprd)
//     x(i, 1) -= box.yprd;

//   if (x(i, 2) < 0.0)
//     x(i, 2) += box.zprd;

//   if (x(i, 2) >= box.zprd)
//     x(i, 2) -= box.zprd;
// }

// void Atom::operator()(TagAtomSort, const int &i, int &sum, bool final) const
// {
//   sum += binpos[i];
//   if (final)
//     binpos[i] = sum;
// }

// void Atom::operator()(TagAtomSort, const int &mybin) const {
//   const int start = mybin > 0 ? binpos[mybin - 1] : 0;
//   const int count = binpos[mybin] - start;
//   for (int k = 0; k < count; k++) {
//     const int new_i = start + k;
//     const int old_i = bins(mybin, k);
//     new_x(new_i, 0) = old_x(old_i, 0);
//     new_x(new_i, 1) = old_x(old_i, 1);
//     new_x(new_i, 2) = old_x(old_i, 2);
//     new_v(new_i, 0) = old_v(old_i, 0);
//     new_v(new_i, 1) = old_v(old_i, 1);
//     new_v(new_i, 2) = old_v(old_i, 2);
//     new_type[new_i] = old_type[old_i];
//   }
// }
#endif
