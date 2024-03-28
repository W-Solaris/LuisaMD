#include "atom.h"
#include "luisa/dsl/func.h"
#include "neighbor.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#define DELTA 20000

Atom::Atom(Device &device, int ntypes_) {
  natoms = 0;
  nlocal = 0;
  nmax = 0;
  copy_size = 0;

  // comm_size = 3;
  // reverse_size = 3;
  // border_size = 4;

  mass = 1;

  ntypes = ntypes_;
}

Atom::~Atom() {}

void Atom::addatom(float x_in, float y_in, float z_in, float vx_in, float vy_in,
                   float vz_in) {
  // store at vector, copy to buffer later
  h_x.push_back(make_float3(x_in, y_in, z_in));
  h_v.push_back(make_float3(vx_in, vy_in, vz_in));
  h_type.push_back(rand() % ntypes);
  nlocal++;
}

void Atom::construct_buf(Stream &stream, Device &device) {
  int max_len = 7 * nlocal;
  // to save ghost atoms : 3 dimensions * 2 directions
  nghost = device.create_buffer<int>(1);
  x = device.create_buffer<float3>(max_len);
  v = device.create_buffer<float3>(nlocal);
  f = device.create_buffer<float3>(nlocal);
  type = device.create_buffer<int>(nlocal);
  xold = device.create_buffer<float3>(nlocal);
  stream << x.copy_from(h_x.data());
  stream << v.copy_from(h_v.data());
  stream << type.copy_from(h_type.data());
}

/* enforce PBC
   order of 2 tests is important to insure lo-bound <= coord < hi-bound
   even with round-off errors where (coord +/- epsilon) +/- period = bound */

void Atom::setup_shader(Device &device) {
  Kernel1D apply_pbc_kernel = [&]() noexcept {
    auto i = dispatch_x();
    Float3 x_ = x->read(i);
    $if(x_[0] < 0.f) { x_[0] += box.xlen; };
    $if(x_[0] >= box.xlen) { x_[0] -= box.xlen; };
    $if(x_[1] < 0.f) { x_[1] += box.ylen; };
    $if(x_[1] >= box.ylen) { x_[1] -= box.ylen; };
    $if(x_[2] < 0.f) { x_[2] += box.zlen; };
    $if(x_[2] >= box.zlen) { x_[2] -= box.zlen; };
    x->write(i, x_);
  };
  apply_pbc = device.compile(apply_pbc_kernel);
}

void Atom::pbc(Stream &stream) {
  stream << apply_pbc().dispatch(nlocal) << synchronize();
}

/* realloc a 2-d float array */

void Atom::sort(Neighbor &neighbor) {
  return;
  // WIP!!
  // neighbor.binatoms(*this, nlocal);

  // Kokkos::fence();

  // binpos = neighbor.bincount;
  // bins = neighbor.bins;

  // const int mbins = neighbor.mbins;

  // Kokkos::parallel_scan(Kokkos::RangePolicy<TagAtomSort>(0, mbins), *this);

  // if (copy_size < nmax) {
  //   x_copy = x_view_type("atom::x_copy", nmax);
  //   v_copy = x_view_type("atom::v_copy", nmax);
  //   type_copy = int_1d_view_type("atom::type_copy", nmax);
  //   copy_size = nmax;
  // }

  // new_x = x_copy;
  // new_v = v_copy;
  // new_type = type_copy;
  // old_x = x;
  // old_v = v;
  // old_type = type;

  // Kokkos::parallel_for(Kokkos::RangePolicy<TagAtomSort>(0, mbins), *this);
  // Kokkos::fence();

  // x_view_type x_tmp = x;
  // x_view_type v_tmp = v;
  // int_1d_view_type type_tmp = type;

  // x = x_copy;
  // v = v_copy;
  // type = type_copy;
  // x_copy = x_tmp;
  // v_copy = v_tmp;
  // type_copy = type_tmp;
}
