#include "integrate.h"
#include "luisa/dsl/builtin.h"
#include "luisa/dsl/sugar.h"
#include "math.h"
#include "stdio.h"

Integrate::Integrate() { sort_every = 20; }
Integrate::~Integrate() {}

void Integrate::setup() { dtforce = 0.5 * dt; }

void Integrate::setup_shader(Device &device, Atom &atom, Force *force) {
  Kernel1D initial_integrate_kernel = [&]() noexcept {
    auto i = dispatch_x();
    Float3 v_ = atom.v->read(i);
    Float3 f_ = force->f->read(i);
    // $if(i == 2000) { device_log("f: ({}, {}, {})", f_[0], f_[1], f_[2]); };
    v_[0] += dtforce * f_[0];
    v_[1] += dtforce * f_[1];
    v_[2] += dtforce * f_[2];
    Float3 x_ = atom.x->read(i);
    x_[0] += dt * v_[0];
    x_[1] += dt * v_[1];
    x_[2] += dt * v_[2];
    atom.v->write(i, v_);
    atom.x->write(i, x_);
  };
  auto o = luisa::compute::ShaderOption{.enable_fast_math = false};
  initial_integrate = device.compile(initial_integrate_kernel, o);

  Kernel1D final_integrate_kernel = [&]() noexcept {
    auto i = dispatch_x();
    Float3 v_ = atom.v->read(i);
    Float3 f_ = force->f->read(i);
    v_[0] += dtforce * f_[0];
    v_[1] += dtforce * f_[1];
    v_[2] += dtforce * f_[2];
    atom.v->write(i, v_);
  };
  final_integrate = device.compile(final_integrate_kernel, o);
}

void Integrate::run(Stream &stream, Atom &atom, Force *force, Comm &comm,
                    Neighbor &neighbor, Thermo &thermo) {
  mass = atom.mass;
  dtforce = dtforce / mass;
  nlocal = atom.nlocal;

  int next_sort = sort_every > 0 ? sort_every : ntimes + 1;

  // Buffer<float> d_max = device.create_buffer<float>(1);
  // Buffer<float> d_record = device.create_buffer<float>(nlocal);

  nlocal = atom.nlocal;
  for (int n = 0; n < ntimes; n++) {
    stream << initial_integrate().dispatch(nlocal) << synchronize();
    if ((n + 1) % neighbor.every != 0) {
      comm.communicate(stream);
    } else {
      atom.pbc(stream);
      // WIP
      if (n + 1 >= next_sort) {
        atom.sort(neighbor);
        next_sort += sort_every;
      }
      comm.borders(stream);
    }

    neighbor.build(stream);

    force->evflag = (n + 1) % thermo.nstat == 0;
    force->compute(stream);

    stream << final_integrate().dispatch(nlocal) << synchronize();

    if (thermo.nstat)
      thermo.compute(stream, n + 1);
  }
}
