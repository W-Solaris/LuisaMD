#include "integrate.h"
#include "luisa/dsl/builtin.h"
#include "luisa/dsl/sugar.h"
#include "math.h"
#include "stdio.h"

Integrate::Integrate() { sort_every = 20; }
Integrate::~Integrate() {}

void Integrate::setup() { dtforce = 0.5 * dt; }

void Integrate::initialIntegrate(Stream &stream, Device &device) {
  Kernel1D initial_integrate_kernel = [&]() noexcept {
    auto i = dispatch_x();
    Float3 v_ = v->read(i);
    Float3 f_ = f->read(i);
    v_[0] += dtforce * f_[0];
    v_[1] += dtforce * f_[1];
    v_[2] += dtforce * f_[2];
    Float3 x_ = x->read(i);
    x_[0] += dt * v_[0];
    x_[1] += dt * v_[1];
    x_[2] += dt * v_[2];
    v->write(i, v_);
    x->write(i, x_);
  };
  auto initial_integrate = device.compile(initial_integrate_kernel);
  stream << initial_integrate().dispatch(nlocal) << synchronize();
}

void Integrate::finalIntegrate(Stream &stream, Device &device) {
  Kernel1D final_integrate_kernel = [&]() noexcept {
    auto i = dispatch_x();
    Float3 v_ = v->read(i);
    Float3 f_ = f->read(i);
    v_[0] += dtforce * f_[0];
    v_[1] += dtforce * f_[1];
    v_[2] += dtforce * f_[2];
    v->write(i, v_);
  };
  auto final_integrate = device.compile(final_integrate_kernel);
  stream << final_integrate().dispatch(nlocal) << synchronize();
}

void Integrate::run(Stream &stream, Device &device, Atom &atom, Force *force,
                    Neighbor &neighbor, Thermo &thermo) {
  mass = atom.mass;
  dtforce = dtforce / mass;
  nlocal = atom.nlocal;

  int next_sort = sort_every > 0 ? sort_every : ntimes + 1;

  Buffer<float> d_max = device.create_buffer<float>(1);
  Buffer<float> d_record = device.create_buffer<float>(nlocal);

  nlocal = atom.nlocal;
  for (int n = 0; n < ntimes; n++) {
    initialIntegrate(stream, device);

    Kernel1D final_integrate_kernel = [&]() noexcept {
      auto i = dispatch_x();
      Float3 x_ = atom.x->read(i);
      Float3 x_old_ = atom.xold->read(i);
      Float dx = x_[0] - x_old_[0];

      $if(dx > atom.box.xlen) { dx -= atom.box.xlen; };
      $if(dx < -atom.box.xlen) { dx += atom.box.xlen; };

      Float dy = x_[1] - x_old_[1];
      $if(dy > atom.box.ylen) { dy -= atom.box.ylen; };
      $if(dy < -atom.box.ylen) { dy += atom.box.ylen; };

      Float dz = x_[2] - x_old_[2];
      $if(dz > atom.box.zlen) { dz -= atom.box.zlen; };
      $if(dz < -atom.box.zlen) { dz += atom.box.zlen; };

      d_record->write(i, dx * dx + dy * dy + dz * dz);
    };

    Kernel1D find_max_d_kernel = [&]() noexcept {
      d_max->write(0, 0.f);
      $for(i, nlocal) {
        Float tmp = d_record->read(i);
        $if(tmp > d_max->read(0)) { d_max->write(i, tmp); };
      };
      d_max->write(0, sqrt(d_max->read(0)));
    };

    auto final_integratete = device.compile(final_integrate_kernel);
    auto find_max_d = device.compile(find_max_d_kernel);
    stream << final_integratete().dispatch(nlocal) << synchronize();
    stream << find_max_d().dispatch(nlocal) << synchronize();

    atom.pbc(stream, device);

    // WIP
    if (n + 1 >= next_sort) {
      atom.sort(neighbor);
      next_sort += sort_every;
    }

    neighbor.build(stream, device, atom);

    force->evflag = (n + 1) % thermo.nstat == 0;
    force->compute(stream, device, atom, neighbor);

    finalIntegrate(stream, device);

    if (thermo.nstat)
      thermo.compute(stream, device, n + 1, atom, neighbor, force);
  }
}
