#include "thermo.h"
#include "force_lj.h"
#include "integrate.h"
#include "stdio.h"
#include "stdlib.h"

Thermo::Thermo() {}
Thermo::~Thermo() {}

void Thermo::setup(Stream &stream, Device &device, float rho_in,
                   Integrate &integrate, Atom &atom, int units) {
  rho = rho_in;
  ntimes = integrate.ntimes;
  nlocal = atom.nlocal;

  int maxstat;

  if (nstat == 0)
    maxstat = 2;
  else
    maxstat = ntimes / nstat + 2;

  steparr = device.create_buffer<int>(maxstat);
  tmparr = device.create_buffer<float>(maxstat);
  engarr = device.create_buffer<float>(maxstat);
  prsarr = device.create_buffer<float>(maxstat);

  t_act = device.create_buffer<float>(1);
  if (units == LJ) {
    mvv2e = 1.0;
    dof_boltz = atom.natoms * 3 - 3;
    t_scale = mvv2e / dof_boltz;
    p_scale = 1.0 / 3.0 / atom.box.xlen / atom.box.ylen / atom.box.zlen;
    e_scale = 0.5;

  } else if (units == METAL) {
    mvv2e = 1.036427e-04;
    dof_boltz = (atom.natoms * 3 - 3) * 8.617343e-05;
    t_scale = mvv2e / dof_boltz;
    p_scale = 1.602176e+06 / 3 / atom.box.xlen / atom.box.ylen / atom.box.zlen;
    e_scale = 524287.985533; // 16.0;
    integrate.dtforce /= mvv2e;
  }
}

void Thermo::compute(Stream &stream, Device &device, int iflag, Atom &atom,
                     Neighbor &neighbor, Force *force) {
  float t, eng, p;

  if (iflag > 0 && iflag % nstat)
    return;

  if (iflag == -1 && nstat > 0 && ntimes % nstat == 0)
    return;

  temperature(stream, device, atom);

  Kernel1D e_p_kernel = [&]() noexcept {
    Int natoms_ = atom.natoms;

    Int mstat_ = mstat->read(0);
    Float t_ = t_act->read(0) * t_scale;
    Float e_ = force->eng_vdwl->read(0) * e_scale / natoms_.cast<float>();
    Float p_ = (t_ * dof_boltz + force->virial->read(0)) * p_scale;

    tmparr->write(mstat_, t_);
    engarr->write(mstat_, e_);
    prsarr->write(mstat_, p_);

    Int istep = 0;
    $if(iflag == -1) { istep = ntimes; };
    $if(iflag == 0) { mstat_ = 0; };
    steparr->write(mstat_, istep);

    mstat_ += 1;
    mstat->write(0, mstat_);
  };

  auto energy_pressure = device.compile(e_p_kernel);
  stream << energy_pressure().dispatch(1) << synchronize();

  // TODO: output
}

void Thermo::temperature(Stream &stream, Device &device, Atom &atom) {

  Kernel1D temperature_kernel = [&]() noexcept {
    t_act->write(0, 0);
    Float t = 0.f;
    auto i = dispatch_x();
    Float3 v_ = atom.v->read(i);

    Float t_ = (v_[0] * v_[0] + v_[1] * v_[1] + v_[2] * v_[2]) * atom.mass;
    t_act->atomic(0).fetch_add(t_);
  };
  auto temperature = device.compile(temperature_kernel);
  stream << temperature().dispatch(nlocal) << synchronize();
}