#include <time.h>

#include "atom.h"
#include "force.h"
#include "integrate.h"
#include "lj.h"
#include "luisa/runtime/stream.h"
#include "math.h"
#include "neighbor.h"
#include "stdio.h"
#include "stdlib.h"
#include "thermo.h"
#include "variant.h"

void stats(int, double *, double *, double *, double *, int, int *);

void output(Stream &stream, Device &device, In &in, Atom &atom, Force *force,
            Neighbor &neighbor, Thermo &thermo, Integrate &integrate,
            int screen_yaml) {
  int i, n;
  int histo[10];
  double tmp, ave, max, min, total;
  FILE *fp;

  /* enforce PBC, then check for lost atoms */

  atom.pbc(stream);

  int natoms = atom.nlocal;
  int nlost = 0;

  /* long-range energy and pressure corrections Whats this???*/

  double engcorr = 8.0 * 3.1415926 * in.rho *
                   (1.0 / (9.0 * pow(force->cutforce, float(9.0))) -
                    1.0 / (3.0 * pow(force->cutforce, float(3.0))));
  double prscorr = 8.0 * 3.1415926 * in.rho * in.rho *
                   (4.0 / (9.0 * pow(force->cutforce, float(9.0))) -
                    2.0 / (3.0 * pow(force->cutforce, float(3.0))));

  /* thermo output */

  double conserve;

  time_t general_time = time(NULL);
  struct tm local_time = *localtime(&general_time);
  char filename[256];

  sprintf(filename, "miniMD-%4d-%02d-%02d-%02d-%02d-%02d.yaml",
          local_time.tm_year + 1900, local_time.tm_mon + 1, local_time.tm_mday,
          local_time.tm_hour, local_time.tm_min, local_time.tm_sec);

  fp = fopen(filename, "w");

  if (screen_yaml) {
    fprintf(stdout, "run_configuration: \n");
    fprintf(stdout, "  variant: " VARIANT_STRING "\n");
    fprintf(stdout, "  datafile: %s\n", in.datafile ? in.datafile : "None");
    fprintf(stdout, "  units: %s\n", in.units == 0 ? "LJ" : "METAL");
    fprintf(stdout, "  atoms: %i\n", atom.natoms);
    fprintf(stdout, "  system_size: %2.2lf %2.2lf %2.2lf\n", atom.box.xlen,
            atom.box.ylen, atom.box.zlen);
    fprintf(stdout, "  unit_cells: %i %i %i\n", in.nx, in.ny, in.nz);
    fprintf(stdout, "  density: %lf\n", in.rho);
    fprintf(stdout, "  force_type: %s\n",
            in.forcetype == FORCELJ ? "LJ" : "EAM");
    fprintf(stdout, "  force_cutoff: %lf\n", force->cutforce);
    fprintf(stdout, "  force_params: %2.2lf %2.2lf\n", force->epsilon_scalar,
            force->sigma_scalar);
    fprintf(stdout, "  neighbor_cutoff: %lf\n", neighbor.cutneigh);
    fprintf(stdout, "  neighbor_type: %i\n", neighbor.halfneigh);
    fprintf(stdout, "  neighbor_team_build: %i\n", neighbor.team_neigh_build);
    fprintf(stdout, "  neighbor_bins: %i %i %i\n", neighbor.nbinx,
            neighbor.nbiny, neighbor.nbinz);
    fprintf(stdout, "  neighbor_frequency: %i\n", neighbor.every);
    fprintf(stdout, "  sort_frequency: %i\n", integrate.sort_every);
    fprintf(stdout, "  timestep_size: %lf\n", integrate.dt);
    fprintf(stdout, "  thermo_frequency: %i\n", thermo.nstat);
    fprintf(stdout, "  use_intrinsics: %i\n", force->use_sse);
    fprintf(stdout, "  float_size: %i\n\n", (int)sizeof(float));
  }

  fprintf(fp, "run_configuration: \n");
  fprintf(fp, "  variant: " VARIANT_STRING "\n");
  fprintf(fp, "  datafile: %s\n", in.datafile ? in.datafile : "None");
  fprintf(fp, "  units: %s\n", in.units == 0 ? "LJ" : "METAL");
  fprintf(fp, "  atoms: %i\n", atom.natoms);
  fprintf(fp, "  system_size: %2.2lf %2.2lf %2.2lf\n", atom.box.xlen,
          atom.box.ylen, atom.box.zlen);
  fprintf(fp, "  unit_cells: %i %i %i\n", in.nx, in.ny, in.nz);
  fprintf(fp, "  density: %lf\n", in.rho);
  fprintf(fp, "  force_type: %s\n", in.forcetype == FORCELJ ? "LJ" : "EAM");
  fprintf(fp, "  force_cutoff: %lf\n", force->cutforce);
  fprintf(fp, "  force_params: %2.2lf %2.2lf\n", force->epsilon_scalar,
          force->sigma_scalar);
  fprintf(fp, "  neighbor_cutoff: %lf\n", neighbor.cutneigh);
  fprintf(fp, "  neighbor_type: %i\n", neighbor.halfneigh);
  fprintf(fp, "  neighbor_team_build: %i\n", neighbor.team_neigh_build);
  fprintf(fp, "  neighbor_bins: %i %i %i\n", neighbor.nbinx, neighbor.nbiny,
          neighbor.nbinz);
  fprintf(fp, "  neighbor_frequency: %i\n", neighbor.every);
  fprintf(fp, "  sort_frequency: %i\n", integrate.sort_every);
  fprintf(fp, "  timestep_size: %lf\n", integrate.dt);
  fprintf(fp, "  thermo_frequency: %i\n", thermo.nstat);
  fprintf(fp, "  ghost_newton: %i\n", neighbor.ghost_newton);
  fprintf(fp, "  use_intrinsics: %i\n", force->use_sse);
  fprintf(fp, "  float_size: %i\n\n", (int)sizeof(float));

  if (screen_yaml) fprintf(stdout, "\n\nthermodynamic_output:\n");

  fprintf(fp, "\n\nthermodynamic_output:\n");

  // copy info from device
  std::vector<int> h_mstat(1);
  stream << thermo.mstat.copy_to(h_mstat.data());
  std::vector<int> h_steparr(h_mstat[0]);
  stream << thermo.steparr.copy_to(h_steparr.data());
  std::vector<float> h_tmparr(h_mstat[0]);
  stream << thermo.tmparr.copy_to(h_tmparr.data());
  std::vector<float> h_engarr(h_mstat[0]);
  stream << thermo.engarr.copy_to(h_engarr.data());
  std::vector<float> h_prsarr(h_mstat[0]);
  stream << thermo.prsarr.copy_to(h_prsarr.data());

  for (i = 0; i < h_mstat[0]; i++) {
    conserve =
        (1.5 * h_tmparr[i] + h_tmparr[i]) / (1.5 * h_tmparr[0] + h_tmparr[0]);

    if (screen_yaml) {
      fprintf(stdout, "  timestep: %d \n", h_steparr[i]);
      fprintf(stdout, "      T*:           %15.10g \n", h_tmparr[i]);
      // fprintf(stdout,"      U*:           %15.10g \n",
      // thermo.engarr[i]+engcorr); fprintf(stdout,"      P*: %15.10g \n",
      // thermo.prsarr[i]+prscorr);
      fprintf(stdout, "      U*:           %15.10g \n", h_engarr[i]);
      fprintf(stdout, "      P*:           %15.10g \n", h_prsarr[i]);
      fprintf(stdout, "      Conservation: %15.10g \n", conserve);
    }

    fprintf(fp, "  timestep: %d \n", h_steparr[i]);
    fprintf(fp, "      T*:           %15.10g \n", h_tmparr[i]);
    // fprintf(fp    ,"      U*:           %15.10g \n",
    // thermo.engarr[i]+engcorr); fprintf(fp    ,"      P*:           %15.10g
    // \n", thermo.prsarr[i]+prscorr);
    fprintf(fp, "      U*:           %15.10g \n", h_engarr[i]);
    fprintf(fp, "      P*:           %15.10g \n", h_prsarr[i]);
    fprintf(fp, "      Conservation: %15.10g \n", conserve);
  }

  /* performance output */

  fprintf(stdout, "\n\n");
  fprintf(fp, "\n\n");

  fprintf(stdout, "\n");
  fprintf(fp, "\n");

  tmp = atom.nlocal;
  stats(1, &tmp, &ave, &max, &min, 10, histo);

  if (screen_yaml)
    fprintf(stdout, "# Nlocal:     %g ave %g max %g min\n", ave, max, min);

  if (screen_yaml) fprintf(stdout, "# Histogram:");

  if (screen_yaml)
    for (i = 0; i < 10; i++) fprintf(stdout, " %d", histo[i]);

  if (screen_yaml) fprintf(stdout, "\n");

  fprintf(fp, "# Nlocal:     %g ave %g max %g min\n", ave, max, min);
  fprintf(fp, "# Histogram:");

  for (i = 0; i < 10; i++) fprintf(fp, " %d", histo[i]);

  fprintf(fp, "\n");

  // tmp = atom.nghost;
  // stats(1, &tmp, &ave, &max, &min, 10, histo);

  // if (screen_yaml)
  //   fprintf(stdout, "# Nghost:     %g ave %g max %g min\n", ave, max, min);

  // if (screen_yaml)
  //   fprintf(stdout, "# Histogram:");

  // if (screen_yaml)
  //   for (i = 0; i < 10; i++)
  //     fprintf(stdout, " %d", histo[i]);

  // if (screen_yaml)
  //   fprintf(stdout, "\n");

  // fprintf(fp, "# Nghost:     %g ave %g max %g min\n", ave, max, min);
  // fprintf(fp, "# Histogram:");

  // for (i = 0; i < 10; i++)
  //   fprintf(fp, " %d", histo[i]);

  // fprintf(fp, "\n");

  n = 0;

  std::vector<int> h_numneigh(neighbor.nmax);
  stream << neighbor.numneigh.copy_to(h_numneigh.data()) << synchronize();
  for (i = 0; i < atom.nlocal; i++) n += h_numneigh[i];

  total = n;
  stats(1, &total, &ave, &max, &min, 10, histo);

  if (screen_yaml)
    fprintf(stdout, "# Neighs:     %g ave %g max %g min\n", ave, max, min);

  if (screen_yaml) fprintf(stdout, "# Histogram:");

  if (screen_yaml)
    for (i = 0; i < 10; i++) fprintf(stdout, " %d", histo[i]);

  if (screen_yaml) fprintf(stdout, "\n");

  fprintf(fp, "# Neighs:     %g ave %g max %g min\n", ave, max, min);
  fprintf(fp, "# Histogram:");

  for (i = 0; i < 10; i++) fprintf(fp, " %d", histo[i]);

  fprintf(fp, "\n");

  if (screen_yaml) fprintf(stdout, "# Total # of neighbors = %g\n", total);

  fprintf(fp, "# Total # of neighbors = %g\n", total);

  if (screen_yaml) fprintf(stdout, "\n");

  fprintf(fp, "\n");

  fclose(fp);
}

void stats(int n, double *data, double *pave, double *pmax, double *pmin,
           int nhisto, int *histo) {
  int i, m;
  int *histotmp;

  double min = 1.0e20;
  double max = -1.0e20;
  double ave = 0.0;

  for (i = 0; i < n; i++) {
    ave += data[i];

    if (data[i] < min) min = data[i];

    if (data[i] > max) max = data[i];
  }

  ave /= n;

  for (i = 0; i < nhisto; i++) histo[i] = 0;

  double del = max - min;

  for (i = 0; i < n; i++) {
    if (del == 0.0)
      m = 0;
    else
      m = static_cast<int>((data[i] - min) / del * nhisto);

    if (m > nhisto - 1) m = nhisto - 1;

    histo[m]++;
  }

  *pave = ave;
  *pmax = max;
  *pmin = min;
}
