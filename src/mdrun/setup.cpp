#include "atom.h"
#include "integrate.h"
#include "neighbor.h"
#include "thermo.h"
#include "types.h"
#include <cmath>
#include <cstdio>

#include <cstdio>
#include <cstring>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

double random(int *);

#define NSECTIONS 3
#define MAXLINE 255
char line[MAXLINE];
char keyword[MAXLINE];
FILE *fp;

void read_lammps_parse_keyword(int first) {
  int eof = 0;
  char buffer[MAXLINE];

  // proc 0 reads upto non-blank line plus 1 following line
  // eof is set to 1 if any read hits end-of-file

  if (!first) {
    if (fgets(line, MAXLINE, fp) == NULL)
      eof = 1;
  }

  while (eof == 0 && strspn(line, " \t\n\r") == strlen(line)) {
    if (fgets(line, MAXLINE, fp) == NULL)
      eof = 1;
  }

  if (fgets(buffer, MAXLINE, fp) == NULL)
    eof = 1;

  // if eof, set keyword empty and return

  if (eof) {
    keyword[0] = '\0';
    return;
  }

  // bcast keyword line to all procs

  // copy non-whitespace portion of line into keyword

  int start = strspn(line, " \t\n\r");
  int stop = strlen(line) - 1;

  while (line[stop] == ' ' || line[stop] == '\t' || line[stop] == '\n' ||
         line[stop] == '\r')
    stop--;

  line[stop + 1] = '\0';
  strcpy(keyword, &line[start]);
}

void read_lammps_header(Atom &atom) {
  int n;
  char *ptr;

  // customize for new sections

  const char *section_keywords[NSECTIONS] = {"Atoms", "Velocities", "Masses"};

  // skip 1st line of file

  char *eof = fgets(line, MAXLINE, fp);

  // customize for new header lines
  int ntypes = 0;

  while (1) {

    if (fgets(line, MAXLINE, fp) == NULL)
      n = 0;
    else
      n = strlen(line) + 1;

    if (n == 0) {
      line[0] = '\0';
      return;
    }

    // trim anything from '#' onward
    // if line is blank, continue

    double xlo, xhi, ylo, yhi, zlo, zhi;

    if ((ptr = strchr(line, '#')))
      *ptr = '\0';

    if (strspn(line, " \t\n\r") == strlen(line))
      continue;

    // search line for header keyword and set corresponding variable

    if (strstr(line, "atoms"))
      sscanf(line, "%i", &atom.natoms);
    else if (strstr(line, "atom types"))
      sscanf(line, "%i", &ntypes);

    // check for these first
    // otherwise "triangles" will be matched as "angles"

    else if (strstr(line, "xlo xhi")) {
      sscanf(line, "%lg %lg", &xlo, &xhi);
      atom.box.xlen = xhi - xlo;
    } else if (strstr(line, "ylo yhi")) {
      sscanf(line, "%lg %lg", &ylo, &yhi);
      atom.box.ylen = yhi - ylo;
    } else if (strstr(line, "zlo zhi")) {
      sscanf(line, "%lg %lg", &zlo, &zhi);
      atom.box.zlen = zhi - zlo;
    } else
      break;
  }

  // error check on total system size

  // check that exiting string is a valid section keyword

  read_lammps_parse_keyword(1);

  for (n = 0; n < NSECTIONS; n++)
    if (strcmp(keyword, section_keywords[n]) == 0)
      break;

  if (n == NSECTIONS) {
    char str[128];
    sprintf(str, "Unknown identifier in data file: %s", keyword);
  }

  // error check on consistency of header values
}

int read_lammps_data(Stream &stream, Device &device, Atom &atom,
                     Neighbor &neighbor, Integrate &integrate, Thermo &thermo,
                     char *file, int units) {
  fp = fopen(file, "r");

  if (fp == NULL) {
    char str[128];
    sprintf(str, "Cannot open file %s", file);
  }

  read_lammps_header(atom);
  int natoms = atom.natoms;

  // if (neighbor.nbinx < 0) {
  //   float volume = atom.box.xlen * atom.box.ylen * atom.box.zlen;
  //   float rho = 1.0 * atom.natoms / volume;
  //   float neigh_bin_size = pow(rho * 16, float(1.0 / 3.0));
  //   neighbor.nbinx = atom.box.xlen / neigh_bin_size;
  //   neighbor.nbiny = atom.box.ylen / neigh_bin_size;
  //   neighbor.nbinz = atom.box.zlen / neigh_bin_size;
  // }

  // if (neighbor.nbinx == 0)
  //   neighbor.nbinx = 1;

  // if (neighbor.nbiny == 0)
  //   neighbor.nbiny = 1;

  // if (neighbor.nbinz == 0)
  //   neighbor.nbinz = 1;

  neighbor.setup(stream, device, atom);

  integrate.setup();
  thermo.setup(stream, device,
               atom.box.xlen * atom.box.ylen * atom.box.zlen / natoms,
               integrate, atom, units);

  std::vector<float3> x_vec(natoms);
  std::vector<float3> v_vec(natoms);

  int atomflag = 0;
  int tmp;
  int i;
  int nread = 0;

  while (strlen(keyword)) {
    if (strcmp(keyword, "Atoms") == 0) {

      atom.nlocal = 0;

      int type;
      double xx, xy, xz;

      while (nread < natoms) {
        fgets(line, MAXLINE, fp);
        sscanf(line, "%i %i %lg %lg %lg", &i, &type, &xx, &xy, &xz);
        i--;
        x_vec.push_back(make_float3(xx, xy, xz));
        nread++;
      }
      atomflag = 1;
    } else if (strcmp(keyword, "Velocities") == 0) {
      if (atomflag == 0)
        printf("Must read Atoms before Velocities\n");

      int i;

      int nread = 0;
      int natoms = atom.natoms;

      double x, y, z;

      while (nread < natoms) {
        fgets(line, MAXLINE, fp);
        sscanf(line, "%i %lg %lg %lg", &i, &x, &y, &z);
        i--;
        v_vec.push_back(make_float3(x, y, z));
        nread++;
      }
    } else if (strcmp(keyword, "Masses") == 0) {
      fgets(line, MAXLINE, fp);

#if PRECISION == 1
      sscanf(line, "%i %g", &tmp, &atom.mass);
#else
      sscanf(line, "%i %g", &tmp, &atom.mass);
#endif
    }

    read_lammps_parse_keyword(0);
  }

  for (int i = 0; i < natoms; i++) {
    if (x_vec[i][0] >= atom.box.xlo && x_vec[i][0] < atom.box.xhi &&
        x_vec[i][1] >= atom.box.ylo && x_vec[i][1] < atom.box.yhi &&
        x_vec[i][2] >= atom.box.zlo && x_vec[i][2] < atom.box.zhi)
      atom.addatom(x_vec[i][0], x_vec[i][1], x_vec[i][2], v_vec[i][0],
                   v_vec[i][1], v_vec[i][2]);
  }

  return 0;
}

/* create simulation box */

void create_box(Atom &atom, int nx, int ny, int nz, double rho) {
  double lattice = pow((4.0 / rho), (1.0 / 3.0));
  atom.box.xlen = nx * lattice;
  atom.box.ylen = ny * lattice;
  atom.box.zlen = nz * lattice;
  atom.box.xlo = 0;
  atom.box.xhi = atom.box.xlen;
  atom.box.ylo = 0;
  atom.box.yhi = atom.box.ylen;
  atom.box.zlo = 0;
  atom.box.zhi = atom.box.zlen;
}

/* initialize atoms on fcc lattice in parallel fashion */

int create_atoms(Atom &atom, int nx, int ny, int nz, double rho) {
  /* total # of atoms */

  atom.natoms = 4 * nx * ny * nz;
  atom.nlocal = 0;

  /* determine loop bounds of lattice subsection that overlaps my sub-box
     insure loop bounds do not exceed nx,ny,nz */

  double alat = pow((4.0 / rho), (1.0 / 3.0));
  int ilo = static_cast<int>(atom.box.xlo / (0.5 * alat) - 1);
  int ihi = static_cast<int>(atom.box.xhi / (0.5 * alat) + 1);
  int jlo = static_cast<int>(atom.box.ylo / (0.5 * alat) - 1);
  int jhi = static_cast<int>(atom.box.yhi / (0.5 * alat) + 1);
  int klo = static_cast<int>(atom.box.zlo / (0.5 * alat) - 1);
  int khi = static_cast<int>(atom.box.zhi / (0.5 * alat) + 1);

  ilo = MAX(ilo, 0);
  ihi = MIN(ihi, 2 * nx - 1);
  jlo = MAX(jlo, 0);
  jhi = MIN(jhi, 2 * ny - 1);
  klo = MAX(klo, 0);
  khi = MIN(khi, 2 * nz - 1);

  double xtmp, ytmp, ztmp, vx, vy, vz;
  int i, j, k, m, n;
  int sx = 0;
  int sy = 0;
  int sz = 0;
  int ox = 0;
  int oy = 0;
  int oz = 0;
  int subboxdim = 8;

  int iflag = 0;

  while (oz * subboxdim <= khi) {
    k = oz * subboxdim + sz;
    j = oy * subboxdim + sy;
    i = ox * subboxdim + sx;

    if (iflag)
      continue;

    if (((i + j + k) % 2 == 0) && (i >= ilo) && (i <= ihi) && (j >= jlo) &&
        (j <= jhi) && (k >= klo) && (k <= khi)) {

      xtmp = 0.5 * alat * i;
      ytmp = 0.5 * alat * j;
      ztmp = 0.5 * alat * k;

      if (xtmp >= atom.box.xlo && xtmp < atom.box.xhi && ytmp >= atom.box.ylo &&
          ytmp < atom.box.yhi && ztmp >= atom.box.zlo && ztmp < atom.box.zhi) {
        n = k * (2 * ny) * (2 * nx) + j * (2 * nx) + i + 1;

        for (m = 0; m < 5; m++)
          random(&n);

        vx = random(&n);

        for (m = 0; m < 5; m++)
          random(&n);

        vy = random(&n);

        for (m = 0; m < 5; m++)
          random(&n);

        vz = random(&n);

        atom.addatom(xtmp, ytmp, ztmp, vx, vy, vz);
      }
    }

    sx++;

    if (sx == subboxdim) {
      sx = 0;
      sy++;
    }

    if (sy == subboxdim) {
      sy = 0;
      sz++;
    }

    if (sz == subboxdim) {
      sz = 0;
      ox++;
    }

    if (ox * subboxdim > ihi) {
      ox = 0;
      oy++;
    }

    if (oy * subboxdim > jhi) {
      oy = 0;
      oz++;
    }
  }

  return 0;
}

/* adjust initial velocities to give desired temperature */

void create_velocity(Stream &stream, Device &device, float t_request,
                     Atom &atom, Thermo &thermo) {

  atom.construct_buf(stream, device);

  int i;

  Buffer<float> vtot = device.create_buffer<float>(3); // x, y,z

  Kernel1D v_add_kernel = [&]() noexcept {
    auto i = dispatch_x();
    vtot->atomic(0).fetch_add(atom.v->read(i)[0]);
    vtot->atomic(1).fetch_add(atom.v->read(i)[1]);
    vtot->atomic(2).fetch_add(atom.v->read(i)[2]);
  };
  Kernel1D v_average_kernel = [&]() noexcept {
    auto i = dispatch_x();
    Int natoms_ = atom.natoms;
    Float3 v_mean;
    v_mean[0] = vtot->read(0) / natoms_.cast<float>();
    v_mean[1] = vtot->read(1) / natoms_.cast<float>();
    v_mean[2] = vtot->read(2) / natoms_.cast<float>();
    Float v_origin = atom.v->read(i);
    atom.v->write(i, v_origin - v_mean);
  };
  auto v_add = device.compile(v_add_kernel);
  auto v_average = device.compile(v_average_kernel);
  stream << v_add().dispatch(atom.nlocal) << synchronize();
  stream << v_average().dispatch(atom.nlocal) << synchronize();

  thermo.temperature(stream, device, atom); // compute t_act

  Kernel1D v_scale_kernel = [&]() noexcept {
    auto i = dispatch_x();
    Float factor = sqrt(t_request / (thermo.t_act->read(0) * thermo.t_scale));
    Float3 v_scaled = atom.v->read(i) * factor;
    atom.v->write(i, v_scaled);
  };
  auto v_scale = device.compile(v_scale_kernel);
  stream << v_scale().dispatch(atom.nlocal) << synchronize();
}

/* Park/Miller RNG w/out MASKING, so as to be like f90s version */

#define IA 16807
#define IM 2147483647
#define AM (1.0 / IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

double random(int *idum) {
  int k;
  double ans;

  k = (*idum) / IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;

  if (*idum < 0)
    *idum += IM;

  ans = AM * (*idum);
  return ans;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK
