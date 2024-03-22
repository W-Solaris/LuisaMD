#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "atom.h"
#include "force.h"
#include "integrate.h"
#include "lj.h"
#include "neighbor.h"
#include "thermo.h"
#include "types.h"

#define MAXLINE 256

int input(In &in, const char *filename) {
  FILE *fp;
  int flag;
  char line[MAXLINE];

  fp = fopen(filename, "r");

  if (fp == NULL)
    flag = 0;
  else
    flag = 1;

  if (flag == 0) {
    printf("ERROR: Cannot open %s\n", filename);
    return 1;
  }

#if PRECISION == 1
  fgets(line, MAXLINE, fp);
  fgets(line, MAXLINE, fp);
  fgets(line, MAXLINE, fp);

  if (strcmp(strtok(line, " \t\n"), "lj") == 0)
    in.units = 0;
  else if (strcmp(strtok(line, " \t\n"), "metal") == 0)
    in.units = 1;
  else {
    printf("Unknown units option in file at line 3 ('%s'). Expecting either "
           "'lj' or 'metal'.\n",
           line);
    exit(0);
  }

  fgets(line, MAXLINE, fp);

  if (strcmp(strtok(line, " \t\n"), "none") == 0)
    in.datafile = NULL;
  else {
    in.datafile = new char[1000];
    char *ptr = strtok(line, " \t");

    if (ptr == NULL)
      ptr = line;

    strcpy(in.datafile, ptr);
  }

  fgets(line, MAXLINE, fp);

  if (strcmp(strtok(line, " \t\n"), "lj") == 0)
    in.forcetype = FORCELJ;
  else if (strcmp(strtok(line, " \t\n"), "eam") == 0)
    in.forcetype = FORCEEAM;
  else {
    printf("Unknown forcetype option in file at line 5 ('%s'). Expecting "
           "either 'lj' or 'eam'.\n",
           line);
    exit(0);
  }

  fgets(line, MAXLINE, fp);
  sscanf(line, "%e %e", &in.epsilon, &in.sigma);
  fgets(line, MAXLINE, fp);
  sscanf(line, "%d %d %d", &in.nx, &in.ny, &in.nz);
  fgets(line, MAXLINE, fp);
  sscanf(line, "%d", &in.ntimes);
  fgets(line, MAXLINE, fp);
  sscanf(line, "%e", &in.dt);
  fgets(line, MAXLINE, fp);
  sscanf(line, "%e", &in.t_request);
  fgets(line, MAXLINE, fp);
  sscanf(line, "%e", &in.rho);
  fgets(line, MAXLINE, fp);
  sscanf(line, "%d", &in.neigh_every);
  fgets(line, MAXLINE, fp);
  sscanf(line, "%e %e", &in.force_cut, &in.neigh_cut);
  fgets(line, MAXLINE, fp);
  sscanf(line, "%d", &in.thermo_nstat);
  fclose(fp);
#else
#if PRECISION == 2
  fgets(line, MAXLINE, fp);
  fgets(line, MAXLINE, fp);
  fgets(line, MAXLINE, fp);

  if (strcmp(strtok(line, " \t\n"), "lj") == 0)
    in.units = 0;
  else if (strcmp(line, "metal") == 0)
    in.units = 1;
  else {
    printf("Unknown units option in file at line 3 ('%s'). Expecting either "
           "'lj' or 'metal'.\n",
           line);
    exit(0);
  }

  fgets(line, MAXLINE, fp);

  if (strcmp(strtok(line, " \t\n"), "none") == 0)
    in.datafile = NULL;
  else {
    in.datafile = new char[1000];
    char *ptr = strtok(line, " \t");

    if (ptr == NULL)
      ptr = line;

    strcpy(in.datafile, ptr);
  }

  fgets(line, MAXLINE, fp);

  if (strcmp(strtok(line, " \t\n"), "lj") == 0)
    in.forcetype = FORCELJ;
  else if (strcmp(line, "eam") == 0)
    in.forcetype = FORCEEAM;
  else {
    printf("Unknown forcetype option in file at line 5 ('%s'). Expecting "
           "either 'lj' or 'eam'.\n",
           line);
    exit(0);
  }

  fgets(line, MAXLINE, fp);
  sscanf(line, "%le %le", &in.epsilon, &in.sigma);
  fgets(line, MAXLINE, fp);
  sscanf(line, "%d %d %d", &in.nx, &in.ny, &in.nz);
  fgets(line, MAXLINE, fp);
  sscanf(line, "%d", &in.ntimes);
  fgets(line, MAXLINE, fp);
  sscanf(line, "%le", &in.dt);
  fgets(line, MAXLINE, fp);
  sscanf(line, "%le", &in.t_request);
  fgets(line, MAXLINE, fp);
  sscanf(line, "%le", &in.rho);
  fgets(line, MAXLINE, fp);
  sscanf(line, "%d", &in.neigh_every);
  fgets(line, MAXLINE, fp);
  sscanf(line, "%le %le", &in.force_cut, &in.neigh_cut);
  fgets(line, MAXLINE, fp);
  sscanf(line, "%d", &in.thermo_nstat);
  fclose(fp);
#else
  printf("Invalid MMD_float size specified: crash imminent.\n");
#endif
#endif

  in.neigh_cut += in.force_cut;

  return 0;
}
