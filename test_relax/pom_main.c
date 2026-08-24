#define _GNU_SOURCE
/* pom_main.c
 *
 * command line driver for the Plate_over_Maxwell C port
 *
 * computes postseismic surface displacements or velocities on a
 * regular grid or at stations read from a file, and writes a plain
 * table (x y ue un uz) suitable for GMT
 *
 * usage:
 *   pom layer  [options]   plate over viscoelastic layer over halfspace
 *   pom maxwell [options]  plate over Maxwell halfspace (point source)
 *
 * options (defaults reproduce run_postseismic_example.m):
 *   -m L W depth dip strike east north ss ds ten   fault patch (10 values)
 *   -grid xmin xmax nx ymin ymax ny                output grid (km)
 *   -sta file                                      stations x y, one per line
 *   -t time                                        years since earthquake
 *   -vel                                           velocities, not displacements
 *   -H1 v -H2 v -nu v -tR1 v -tR2 v                layer model parameters
 *   -H v -mulam v -tR v                            maxwell model parameters
 *   -nl n -nw n                                    quadrature override (layer)
 *   -o file                                        output file (default stdout)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pom.h"

static void usage(void)
{
  fprintf(stderr,
	  "usage: pom layer|maxwell [-m L W depth dip strike east north ss ds ten]\n"
	  "  [-grid xmin xmax nx ymin ymax ny] [-sta file] [-t time] [-vel]\n"
	  "  [-H1 v] [-H2 v] [-nu v] [-tR1 v] [-tR2 v]\n"
	  "  [-H v] [-mulam v] [-tR v] [-nl n] [-nw n] [-o file]\n");
  exit(1);
}

int main(int argc, char **argv)
{
  /* defaults: the example in run_postseismic_example.m */
  double m[10] = { 200, 20, 20, 90, 0, 0, 0, 1, 0, 0 };
  double xmin = -300, xmax = 300, ymin = -300, ymax = 300;
  int nx = 24, ny = 24;
  double H1 = 20, H2 = 40, nu = 0.25, tR1 = 25, tR2 = 25;
  double H = 20, mulam = 1.0, tR = 25;
  double zr = 0.0;
  double *szz = NULL, *sxz = NULL, *syz = NULL;
  double t = 1.0;
  int mode = POM_DISP, nl = 0, nw = 0;
  const char *stafile = NULL, *ofile = NULL;
  int use_maxwell = 0;
  int i, j, n, rc;
  double *xs, *ys, *ue, *un, *uz;
  FILE *out;

  if (argc < 2)
    usage();
  if (strcmp(argv[1], "maxwell") == 0)
    use_maxwell = 1;
  else if (strcmp(argv[1], "layer") != 0)
    usage();
  for (i = 2; i < argc; i++) {
    if (strcmp(argv[i], "-m") == 0 && i + 10 < argc) {
      for (j = 0; j < 10; j++)
	m[j] = atof(argv[++i]);
    } else if (strcmp(argv[i], "-grid") == 0 && i + 6 < argc) {
      xmin = atof(argv[++i]);
      xmax = atof(argv[++i]);
      nx = atoi(argv[++i]);
      ymin = atof(argv[++i]);
      ymax = atof(argv[++i]);
      ny = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-sta") == 0 && i + 1 < argc) {
      stafile = argv[++i];
    } else if (strcmp(argv[i], "-t") == 0 && i + 1 < argc) {
      t = atof(argv[++i]);
    } else if (strcmp(argv[i], "-vel") == 0) {
      mode = POM_VEL;
    } else if (strcmp(argv[i], "-zr") == 0 && i + 1 < argc) {
      zr = atof(argv[++i]);
    } else if (strcmp(argv[i], "-H1") == 0 && i + 1 < argc) {
      H1 = atof(argv[++i]);
    } else if (strcmp(argv[i], "-H2") == 0 && i + 1 < argc) {
      H2 = atof(argv[++i]);
    } else if (strcmp(argv[i], "-nu") == 0 && i + 1 < argc) {
      nu = atof(argv[++i]);
    } else if (strcmp(argv[i], "-tR1") == 0 && i + 1 < argc) {
      tR1 = atof(argv[++i]);
    } else if (strcmp(argv[i], "-tR2") == 0 && i + 1 < argc) {
      tR2 = atof(argv[++i]);
    } else if (strcmp(argv[i], "-H") == 0 && i + 1 < argc) {
      H = atof(argv[++i]);
    } else if (strcmp(argv[i], "-mulam") == 0 && i + 1 < argc) {
      mulam = atof(argv[++i]);
    } else if (strcmp(argv[i], "-tR") == 0 && i + 1 < argc) {
      tR = atof(argv[++i]);
    } else if (strcmp(argv[i], "-nl") == 0 && i + 1 < argc) {
      nl = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-nw") == 0 && i + 1 < argc) {
      nw = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
      ofile = argv[++i];
    } else {
      usage();
    }
  }

  if (stafile) {
    double x, y;
    FILE *f;

    f = fopen(stafile, "r");
    if (!f) {
      fprintf(stderr, "pom: cannot open %s\n", stafile);
      return 1;
    }
    n = 0;
    while (fscanf(f, "%lf %lf", &x, &y) == 2)
      n++;
    rewind(f);
    xs = malloc(sizeof(double) * n);
    ys = malloc(sizeof(double) * n);
    for (i = 0; i < n; i++)
      if (fscanf(f, "%lf %lf", &xs[i], &ys[i]) != 2) {
	fprintf(stderr, "pom: read error in %s\n", stafile);
	return 1;
      }
    fclose(f);
  } else {
    n = nx * ny;
    xs = malloc(sizeof(double) * n);
    ys = malloc(sizeof(double) * n);
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++) {
	xs[j * nx + i] = (nx > 1) ? xmin + (xmax - xmin) * i / (nx - 1) : xmin;
	ys[j * nx + i] = (ny > 1) ? ymin + (ymax - ymin) * j / (ny - 1) : ymin;
      }
  }
  ue = malloc(sizeof(double) * n);
  un = malloc(sizeof(double) * n);
  uz = malloc(sizeof(double) * n);
  if (!xs || !ys || !ue || !un || !uz) {
    fprintf(stderr, "pom: out of memory\n");
    return 1;
  }

  if (use_maxwell) {
    if (nl > 0 || nw > 0)
      fprintf(stderr, "pom: warning: -nl/-nw have no effect for maxwell "
	      "(single point source by construction)\n");
    if (mode == POM_VEL)
      fprintf(stderr, "pom: -vel not available for maxwell; "
	      "computing displacements\n");
    if (zr > 0.0) {
      szz = malloc(sizeof(double) * n);
      sxz = malloc(sizeof(double) * n);
      syz = malloc(sizeof(double) * n);
      rc = pom_maxwell_z(m, xs, ys, n, H, mulam, &t, 1, tR, zr,
			 ue, un, uz, szz, sxz, syz);
    } else {
      rc = pom_maxwell(m, xs, ys, n, H, mulam, &t, 1, tR, ue, un, uz);
    }
  } else {
    rc = pom_layer(m, xs, ys, n, H1, H2, nu, t, tR1, tR2, mode,
		   nl, nw, ue, un, uz);
  }
  if (rc) {
    fprintf(stderr, "pom: computation failed\n");
    return 1;
  }

  out = ofile ? fopen(ofile, "w") : stdout;
  if (!out) {
    fprintf(stderr, "pom: cannot open %s\n", ofile);
    return 1;
  }
  if (szz) {
    fprintf(out, "# x_km y_km u_east_m u_north_m u_up_m szz sxz syz"
	    "  (t = %g yr, zr = %g km, transient, stress in G units)\n",
	    t, zr);
    for (i = 0; i < n; i++)
      fprintf(out, "%g %g %19.12e %19.12e %19.12e %19.12e %19.12e %19.12e\n",
	      xs[i], ys[i], ue[i], un[i], uz[i], szz[i], sxz[i], syz[i]);
  } else {
    fprintf(out, "# x_km y_km u_east_m u_north_m u_up_m  (t = %g yr, %s)\n",
	    t, (mode == POM_VEL) ? "velocity m/yr" : "displacement m");
    for (i = 0; i < n; i++)
      fprintf(out, "%g %g %19.12e %19.12e %19.12e\n",
	      xs[i], ys[i], ue[i], un[i], uz[i]);
  }
  if (ofile)
    fclose(out);
  free(xs);
  free(ys);
  free(ue);
  free(un);
  free(uz);
  return 0;
}
