/* fd3d_test: TPV3 (Day et al. 2005) with the 3-D C marcher, fixed rake,
   compared step by step with the prototype (run_tpv3.py, FIXED_RAKE=1).
   usage: fd3d_test tpv3 dx [ref_file] */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fd3d.h"
/* cells for a station at x on a grid starting at edge lo: the cell containing x,
   or the two cells sharing an edge at x */
static int near2(double x, double lo, double dx, int n, int *c)
{
  double f = (x - lo)/dx; int k = (int)floor(f + 0.5);
  if(fabs(f - k) < 1e-9){ c[0] = k - 1; c[1] = k; return 2; }
  c[0] = (int)floor(f); return 1;
}
static int compare(const char *ref, int nst, const char **nm, double *rec, int nsteps, double dt, double tol, const char *what)
{
  FILE *fr = fopen(ref, "r"); char name[32]; double r[4], emax[3] = {0, 0, 0}, smax[3] = {0, 0, 0};
  int cnt = 0, s, n, q, fail;
  if(!fr){ fprintf(stderr, "cannot open %s\n", ref); return 1; }
  while(fscanf(fr, "%31s %lf %lf %lf %lf", name, r, r + 1, r + 2, r + 3) == 5){
    for(s = 0; s < nst; s++) if(!strcmp(name, nm[s])) break;
    if(s == nst) continue;
    n = (int)floor(r[0]/dt + 0.5) - 1;
    if(n < 0 || n >= nsteps) continue;
    for(q = 0; q < 3; q++){
      double e = fabs(rec[(s*nsteps + n)*4 + 1 + q] - r[1 + q]);
      if(e > emax[q]) emax[q] = e;
      if(fabs(r[1 + q]) > smax[q]) smax[q] = fabs(r[1 + q]);
    }
    cnt++;
  }
  fclose(fr);
  printf("fd3d_test %s vs %s: %d samples; max |diff| slip %.2e (of %.3f), slip rate %.2e (of %.3f), stress %.2e (of %.3e)\n",
         what, ref, cnt, emax[0], smax[0], emax[1], smax[1], emax[2], smax[2]);
  fail = (cnt == 0 || emax[0] > tol*smax[0] || emax[1] > tol*smax[1] || emax[2] > tol*smax[2]);
  printf("fd3d_test %s: %s\n", what, fail ? "FAIL" : "PASS");
  return fail;
}

/* SCEC TPV101 (rate and state, aging law, whole space), fixed rake, as in
   fully_dynamic/run_tpv101.py: fault x in [-21, 21] km, y in [-6, 21] km */
static int run_tpv101(int argc, char **argv)
{
  double cL = 6000, cT = 3464, rho = 2670, dx = atof(argv[2]), tend = argc > 4 ? atof(argv[4]) : 12.0, dt;
  double xlo = -21e3, ylo = -6e3, xlen = 42e3, ylen = 27e3, W = 15e3, w = 3e3, *load, *Fr, R = 3e3, T = 1.0, *rt;
  int nx, nz, i, j, n, nsteps, s, k, N;
  const char *nm[3] = {"st000dp075", "st090dp075", "st000dp030"};
  double sx[3] = {0, 9e3, 0}, sy[3] = {7.5e3, 7.5e3, 3e3}, *rec;
  int sc[3][4]; double sw[3][4];
  struct fd3d_tab *t; struct fd3d *f;
  if(getenv("FD3D_SMALL")){ xlo = -19.5e3; ylo = -4.5e3; xlen = 39e3; ylen = 24e3; }   /* 1.5 km of strengthening fault beyond the transition */
  nx = (int)floor(xlen/dx + 0.5); nz = (int)floor(ylen/dx + 0.5); dt = 0.5*dx/cT; N = nx*nz;
  t = fd3d_tab_create(nx, nz, dx, dt, cL, cT, rho, 0);
  fprintf(stderr, "fd3d_test tpv101: %d x %d elements, dt %.4f s, window %d lags, FFT %d x %d\n", nx, nz, dt, t->L, t->Mx, t->Mz);
  f = fd3d_create(t); f->law = 1; f->f0 = 0.6; f->V0 = 1e-6;
  f->a = malloc(sizeof(double)*N); f->b = malloc(sizeof(double)*N); f->Lc = malloc(sizeof(double)*N); f->theta = malloc(sizeof(double)*N);
  load = malloc(sizeof(double)*N); Fr = malloc(sizeof(double)*N);
  for(i = 0; i < nx; i++) for(j = 0; j < nz; j++){
    double x = xlo + (i + 0.5)*dx, y = ylo + (j + 0.5)*dx, Bx, By, r2;
    double ax = fabs(x), ay = fabs(y - 7.5e3);
    Bx = ax <= W ? 1.0 : (ax < W + w ? 0.5*(1.0 + tanh(w/(ax - W - w) + w/(ax - W))) : 0.0);
    By = ay <= W/2 ? 1.0 : (ay < W/2 + w ? 0.5*(1.0 + tanh(w/(ay - W/2 - w) + w/(ay - W/2))) : 0.0);
    k = i*nz + j;
    f->a[k] = 0.008 + 0.008*(1.0 - Bx*By); f->b[k] = 0.012; f->Lc[k] = 0.02;
    f->sigma[k] = 120e6; f->tau0[k] = 75e6; f->V[k] = 1e-12;
    f->theta[k] = f->Lc[k]/f->V0*exp((f->a[k]*log(2.0*sinh(75e6/(f->a[k]*120e6))) - f->f0 - f->a[k]*log(1e-12/f->V0))/f->b[k]);
    r2 = x*x + (y - 7.5e3)*(y - 7.5e3);
    Fr[k] = r2 < R*R ? exp(r2/(r2 - R*R)) : 0.0;
  }
  for(s = 0; s < 3; s++){
    double fi = (sx[s] - xlo)/dx - 0.5, fj = (sy[s] - ylo)/dx - 0.5; int i0 = (int)floor(fi), j0 = (int)floor(fj); double wi = fi - i0, wj = fj - j0;
    sc[s][0] = i0*nz + j0; sw[s][0] = (1 - wi)*(1 - wj); sc[s][1] = (i0 + 1)*nz + j0; sw[s][1] = wi*(1 - wj);
    sc[s][2] = i0*nz + j0 + 1; sw[s][2] = (1 - wi)*wj; sc[s][3] = (i0 + 1)*nz + j0 + 1; sw[s][3] = wi*wj;
  }
  nsteps = (int)floor(tend/dt + 0.5); rec = malloc(sizeof(double)*3*4*nsteps);
  rt = malloc(sizeof(double)*N); for(k = 0; k < N; k++) rt[k] = -1.0;
  for(n = 0; n < nsteps; n++){
    double t1 = f->time + dt, G = (t1 > 0 && t1 < T) ? exp((t1 - T)*(t1 - T)/(t1*(t1 - 2*T))) : (t1 >= T ? 1.0 : 0.0);
    for(k = 0; k < N; k++) load[k] = 25e6*Fr[k]*G;
    fd3d_step(f, load);
    for(k = 0; k < N; k++) if(rt[k] < 0 && fabs(f->V[k]) > 1e-3) rt[k] = f->time;
    if(n % 50 == 0){ double vm = 0; for(k = 0; k < N; k++) if(fabs(f->V[k]) > vm) vm = fabs(f->V[k]); fprintf(stderr, "  t %6.3f max V %.3e\n", f->time, vm); }
    for(s = 0; s < 3; s++){
      double d = 0, v = 0, ta = 0;
      for(k = 0; k < 4; k++){ d += sw[s][k]*f->delta[sc[s][k]]; v += sw[s][k]*f->V[sc[s][k]]; ta += sw[s][k]*f->tau[sc[s][k]]; }
      rec[(s*nsteps + n)*4] = f->time; rec[(s*nsteps + n)*4 + 1] = d; rec[(s*nsteps + n)*4 + 2] = v; rec[(s*nsteps + n)*4 + 3] = ta;
    }
  }
  for(s = 0; s < 3; s++){
    double on = -1, pk = 0;
    for(n = 0; n < nsteps; n++){ double v = rec[(s*nsteps + n)*4 + 2]; if(on < 0 && v > 1e-3) on = rec[(s*nsteps + n)*4]; if(v > pk) pk = v; }
    printf("%s: rupture time %.3f s, peak slip rate %.3f m/s, slip %.3f m\n", nm[s], on, pk, rec[(s*nsteps + nsteps - 1)*4 + 1]);
  }
  { FILE *fo = fopen("tpv101_rt.dat", "w");       /* x y rupture_time final_slip (km, s, m) */
    for(i = 0; i < nx; i++) for(j = 0; j < nz; j++)
      fprintf(fo, "%.4f %.4f %.5f %.5f\n", (xlo + (i + 0.5)*dx)/1e3, (ylo + (j + 0.5)*dx)/1e3, rt[i*nz + j], f->delta[i*nz + j]);
    fclose(fo); }
  i = (argc > 3 && strcmp(argv[3], "-")) ? compare(argv[3], 3, nm, rec, nsteps, dt, 1e-8, "tpv101") : 0;
  free(rt);
  fd3d_free(f); fd3d_tab_free(t); free(rec); free(load); free(Fr);
  return i;
}

/* mirror check: a fault reaching the surface with its image (fd3d_create2, mirror = 1)
   against the explicit doubled fault (2 nzp rows, symmetric initial stress) in a whole
   space; both use the same tables, so the solutions must agree to roundoff */
static int run_mirror(int argc, char **argv)
{
  double cL = 6000, cT = 3464, rho = 2670, dx = argc > 2 ? atof(argv[2]) : 500.0, dt = 0.5*dx/cT, tend = argc > 3 ? atof(argv[3]) : 6.0;
  int nx = (int)floor(20e3/dx + 0.5), nzp = (int)floor(5e3/dx + 0.5), i, j, n, nsteps = (int)floor(tend/dt + 0.5), fail;
  struct fd3d_tab *t = fd3d_tab_create(nx, 2*nzp, dx, dt, cL, cT, rho, 0);
  struct fd3d *fm = fd3d_create2(t, nzp, 1), *fd = fd3d_create(t);
  double emax = 0, vmax = 0, dmax = 0, ed = 0;
  fd->law = fm->law = 0; fd->mu_s = fm->mu_s = 0.677; fd->mu_d = fm->mu_d = 0.525; fd->d0 = fm->d0 = 0.40;
  for(i = 0; i < nx; i++) for(j = 0; j < nzp; j++){
    double x = -10e3 + (i + 0.5)*dx, z = (j + 0.5)*dx, tau = (fabs(x) < 1.5e3 && z < 3e3) ? 81.6e6 : 70e6;
    fm->tau0[i*nzp + j] = tau; fm->sigma[i*nzp + j] = 120e6;
    fd->tau0[i*2*nzp + nzp + j] = tau; fd->tau0[i*2*nzp + nzp - 1 - j] = tau;
    fd->sigma[i*2*nzp + nzp + j] = 120e6; fd->sigma[i*2*nzp + nzp - 1 - j] = 120e6;
  }
  for(n = 0; n < nsteps; n++){
    fd3d_step(fm, NULL); fd3d_step(fd, NULL);
    for(i = 0; i < nx; i++) for(j = 0; j < nzp; j++){
      double a = fm->V[i*nzp + j], b = fd->V[i*2*nzp + nzp + j];
      if(fabs(a - b) > emax) emax = fabs(a - b);
      if(fabs(b) > vmax) vmax = fabs(b);
      if(fabs(fm->delta[i*nzp + j] - fd->delta[i*2*nzp + nzp + j]) > ed) ed = fabs(fm->delta[i*nzp + j] - fd->delta[i*2*nzp + nzp + j]);
      if(fabs(fd->delta[i*2*nzp + nzp + j]) > dmax) dmax = fabs(fd->delta[i*2*nzp + nzp + j]);
    }
  }
  fail = !(emax <= 1e-8*vmax && ed <= 1e-8*dmax && dmax > 0);
  printf("fd3d_test mirror: %d x %d elements of %g m, %d steps: max |V diff| %.2e (of %.3f m/s), max |slip diff| %.2e (of %.3f m)\n",
         nx, nzp, dx, nsteps, emax, vmax, ed, dmax);
  printf("fd3d_test mirror: %s\n", fail ? "FAIL" : "PASS");
  fd3d_free(fm); fd3d_free(fd); fd3d_tab_free(t);
  return fail;
}

int main(int argc, char **argv)
{
  double cL = 6000, cT = 3464, rho = 2670, dx = argc > 2 ? atof(argv[2]) : 500.0, dt;
  int nx, nz, i, j, n, nsteps, s, ci[2], cj[2], ni, nj, cnt = 0, fail = 0;
  struct fd3d_tab *t; struct fd3d *f; FILE *fr = NULL;
  const char *nm[2] = {"PI", "PA"}; double sx[2] = {7.5e3, 0.0}, sz[2] = {0.0, 6.0e3};
  int scell[2][4][2], sn[2];
  double emax[3] = {0, 0, 0}, smax[3] = {0, 0, 0};
  double *rec = NULL; int nrec = 0;
  if(argc >= 3 && !strcmp(argv[1], "tpv101")) return run_tpv101(argc, argv);
  if(argc >= 2 && !strcmp(argv[1], "mirror")) return run_mirror(argc, argv);
  if(argc < 3 || strcmp(argv[1], "tpv3")){ fprintf(stderr, "usage: %s tpv3|tpv101 dx [ref] [tend]\n", argv[0]); return 2; }
  nx = (int)floor(30e3/dx + 0.5); nz = (int)floor(15e3/dx + 0.5); dt = 0.5*dx/cT;
  t = fd3d_tab_create(nx, nz, dx, dt, cL, cT, rho, 0);
  fprintf(stderr, "fd3d_test tpv3: %d x %d elements, dt %.4f s, window %d lags, FFT %d x %d, tables %.0f MB, k1 self %.6e, ks self %.6e\n",
          nx, nz, dt, t->L, t->Mx, t->Mz, 16.0*t->L*t->Mx*t->Mz/1e6, t->k1[1][1], t->ks_self);
  f = fd3d_create(t); f->law = 0; f->mu_s = 0.677; f->mu_d = 0.525; f->d0 = 0.40;
  for(i = 0; i < nx; i++) for(j = 0; j < nz; j++){
    double x = -15e3 + (i + 0.5)*dx, z = -7.5e3 + (j + 0.5)*dx;
    f->tau0[i*nz + j] = (fabs(x) < 1.5e3 && fabs(z) < 1.5e3) ? 81.6e6 : 70e6;
    f->sigma[i*nz + j] = 120e6;
  }
  for(s = 0; s < 2; s++){
    ni = near2(sx[s], -15e3, dx, nx, ci); nj = near2(sz[s], -7.5e3, dx, nz, cj);
    sn[s] = 0;
    for(i = 0; i < ni; i++) for(j = 0; j < nj; j++){ scell[s][sn[s]][0] = ci[i]; scell[s][sn[s]][1] = cj[j]; sn[s]++; }
  }
  nsteps = (int)floor(12.0/dt + 0.5);
  rec = malloc(sizeof(double)*2*4*nsteps);
  for(n = 0; n < nsteps; n++){
    fd3d_step(f, NULL);
    for(s = 0; s < 2; s++){
      double d = 0, v = 0, ta = 0; int k;
      for(k = 0; k < sn[s]; k++){ int idx = scell[s][k][0]*nz + scell[s][k][1]; d += f->delta[idx]; v += f->V[idx]; ta += f->tau[idx]; }
      rec[(s*nsteps + n)*4 + 0] = f->time; rec[(s*nsteps + n)*4 + 1] = d/sn[s]; rec[(s*nsteps + n)*4 + 2] = v/sn[s]; rec[(s*nsteps + n)*4 + 3] = ta/sn[s];
    }
  }
  nrec = nsteps;
  { FILE *fo = fopen("c_series.txt", "w"); for(s = 0; s < 2; s++) for(n = 0; n < nrec; n++) fprintf(fo, "%s %.10e %.12e %.12e %.12e\n", nm[s], rec[(s*nsteps+n)*4], rec[(s*nsteps+n)*4+1], rec[(s*nsteps+n)*4+2], rec[(s*nsteps+n)*4+3]); fclose(fo); }
  for(s = 0; s < 2; s++){
    double on = -1, pk = 0;
    for(n = 0; n < nrec; n++){ double v = rec[(s*nsteps + n)*4 + 2]; if(on < 0 && v > 0.1) on = rec[(s*nsteps + n)*4]; if(v > pk) pk = v; }
    printf("%s: V > 0.1 m/s at %.3f s, peak slip rate %.3f m/s, final slip %.3f m\n", nm[s], on, pk, rec[(s*nsteps + nrec - 1)*4 + 1]);
  }
  if(argc > 3 && (fr = fopen(argv[3], "r"))){
    char name[8]; double r[4];
    while(fscanf(fr, "%7s %lf %lf %lf %lf", name, r, r + 1, r + 2, r + 3) == 5){
      int q; s = strcmp(name, "PI") ? 1 : 0;
      n = (int)floor(r[0]/dt + 0.5) - 1;
      if(n < 0 || n >= nrec) continue;
      for(q = 0; q < 3; q++){
        double c = rec[(s*nsteps + n)*4 + 1 + q], e = fabs(c - r[1 + q]);
        if(e > emax[q]) emax[q] = e;
        if(fabs(r[1 + q]) > smax[q]) smax[q] = fabs(r[1 + q]);
      }
      cnt++;
    }
    fclose(fr);
    printf("fd3d_test tpv3 dx %g vs %s: %d samples; max |diff| slip %.2e (of %.3f), slip rate %.2e (of %.3f), stress %.2e (of %.3e)\n",
           dx, argv[3], cnt, emax[0], smax[0], emax[1], smax[1], emax[2], smax[2]);
    fail = (cnt == 0 || emax[0] > 1e-8*smax[0] || emax[1] > 1e-8*smax[1] || emax[2] > 1e-8*smax[2]);
    printf("fd3d_test tpv3: %s\n", fail ? "FAIL" : "PASS");
  }
  fd3d_free(f); fd3d_tab_free(t); free(rec);
  return fail;
}
