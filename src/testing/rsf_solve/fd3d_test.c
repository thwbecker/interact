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
  (void)n;
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
  /* the first three stations are those of the prototype reference; st120dp030 is a SCEC station
     (Kaneko et al.'s SPECFEM3D series) added for the external comparison */
  const char *nm[4] = {"st000dp075", "st090dp075", "st000dp030", "st120dp030"};
  double sx[4] = {0, 9e3, 0, 12e3}, sy[4] = {7.5e3, 7.5e3, 3e3, 3e3}, *rec;
  int sc[4][4]; double sw[4][4];
  struct fd3d_tab *t; struct fd3d *f;
  if(getenv("FD3D_SMALL")){ xlo = -19.5e3; ylo = -4.5e3; xlen = 39e3; ylen = 24e3; }   /* 1.5 km of strengthening fault beyond the transition */
  nx = (int)floor(xlen/dx + 0.5); nz = (int)floor(ylen/dx + 0.5); dt = 0.5*dx/cT; N = nx*nz;
  t = fd3d_tab_create(nx, nz, dx, dt, cL, cT, rho, 0);
  fprintf(stderr, "fd3d_test tpv101: %d x %d elements, dt %.4f s, window %d lags, FFT %d x %d\n", nx, nz, dt, t->L, t->Mx, t->Mz);
  f = fd3d_create(t); f->law = 1; f->f0 = 0.6; f->V0 = 1e-6;
  f->nl = getenv("FD3D_NL") ? 1 : 0;
  f->generic = getenv("FD3D_GENERIC") ? 1 : 0;
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
  for(s = 0; s < 4; s++){
    double fi = (sx[s] - xlo)/dx - 0.5, fj = (sy[s] - ylo)/dx - 0.5; int i0 = (int)floor(fi), j0 = (int)floor(fj); double wi = fi - i0, wj = fj - j0;
    sc[s][0] = i0*nz + j0; sw[s][0] = (1 - wi)*(1 - wj); sc[s][1] = (i0 + 1)*nz + j0; sw[s][1] = wi*(1 - wj);
    sc[s][2] = i0*nz + j0 + 1; sw[s][2] = (1 - wi)*wj; sc[s][3] = (i0 + 1)*nz + j0 + 1; sw[s][3] = wi*wj;
  }
  { int K = getenv("FD3D_K") ? atoi(getenv("FD3D_K")) : 1;    /* fixed substeps per step (test) */
  nsteps = (int)floor(tend/(dt*K) + 0.5); rec = malloc(sizeof(double)*4*4*nsteps);
  rt = malloc(sizeof(double)*N); for(k = 0; k < N; k++) rt[k] = -1.0;
  for(n = 0; n < nsteps; n++){
    double t1 = f->time + K*dt, G = (t1 > 0 && t1 < T) ? exp((t1 - T)*(t1 - T)/(t1*(t1 - 2*T))) : (t1 >= T ? 1.0 : 0.0);
    for(k = 0; k < N; k++) load[k] = 25e6*Fr[k]*G;
    if(fd3d_step_k(f, load, K) < 0){ fprintf(stderr, "fd3d_step_k failed\n"); return 1; }
    for(k = 0; k < N; k++) if(rt[k] < 0 && fabs(f->V[k]) > 1e-3) rt[k] = f->time;
    if(n % 50 == 0){ double vm = 0; for(k = 0; k < N; k++) if(fabs(f->V[k]) > vm) vm = fabs(f->V[k]); fprintf(stderr, "  t %6.3f max V %.3e\n", f->time, vm); }
    for(s = 0; s < 4; s++){
      double d = 0, v = 0, ta = 0;
      for(k = 0; k < 4; k++){ d += sw[s][k]*f->delta[sc[s][k]]; v += sw[s][k]*f->V[sc[s][k]]; ta += sw[s][k]*f->tau[sc[s][k]]; }
      rec[(s*nsteps + n)*4] = f->time; rec[(s*nsteps + n)*4 + 1] = d; rec[(s*nsteps + n)*4 + 2] = v; rec[(s*nsteps + n)*4 + 3] = ta;
    }
  }
  for(s = 0; s < 4; s++){
    double on = -1, pk = 0;
    for(n = 0; n < nsteps; n++){ double v = rec[(s*nsteps + n)*4 + 2]; if(on < 0 && v > 1e-3) on = rec[(s*nsteps + n)*4]; if(v > pk) pk = v; }
    printf("%s: rupture time %.3f s, peak slip rate %.3f m/s, slip %.3f m\n", nm[s], on, pk, rec[(s*nsteps + nsteps - 1)*4 + 1]);
    if(getenv("FD3D_STFILES")){                       /* time, slip, slip rate, shear stress (MPa) */
      char fn[64]; FILE *fs; sprintf(fn, "tpv101_%s.dat", nm[s]); fs = fopen(fn, "w");
      for(n = 0; n < nsteps; n++) fprintf(fs, "%.6e %.6e %.6e %.6e\n", rec[(s*nsteps + n)*4], rec[(s*nsteps + n)*4 + 1], rec[(s*nsteps + n)*4 + 2], rec[(s*nsteps + n)*4 + 3]/1e6);
      fclose(fs);
    }
  }
  { FILE *fo = fopen("tpv101_rt.dat", "w");       /* x y rupture_time final_slip (km, s, m) */
    for(i = 0; i < nx; i++) for(j = 0; j < nz; j++)
      fprintf(fo, "%.4f %.4f %.5f %.5f\n", (xlo + (i + 0.5)*dx)/1e3, (ylo + (j + 0.5)*dx)/1e3, rt[i*nz + j], f->delta[i*nz + j]);
    fclose(fo); }
  i = (argc > 3 && strcmp(argv[3], "-")) ? compare(argv[3], 3, nm, rec, nsteps, dt*K, 1e-8, "tpv101") : 0;
  }
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
  fd->nl = fm->nl = getenv("FD3D_NL") ? 1 : 0;
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

/* free-surface layer check: vertical fault 2 x 1 km of 500 m reaching the surface, mirror
   plus a 16 x 16 layer of 250 m opening elements, prescribed slip rate 1 for 0.4 s, NL;
   fault traction history against fully_dynamic/fs3d_dyn.py */
static int run_layer(const char *ref)
{
  double cL = 6000, cT = 3464, rho = 2670, dx = 500.0, dt = 0.5*dx/cT, *vp, emax = 0.0, smax = 0.0, r;
  int nx = 4, nzp = 2, N = nx*nzp, n, i, nsteps = (int)floor(3.5/dt + 0.5), cnt = 0, fail;
  struct fd3d_tab *t = fd3d_tab_create(nx, 2*nzp, dx, dt, cL, cT, rho, 45);
  struct fd3d *f = fd3d_create2(t, nzp, 1);
  FILE *fr = fopen(ref, "r");
  f->nl = 1; f->law = 0; f->mu_s = f->mu_d = 0.0; f->d0 = 1.0;
  if(getenv("FD3D_LAYER_FFT")) f->layer = fd3d_layer_create_fft(f, -1000.0, 2, 16, 16, -2000.0, 45);
  else f->layer = fd3d_layer_create(f, -1000.0, 16, 16, 250.0, 250.0, -2000.0, -2000.0, 45);
  if(!f->layer || !fr){ fprintf(stderr, "layer setup or reference failed\n"); return 1; }
  vp = malloc(sizeof(double)*N); f->vpres = vp;
  for(n = 0; n < nsteps; n++){
    for(i = 0; i < N; i++) vp[i] = (f->time + dt <= 0.4 + 1e-12) ? 1.0 : 0.0;
    fd3d_step(f, NULL);
    for(i = 0; i < N; i++){
      if(fscanf(fr, "%lf", &r) != 1) break;
      if(fabs(f->tau[i] - r) > emax) emax = fabs(f->tau[i] - r);
      if(fabs(r) > smax) smax = fabs(r);
      cnt++;
    }
  }
  fclose(fr);
  { extern double fd3d_layer_prof[4]; if(getenv("FD3D_PROF")) fprintf(stderr, "layer profile: known %.3f s, fault solves %.3f s, layer solves %.3f s, finish %.3f s\n",
       fd3d_layer_prof[0], fd3d_layer_prof[1], fd3d_layer_prof[2], fd3d_layer_prof[3]); }
  fail = !(cnt == nsteps*N && emax <= 1e-8*smax);
  printf("fd3d_test layer: %d values, max |tau diff| %.2e of %.3e\n", cnt, emax, smax);
  printf("fd3d_test layer: %s\n", fail ? "FAIL" : "PASS");
  fd3d_free(f); fd3d_tab_free(t); free(vp);
  return fail;
}

/* full-surface check: 2 x 1 km fault of 500 m dipping 30 degrees, top at 500 m, dip slip, a
   24 x 35 surface of 250 m with three unknowns per element, prescribed slip rate 1 for 0.4 s,
   NL; shear and normal traction histories against fully_dynamic/fs3d_dip.py */
static int run_surf(const char *ref)
{
  double cL = 6000, cT = 3464, rho = 2670, dx = 500.0, dt = 0.5*dx/cT, *vp, es = 0.0, ss = 0.0, en = 0.0, sn = 0.0, r;
  double dip = 30.0, ny = 35, ymid = -0.5*1000.0*cos(dip*M_PI/180.0);
  int nx = 4, nz = 2, N = nx*nz, n, i, nsteps = (int)floor(3.5/dt + 0.5), cnt = 0, fail;
  struct fd3d_tab *t = fd3d_tab_create_c(nx, nz, dx, dt, cL, cT, rho, 45, 2);     /* dip slip along grid axis 2 */
  struct fd3d *f = fd3d_create(t);
  FILE *fr = fopen(ref, "r");
  f->nl = 1; f->law = 0; f->mu_s = f->mu_d = 0.0; f->d0 = 1.0;
  f->layer = fd3d_surf_create(f, -1000.0, 0.0, -500.0, dip, 2, 2, 24, (int)ny, -3000.0, ymid - ny*125.0, 45, 1.0);
  if(!f->layer || !fr){ fprintf(stderr, "surface setup or reference failed\n"); return 1; }
  vp = malloc(sizeof(double)*N); f->vpres = vp;
  for(n = 0; n < nsteps; n++){
    for(i = 0; i < N; i++) vp[i] = (f->time + dt <= 0.4 + 1e-12) ? 1.0 : 0.0;
    fd3d_step(f, NULL);
    for(i = 0; i < N; i++){ if(fscanf(fr, "%lf", &r) != 1) break; if(fabs(f->tau[i] - r) > es) es = fabs(f->tau[i] - r); if(fabs(r) > ss) ss = fabs(r); cnt++; }
    for(i = 0; i < N; i++){ if(fscanf(fr, "%lf", &r) != 1) break; if(fabs(f->layer->tn[i] - r) > en) en = fabs(f->layer->tn[i] - r); if(fabs(r) > sn) sn = fabs(r); cnt++; }
  }
  fclose(fr);
  fail = !(cnt == 2*nsteps*N && es <= 1e-8*ss && en <= 1e-8*sn);
  printf("fd3d_test surf: %d values, max |shear diff| %.2e of %.3e, max |normal diff| %.2e of %.3e\n", cnt, es, ss, en, sn);
  printf("fd3d_test surf: %s\n", fail ? "FAIL" : "PASS");
  fd3d_free(f); fd3d_tab_free(t); free(vp);
  return fail;
}

/* stability of the virtual surface alone (fault at rest, random initial surface rates, NL):
   the largest surface rate must decay by a factor 1e3 over 300 steps with surface elements
   of dx/m; m = 4 (the rsf_solve default) passes, m = 2 grows (fd_plan.md) */
static int run_surfstab(int m)
{
  double dx = 1000.0, dt = 0.5*dx/3464.0, dip = 30.0, c = cos(dip*M_PI/180.0), s = sin(dip*M_PI/180.0), v0 = 0.0, v1 = 0.0;
  int nx = 12, nz = 6, mg = 4, i, n, fail;
  double hs = dx/m, ytop = 0.5*dx*c, ztop = -750.0 + 0.5*dx*s, Wd = nz*dx*c, ys0 = ytop - Wd - mg*dx;
  int nys, nxs = (nx + 2*mg)*m;
  struct fd3d_tab *t = fd3d_tab_create_c(nx, nz, dx, dt, 6000.0, 3464.0, 2670.0, 50, 2);
  struct fd3d *f = fd3d_create(t);
  double *vp = calloc(nx*nz, sizeof(double));
  ys0 = ytop - ceil((ytop - ys0)/hs)*hs; nys = (int)ceil((ytop + mg*dx - ys0)/hs);
  f->nl = 1; f->law = 0; f->mu_s = f->mu_d = 0.0; f->d0 = 1.0; f->vpres = vp;
  f->layer = fd3d_surf_create(f, -0.5*dx, ytop, ztop, dip, 2, m, nxs, nys, -0.5*dx - mg*dx, ys0, 50, 1.0);
  srand(1);
  for(i = 0; i < 3*f->layer->ns; i++){ f->layer->V3[i] = rand()/(double)RAND_MAX - 0.5; if(fabs(f->layer->V3[i]) > v0) v0 = fabs(f->layer->V3[i]); }
  for(n = 0; n < 300; n++) fd3d_step(f, NULL);
  for(i = 0; i < 3*f->layer->ns; i++) if(fabs(f->layer->V3[i]) > v1 || !isfinite(f->layer->V3[i])) v1 = isfinite(f->layer->V3[i]) ? fabs(f->layer->V3[i]) : 1e300;
  fail = !(v1 <= 1e-3*v0);
  printf("fd3d_test surfstab: surface elements dx/%d, max rate %.2e -> %.2e after 300 steps\n", m, v0, v1);
  printf("fd3d_test surfstab: %s\n", fail ? "FAIL" : "PASS");
  fd3d_free(f); fd3d_tab_free(t); free(vp);
  return fail;
}

/* surface stability analysis data: for a surface with elements dx/m at the given CFL, write the
   symbols of the x-slip -> xz-traction tables (Ks dt and W(l), l = 1..L) at wavenumbers
   (i/8, j/8) of Nyquist, i = 0..8, j = 0..4, for fully_dynamic/fs3d_surf_stability.py */
static int run_surfsym(double cfl, int m, const char *out)
{
  double dx = 1000.0, dt = cfl*dx/3464.0, dip = 30.0, c = cos(dip*M_PI/180.0), s = sin(dip*M_PI/180.0);
  int nx = 12, nz = 6, mg = 4, i, j, l;
  double hs = dx/m, ytop = 0.5*dx*c, ztop = -750.0 + 0.5*dx*s, Wd = nz*dx*c, ys0 = ytop - Wd - mg*dx;
  int nys, nxs = (nx + 2*mg)*m, L = 60;
  struct fd3d_tab *t = fd3d_tab_create_c(nx, nz, dx, dt, 6000.0, 3464.0, 2670.0, L, 2);
  struct fd3d *f = fd3d_create(t);
  struct fd3d_layer *ly;
  FILE *fo = fopen(out, "w");
  ys0 = ytop - ceil((ytop - ys0)/hs)*hs; nys = (int)ceil((ytop + mg*dx - ys0)/hs);
  f->nl = 1; f->law = 0; f->mu_s = f->mu_d = 0.0; f->d0 = 1.0;
  ly = fd3d_surf_create(f, -0.5*dx, ytop, ztop, dip, 2, m, nxs, nys, -0.5*dx - mg*dx, ys0, L, 1.0);
  if(!ly || !fo){ fprintf(stderr, "surfsym: setup failed\n"); return 1; }
  fprintf(fo, "# cS dt/hs %.4f, dt %.6e, eta %.6e, L %d; columns: i j Ks*dt W(1..L)  (wavenumber i/8, j/8 of Nyquist)\n", m*cfl, dt, 0.5*2670.0*3464.0, L);
  for(i = 0; i <= 8; i++) for(j = 0; j <= 4; j++){
    size_t pq = (size_t)((i*ly->Mx2/16) % ly->Mx2)*ly->Mh2 + (j*ly->My2/8) % ly->Mh2;
    fprintf(fo, "%d %d %.10e", i, j, creal(ly->T2[0][1][pq])*dt);
    for(l = 0; l < L; l++) fprintf(fo, " %.10e", creal(ly->T2[0][0][(size_t)l*ly->Mx2*ly->Mh2 + pq]));
    fprintf(fo, "\n");
  }
  fclose(fo); f->layer = ly; fd3d_free(f); fd3d_tab_free(t);
  printf("fd3d_test surfsym: wrote %s\n", out);
  return 0;
}

/* quasi-2-D check of the dipping fault with the full surface and normal-stress coupling: a 16 km
   long thrust dipping 30 degrees (6 km down dip, 1 km squares, top at 500 m), slip weakening,
   nucleation on the bottom 2 km along the whole strike; the central column's final slip and
   rupture times against the 2-D plane-strain prototype (ref_fd3d_q2d.txt), 2 percent in slip
   (the remaining along-strike edge effect at 16 km) and one step in time */
static int run_q2d(const char *ref)
{
  double dx = 1000.0, dt = 0.5*dx/3464.0, dip = 30.0, c = cos(dip*M_PI/180.0), top = 500.0, T = 6.0;
  int nx = 16, nz = 6, N = nx*nz, i, n, m = 4, mg = 2, nsteps = (int)(T/dt), fail = 0, cnt = 0;
  double hs = dx/m, ytop = 0.5*dx*c, ztop = -top, Wd = nz*dx*c, ys0 = ytop - Wd - mg*dx, emax = 0.0, tmax = 0.0;
  int nys, nxs = (nx + 2*mg)*m, L;
  struct fd3d_tab *t; struct fd3d *f; double *rt = malloc(sizeof(double)*N); FILE *fr = fopen(ref, "r"); char line[512];
  ys0 = ytop - ceil((ytop - ys0)/hs)*hs; nys = (int)ceil((ytop + mg*dx - ys0)/hs);
  L = (int)(sqrt(pow((nx + 2*mg)*dx, 2) + pow(Wd + 2*mg*dx, 2))/(3464.0*dt)) + 4;
  t = fd3d_tab_create_c(nx, nz, dx, dt, 6000.0, 3464.0, 2670.0, L, 2); f = fd3d_create(t);
  f->nl = 1; f->law = 0; f->mu_s = 0.677; f->mu_d = 0.525; f->d0 = 0.40;
  for(i = 0; i < N; i++){ f->sigma[i] = 120e6; f->tau0[i] = ((i % nz) >= nz - 2) ? 81.6e6 : 70e6; rt[i] = -1.0; }   /* thrust: + slip up dip */
  f->layer = fd3d_surf_create(f, -0.5*dx, ytop, ztop, dip, 2, m, nxs, nys, -0.5*dx - mg*dx, ys0, L, 1.0);
  if(!f->layer || !fr){ fprintf(stderr, "q2d: setup or reference failed\n"); return 1; }
  for(i = 0; i < N; i++) f->layer->sig0[i] = 120e6;
  for(n = 0; n < nsteps; n++){ fd3d_step(f, NULL); for(i = 0; i < N; i++) if(rt[i] < 0 && fabs(f->V[i]) > 0.1) rt[i] = f->time; }
  while(fgets(line, sizeof(line), fr)){
    double sl, tr; int jj;
    if(line[0] == '#' || sscanf(line, "%d %lf %lf", &jj, &sl, &tr) != 3) continue;
    i = (nx/2)*nz + jj;
    if(fabs(fabs(f->delta[i]) - sl)/sl > emax) emax = fabs(fabs(f->delta[i]) - sl)/sl;
    if(fabs(rt[i] - tr) > tmax) tmax = fabs(rt[i] - tr);
    cnt++;
  }
  fclose(fr);
  fail = !(cnt == nz && emax <= 0.02 && tmax <= 1.01*dt);
  printf("fd3d_test q2d: central column against the 2-D prototype: max relative slip difference %.2e, max rupture-time difference %.3f s (step %.3f s)\n", emax, tmax, dt);
  printf("fd3d_test q2d: %s\n", fail ? "FAIL" : "PASS");
  fd3d_free(f); fd3d_tab_free(t); free(rt);
  return fail;
}

/* SCEC TPV10 (60 degree normal fault reaching the surface, slip weakening with cohesion, tension
   as zero normal stress, slip path) with the full virtual surface, the extension row above the
   trace and normal-stress coupling; q2d = 1 nucleates along the whole strike (a plane-strain
   analogue for comparison with Dunham's 2-D runs). Writes <out>_st<down-dip m>.dat for the centre
   column (time, up-dip slip, up-dip slip rate, up-dip shear stress, normal stress; SCEC signs:
   normal stress negative in compression, slip positive up dip) at the requested down-dip
   distances, interpolated between element centres. */
static int run_tpv10(double dx, double Ls, int q2d, int mg, const char *out)
{
  double cL = 5716.0, cT = 3300.0, rho = 2700.0, dt = 0.5*dx/cT, dip = 60.0, T = 15.0;
  double c = cos(dip*M_PI/180.0), s = sin(dip*M_PI/180.0), W = 15000.0, mus = 0.760, mud = 0.448, d0 = 0.5, coh = 0.2e6;
  int nx = (int)(Ls/dx + 0.5), nz = (int)(W/dx + 0.5), N = nx*nz, i, j, n, m = 1, nsteps = (int)(T/dt + 0.5), q, nst = 4;
  double ytop = 0.0, ztop = 0.0, Wd = nz*dx*c, hs = dx/m, ys0, *yc, *wy;
  double st_dd[4] = {0.0, 4500.0, 7500.0, 12000.0};
  int nys, nxs = (nx + 2*mg)*m, L;
  struct fd3d_tab *t; struct fd3d *f; FILE *fo[4]; char fn[256];
  /* surface rows: uniform hs over the projection plus mg fault elements either side, the trace on
     a row boundary */
  ys0 = ytop - Wd - mg*dx; nys = (int)((ytop + mg*dx - ys0)/hs + 0.5);
  /* window: the tables vanish once the S wave has crossed the farthest pair */
  L = (int)ceil(sqrt(pow((nx + 2*mg)*dx, 2) + pow(Wd + 2*mg*dx, 2) + pow(W*s, 2))/(cT*dt)) + 4;
  if(L > nsteps + 2) L = nsteps + 2;
  yc = malloc(sizeof(double)*nys); wy = malloc(sizeof(double)*nys);
  for(q = 0; q < nys; q++){ yc[q] = ys0 + (q + 0.5)*hs; wy[q] = hs; }
  if(getenv("TPV10_NOSURF")){ t = fd3d_tab_create_c(nx, nz, dx, dt, cL, cT, rho, L, 2); f = fd3d_create(t); }
  else{ t = fd3d_tab_create_c(nx, nz + 1, dx, dt, cL, cT, rho, L, 2); f = fd3d_create3(t, nz, 0, 1); }
  f->nl = 1; f->law = 0; f->mu_s = mus; f->mu_d = mud; f->d0 = d0;
  f->coh = malloc(sizeof(double)*N); f->path = 1; f->dpath = calloc(N, sizeof(double)); f->sw_saturate = 1;
  if(!getenv("TPV10_NOSURF")){
    f->layer = fd3d_surf_create_rows(f, -0.5*dx, ytop, ztop, dip, 2, m, nxs, nys, -0.5*dx - mg*dx, yc, wy, L, 1.0);
    if(!f->layer){ fprintf(stderr, "tpv10: surface setup failed\n"); return 1; }
    if(getenv("TPV10_NONC")) f->layer->normal_coupling = 0;
  }
  fprintf(stderr, "tpv10: %d x %d elements of %g m, %d steps of %.4f s, surface %d x %d of %g m, %d lags\n", nx, nz, dx, nsteps, dt, nxs, nys, hs, L);
  for(i = 0; i < nx; i++) for(j = 0; j < nz; j++){
    int e = i*nz + j;
    double sdd = (j + 0.5)*dx, xs = (i + 0.5)*dx - 0.5*nx*dx, sig = 7378.0*sdd, tau;
    int nucl = (fabs(sdd - 12000.0) <= 1500.0 + 1e-6) && (q2d || fabs(xs) <= 1500.0 + 1e-6);
    tau = nucl ? coh + (mus + 0.0057)*sig : 0.55*sig;
    f->tau0[e] = -tau;                 /* normal faulting: the hanging wall moves down dip (negative up-dip slip here) */
    f->sigma[e] = sig; if(f->layer) f->layer->sig0[e] = sig; f->coh[e] = coh;
  }
  for(q = 0; q < nst; q++){ sprintf(fn, "%s_st%05.0f.dat", out, st_dd[q]); fo[q] = fopen(fn, "w"); }
  for(n = 0; n < nsteps; n++){
    fd3d_step(f, NULL);
    for(q = 0; q < nst; q++){
      /* centre column (average of the two middle columns when nx is even), interpolated along dip */
      double sdd = st_dd[q], jj = sdd/dx - 0.5, w1; int j0, j1, ic[2], k, ncol;
      double v[4] = {0, 0, 0, 0};
      if(jj < 0) jj = 0;
      if(jj > nz - 1) jj = nz - 1;
      j0 = (int)floor(jj); j1 = j0 + 1 < nz ? j0 + 1 : j0; w1 = jj - j0;
      if(nx % 2 == 0){ ic[0] = nx/2 - 1; ic[1] = nx/2; ncol = 2; } else { ic[0] = nx/2; ncol = 1; }
      for(k = 0; k < ncol; k++){
        int e0 = ic[k]*nz + j0, e1 = ic[k]*nz + j1;
        v[0] += ((1 - w1)*f->delta[e0] + w1*f->delta[e1])/ncol;
        v[1] += ((1 - w1)*f->V[e0] + w1*f->V[e1])/ncol;
        v[2] += ((1 - w1)*f->tau[e0] + w1*f->tau[e1])/ncol;
        v[3] += ((1 - w1)*f->sigma[e0] + w1*f->sigma[e1])/ncol;
      }
      fprintf(fo[q], "%.6e %.6e %.6e %.6e %.6e\n", f->time, -v[0], -v[1], -v[2]/1e6, -v[3]/1e6);
    }
  }
  for(q = 0; q < nst; q++) fclose(fo[q]);
  if(f->layer) fprintf(stderr, "tpv10: done; outer/inner iterations per step %.1f / %.1f\n", (double)f->layer->n_outer/f->layer->n_steps, (double)f->layer->n_inner/f->layer->n_steps);
  fd3d_free(f); fd3d_tab_free(t); free(yc); free(wy);
  return 0;
}

/* kernel check: full stress tensors (slip and opening, rectangles) against tada_quad.py values */
static int run_kernels(const char *ref)
{
  FILE *fr = fopen(ref, "r"); double v[13], s[6], emax = 0.0, smax = 0.0; int n = 0, k;
  if(!fr){ fprintf(stderr, "cannot open %s\n", ref); return 1; }
  while(fscanf(fr, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", v, v + 1, v + 2, v + 3, v + 4, v + 5, v + 6,
                v + 7, v + 8, v + 9, v + 10, v + 11, v + 12) == 13){
    double sc = 1.0;
    fd3d_rect_ramp_tensor(v[1], v[2], v[3], v[4], v[5], v[6], 6000.0, 3464.0, 2670.0, (int)v[0], s);
    for(k = 0; k < 6; k++) if(fabs(v[7 + k]) > sc) sc = fabs(v[7 + k]);
    for(k = 0; k < 6; k++){
      double e = fabs(s[k] - v[7 + k])/sc;          /* relative to the point's largest component */
      if(!(e <= emax)) emax = isfinite(e) ? e : 1e300;
    }
    if(sc > smax) smax = sc;
    n++;
  }
  fclose(fr);
  printf("fd3d_test kernels: %d points (slip along x1 and x2, opening; rectangles, edge cases), max relative diff %.2e\n", n, emax);
  /* 1e-7: the closed forms cancel large terms at some points; Python and C (with fused
     multiply-add) differ at the 1e-8 level there */
  printf("fd3d_test kernels: %s\n", (n > 0 && emax <= 1e-7) ? "PASS" : "FAIL");
  return !(n > 0 && emax <= 1e-7);
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
  if(argc >= 3 && !strcmp(argv[1], "kernels")) return run_kernels(argv[2]);
  if(argc >= 3 && !strcmp(argv[1], "layer")) return run_layer(argv[2]);
  if(argc >= 3 && !strcmp(argv[1], "surf")) return run_surf(argv[2]);
  if(argc >= 2 && !strcmp(argv[1], "surfstab")) return run_surfstab(argc >= 3 ? atoi(argv[2]) : 4);
  if(argc >= 5 && !strcmp(argv[1], "surfsym")) return run_surfsym(atof(argv[2]), atoi(argv[3]), argv[4]);
  if(argc >= 3 && !strcmp(argv[1], "q2d")) return run_q2d(argv[2]);
  if(argc >= 7 && !strcmp(argv[1], "tpv10")) return run_tpv10(atof(argv[2]), atof(argv[3]), atoi(argv[4]), atoi(argv[5]), argv[6]);
  if(argc < 3 || strcmp(argv[1], "tpv3")){ fprintf(stderr, "usage: %s tpv3|tpv101 dx [ref] [tend]\n", argv[0]); return 2; }
  nx = (int)floor(30e3/dx + 0.5); nz = (int)floor(15e3/dx + 0.5); dt = 0.5*dx/cT;
  t = fd3d_tab_create(nx, nz, dx, dt, cL, cT, rho, 0);
  fprintf(stderr, "fd3d_test tpv3: %d x %d elements, dt %.4f s, window %d lags, FFT %d x %d, tables %.0f MB, k1 self %.6e, ks self %.6e\n",
          nx, nz, dt, t->L, t->Mx, t->Mz, 16.0*t->L*t->Mx*t->Mz/1e6, t->k1[1][1], t->ks_self);
  f = fd3d_create(t); f->law = 0; f->mu_s = 0.677; f->mu_d = 0.525; f->d0 = 0.40;
  f->nl = getenv("FD3D_NL") ? 1 : 0;
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
