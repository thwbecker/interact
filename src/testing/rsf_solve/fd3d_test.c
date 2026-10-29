/* fd3d_test: TPV3 (Day et al. 2005) with the 3-D C marcher, fixed rake,
   compared step by step with the prototype (run_tpv3.py, FIXED_RAKE=1).
   usage: fd3d_test tpv3 dx [ref_file] */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "fd3d.h"
/* cells for a station at x on a grid starting at edge lo: the cell containing x,
   or the two cells sharing an edge at x */
static int near2(double x, double lo, double dx, int n, int *c)
{
  (void)n;
  double f = (x - lo)/dx;
  int k = (int)floor(f + 0.5);
  if(fabs(f - k) < 1e-9){
    c[0] = k - 1;
    c[1] = k;
    return 2;
  }
  c[0] = (int)floor(f);
  return 1;
}
static int compare(const char *ref, int nst, const char **nm, double *rec, int nsteps, double dt, double tol, const char *what)
{
  FILE *fr = fopen(ref, "r");
  char name[32];
  double r[4], emax[3] = {0, 0, 0}, smax[3] = {0, 0, 0};
  int cnt = 0, s, n, q, fail;
  double e;
  if(!fr){
    fprintf(stderr, "cannot open %s\n", ref);
    return 1;
  }
  while(fscanf(fr, "%31s %lf %lf %lf %lf", name, r, r + 1, r + 2, r + 3) == 5){
    for(s = 0; s < nst; s++)
      if(!strcmp(name, nm[s]))
        break;
    if(s == nst)
      continue;
    n = (int)floor(r[0]/dt + 0.5) - 1;
    if(n < 0 || n >= nsteps)
      continue;
    for(q = 0; q < 3; q++){
      e = fabs(rec[(s*nsteps + n)*4 + 1 + q] - r[1 + q]);
      if(e > emax[q])
        emax[q] = e;
      if(fabs(r[1 + q]) > smax[q])
        smax[q] = fabs(r[1 + q]);
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
  int sc[4][4];
  double sw[4][4];
  struct fd3d_tab *t;
  struct fd3d *f;
  double x;
  double y;
  double Bx;
  double By;
  double r2;
  double ax;
  double ay;
  double fi;
  double fj;
  int i0;
  int j0;
  double wi;
  double wj;
  int K;
  double t1;
  double G;
  double vm;
  double d;
  double v;
  double ta;
  double on;
  double pk;
  double v_2;
  char fn[64];
  FILE *fs;
  extern double fd3d_prof3[4];
  extern long fd3d_cnt[3];
  FILE *fo;
  if(getenv("FD3D_SMALL")){
    xlo = -19.5e3;
    ylo = -4.5e3;
    xlen = 39e3;
    ylen = 24e3;
  }   /* 1.5 km of strengthening fault beyond the transition */
  nx = (int)floor(xlen/dx + 0.5);
  nz = (int)floor(ylen/dx + 0.5);
  dt = 0.5*dx/cT;
  N = nx*nz;
  t = fd3d_tab_create(nx, nz, dx, dt, cL, cT, rho, 0);
  fprintf(stderr, "fd3d_test tpv101: %d x %d elements, dt %.4f s, window %d lags, FFT %d x %d\n", nx, nz, dt, t->L, t->Mx, t->Mz);
  f = fd3d_create(t);
  f->law = 1;
  f->f0 = 0.6;
  f->V0 = 1e-6;
  f->nl = getenv("FD3D_NL") ? 1 : 0;
  f->generic = getenv("FD3D_GENERIC") ? 1 : 0;
  f->a = malloc(sizeof(double)*N);
  f->b = malloc(sizeof(double)*N);
  f->Lc = malloc(sizeof(double)*N);
  f->theta = malloc(sizeof(double)*N);
  load = malloc(sizeof(double)*N);
  Fr = malloc(sizeof(double)*N);
  for(i = 0; i < nx; i++)
    for(j = 0; j < nz; j++){
      x = xlo + (i + 0.5)*dx;
      y = ylo + (j + 0.5)*dx;
      ax = fabs(x);
      ay = fabs(y - 7.5e3);
      Bx = ax <= W ? 1.0 : (ax < W + w ? 0.5*(1.0 + tanh(w/(ax - W - w) + w/(ax - W))) : 0.0);
      By = ay <= W/2 ? 1.0 : (ay < W/2 + w ? 0.5*(1.0 + tanh(w/(ay - W/2 - w) + w/(ay - W/2))) : 0.0);
      k = i*nz + j;
      f->a[k] = 0.008 + 0.008*(1.0 - Bx*By);
      f->b[k] = 0.012;
      f->Lc[k] = 0.02;
      f->sigma[k] = 120e6;
      f->tau0[k] = 75e6;
      f->V[k] = 1e-12;
      f->theta[k] = f->Lc[k]/f->V0*exp((f->a[k]*log(2.0*sinh(75e6/(f->a[k]*120e6))) - f->f0 - f->a[k]*log(1e-12/f->V0))/f->b[k]);
      r2 = x*x + (y - 7.5e3)*(y - 7.5e3);
      Fr[k] = r2 < R*R ? exp(r2/(r2 - R*R)) : 0.0;
    }
  for(s = 0; s < 4; s++){
    fi = (sx[s] - xlo)/dx - 0.5;
    fj = (sy[s] - ylo)/dx - 0.5;
    i0 = (int)floor(fi);
    j0 = (int)floor(fj);
    wi = fi - i0;
    wj = fj - j0;
    sc[s][0] = i0*nz + j0;
    sw[s][0] = (1 - wi)*(1 - wj);
    sc[s][1] = (i0 + 1)*nz + j0;
    sw[s][1] = wi*(1 - wj);
    sc[s][2] = i0*nz + j0 + 1;
    sw[s][2] = (1 - wi)*wj;
    sc[s][3] = (i0 + 1)*nz + j0 + 1;
    sw[s][3] = wi*wj;
  }
  {
    K = getenv("FD3D_K") ? atoi(getenv("FD3D_K")) : 1;    /* fixed substeps per step (test) */
    nsteps = (int)floor(tend/(dt*K) + 0.5);
    rec = malloc(sizeof(double)*4*4*nsteps);
    rt = malloc(sizeof(double)*N);
    for(k = 0; k < N; k++)
      rt[k] = -1.0;
    for(n = 0; n < nsteps; n++){
      t1 = f->time + K*dt;
      G = (t1 > 0 && t1 < T) ? exp((t1 - T)*(t1 - T)/(t1*(t1 - 2*T))) : (t1 >= T ? 1.0 : 0.0);
      for(k = 0; k < N; k++)
        load[k] = 25e6*Fr[k]*G;
      if(fd3d_step_k(f, load, K) < 0){
        fprintf(stderr, "fd3d_step_k failed\n");
        return 1;
      }
      for(k = 0; k < N; k++)
        if(rt[k] < 0 && fabs(f->V[k]) > 1e-3)
          rt[k] = f->time;
      if(n % 50 == 0){
        vm = 0;
        for(k = 0; k < N; k++)
          if(fabs(f->V[k]) > vm)
            vm = fabs(f->V[k]);
        fprintf(stderr, "  t %6.3f max V %.3e\n", f->time, vm);
      }
      for(s = 0; s < 4; s++){
        d = 0;
        v = 0;
        ta = 0;
        for(k = 0; k < 4; k++){
          d += sw[s][k]*f->delta[sc[s][k]];
          v += sw[s][k]*f->V[sc[s][k]];
          ta += sw[s][k]*f->tau[sc[s][k]];
        }
        rec[(s*nsteps + n)*4] = f->time;
        rec[(s*nsteps + n)*4 + 1] = d;
        rec[(s*nsteps + n)*4 + 2] = v;
        rec[(s*nsteps + n)*4 + 3] = ta;
      }
    }
    for(s = 0; s < 4; s++){
      on = -1;
      pk = 0;
      for(n = 0; n < nsteps; n++){
        v_2 = rec[(s*nsteps + n)*4 + 2];
        if(on < 0 && v_2 > 1e-3)
          on = rec[(s*nsteps + n)*4];
        if(v_2 > pk)
          pk = v_2;
      }
      printf("%s: rupture time %.3f s, peak slip rate %.3f m/s, slip %.3f m\n", nm[s], on, pk, rec[(s*nsteps + nsteps - 1)*4 + 1]);
      if(getenv("FD3D_STFILES")){                       /* time, slip, slip rate, shear stress (MPa) */
        sprintf(fn, "tpv101_%s.dat", nm[s]);
        fs = fopen(fn, "w");
        for(n = 0; n < nsteps; n++)
          fprintf(fs, "%.6e %.6e %.6e %.6e\n", rec[(s*nsteps + n)*4], rec[(s*nsteps + n)*4 + 1], rec[(s*nsteps + n)*4 + 2], rec[(s*nsteps + n)*4 + 3]/1e6);
        fclose(fs);
      }
    }
    {
      if(getenv("FD3D_PROF"))
        fprintf(stderr, "profile: %d steps; per step known part %.3f ms, solve %.3f ms, finish %.3f ms; per step %.1f Jacobi sweeps, %.0f local solves, %.2f Newton iterations each\n",
            (int)fd3d_prof3[3], 1e3*fd3d_prof3[0]/fd3d_prof3[3], 1e3*fd3d_prof3[1]/fd3d_prof3[3], 1e3*fd3d_prof3[2]/fd3d_prof3[3],
            fd3d_cnt[0]/fd3d_prof3[3], fd3d_cnt[1]/fd3d_prof3[3], (double)fd3d_cnt[2]/(fd3d_cnt[1] > 0 ? fd3d_cnt[1] : 1));
    }
    {
      fo = fopen("tpv101_rt.dat", "w");       /* x y rupture_time final_slip (km, s, m) */
      for(i = 0; i < nx; i++)
        for(j = 0; j < nz; j++)
          fprintf(fo, "%.4f %.4f %.5f %.5f\n", (xlo + (i + 0.5)*dx)/1e3, (ylo + (j + 0.5)*dx)/1e3, rt[i*nz + j], f->delta[i*nz + j]);
      fclose(fo);
    }
    i = (argc > 3 && strcmp(argv[3], "-")) ? compare(argv[3], 3, nm, rec, nsteps, dt*K, 1e-8, "tpv101") : 0;
  }
  free(rt);
  fd3d_free(f);
  fd3d_tab_free(t);
  free(rec);
  free(load);
  free(Fr);
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
  double x;
  double z;
  double tau;
  double a;
  double b;
  fd->nl = fm->nl = getenv("FD3D_NL") ? 1 : 0;
  fd->law = fm->law = 0;
  fd->mu_s = fm->mu_s = 0.677;
  fd->mu_d = fm->mu_d = 0.525;
  fd->d0 = fm->d0 = 0.40;
  for(i = 0; i < nx; i++)
    for(j = 0; j < nzp; j++){
      x = -10e3 + (i + 0.5)*dx;
      z = (j + 0.5)*dx;
      tau = (fabs(x) < 1.5e3 && z < 3e3) ? 81.6e6 : 70e6;
      fm->tau0[i*nzp + j] = tau;
      fm->sigma[i*nzp + j] = 120e6;
      fd->tau0[i*2*nzp + nzp + j] = tau;
      fd->tau0[i*2*nzp + nzp - 1 - j] = tau;
      fd->sigma[i*2*nzp + nzp + j] = 120e6;
      fd->sigma[i*2*nzp + nzp - 1 - j] = 120e6;
    }
  for(n = 0; n < nsteps; n++){
    fd3d_step(fm, NULL);
    fd3d_step(fd, NULL);
    for(i = 0; i < nx; i++)
      for(j = 0; j < nzp; j++){
        a = fm->V[i*nzp + j];
        b = fd->V[i*2*nzp + nzp + j];
        if(fabs(a - b) > emax)
          emax = fabs(a - b);
        if(fabs(b) > vmax)
          vmax = fabs(b);
        if(fabs(fm->delta[i*nzp + j] - fd->delta[i*2*nzp + nzp + j]) > ed)
          ed = fabs(fm->delta[i*nzp + j] - fd->delta[i*2*nzp + nzp + j]);
        if(fabs(fd->delta[i*2*nzp + nzp + j]) > dmax)
          dmax = fabs(fd->delta[i*2*nzp + nzp + j]);
      }
  }
  fail = !(emax <= 1e-8*vmax && ed <= 1e-8*dmax && dmax > 0);
  printf("fd3d_test mirror: %d x %d elements of %g m, %d steps: max |V diff| %.2e (of %.3f m/s), max |slip diff| %.2e (of %.3f m)\n",
      nx, nzp, dx, nsteps, emax, vmax, ed, dmax);
  printf("fd3d_test mirror: %s\n", fail ? "FAIL" : "PASS");
  fd3d_free(fm);
  fd3d_free(fd);
  fd3d_tab_free(t);
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
  extern double fd3d_layer_prof[4];
  f->nl = 1;
  f->law = 0;
  f->mu_s = f->mu_d = 0.0;
  f->d0 = 1.0;
  if(getenv("FD3D_LAYER_FFT"))
    f->layer = fd3d_layer_create_fft(f, -1000.0, 2, 16, 16, -2000.0, 45);
  else
    f->layer = fd3d_layer_create(f, -1000.0, 16, 16, 250.0, 250.0, -2000.0, -2000.0, 45);
  if(!f->layer || !fr){
    fprintf(stderr, "layer setup or reference failed\n");
    return 1;
  }
  vp = malloc(sizeof(double)*N);
  f->vpres = vp;
  for(n = 0; n < nsteps; n++){
    for(i = 0; i < N; i++)
      vp[i] = (f->time + dt <= 0.4 + 1e-12) ? 1.0 : 0.0;
    fd3d_step(f, NULL);
    for(i = 0; i < N; i++){
      if(fscanf(fr, "%lf", &r) != 1)
        break;
      if(fabs(f->tau[i] - r) > emax)
        emax = fabs(f->tau[i] - r);
      if(fabs(r) > smax)
        smax = fabs(r);
      cnt++;
    }
  }
  fclose(fr);
  {
    if(getenv("FD3D_PROF"))
      fprintf(stderr, "layer profile: known %.3f s, fault solves %.3f s, layer solves %.3f s, finish %.3f s\n",
          fd3d_layer_prof[0], fd3d_layer_prof[1], fd3d_layer_prof[2], fd3d_layer_prof[3]);
  }
  fail = !(cnt == nsteps*N && emax <= 1e-8*smax);
  printf("fd3d_test layer: %d values, max |tau diff| %.2e of %.3e\n", cnt, emax, smax);
  printf("fd3d_test layer: %s\n", fail ? "FAIL" : "PASS");
  fd3d_free(f);
  fd3d_tab_free(t);
  free(vp);
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
  f->nl = 1;
  f->law = 0;
  f->mu_s = f->mu_d = 0.0;
  f->d0 = 1.0;
  f->layer = fd3d_surf_create(f, -1000.0, 0.0, -500.0, dip, 2, 2, 24, (int)ny, -3000.0, ymid - ny*125.0, 45, 1.0);
  if(!f->layer || !fr){
    fprintf(stderr, "surface setup or reference failed\n");
    return 1;
  }
  vp = malloc(sizeof(double)*N);
  f->vpres = vp;
  for(n = 0; n < nsteps; n++){
    for(i = 0; i < N; i++)
      vp[i] = (f->time + dt <= 0.4 + 1e-12) ? 1.0 : 0.0;
    fd3d_step(f, NULL);
    for(i = 0; i < N; i++){
      if(fscanf(fr, "%lf", &r) != 1)
        break;
      if(fabs(f->tau[i] - r) > es)
        es = fabs(f->tau[i] - r);
      if(fabs(r) > ss)
        ss = fabs(r);
      cnt++;
    }
    for(i = 0; i < N; i++){
      if(fscanf(fr, "%lf", &r) != 1)
        break;
      if(fabs(f->layer->tn[i] - r) > en)
        en = fabs(f->layer->tn[i] - r);
      if(fabs(r) > sn)
        sn = fabs(r);
      cnt++;
    }
  }
  fclose(fr);
  fail = !(cnt == 2*nsteps*N && es <= 1e-8*ss && en <= 1e-8*sn);
  printf("fd3d_test surf: %d values, max |shear diff| %.2e of %.3e, max |normal diff| %.2e of %.3e\n", cnt, es, ss, en, sn);
  printf("fd3d_test surf: %s\n", fail ? "FAIL" : "PASS");
  fd3d_free(f);
  fd3d_tab_free(t);
  free(vp);
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
  ys0 = ytop - ceil((ytop - ys0)/hs)*hs;
  nys = (int)ceil((ytop + mg*dx - ys0)/hs);
  f->nl = 1;
  f->law = 0;
  f->mu_s = f->mu_d = 0.0;
  f->d0 = 1.0;
  f->vpres = vp;
  f->layer = fd3d_surf_create(f, -0.5*dx, ytop, ztop, dip, 2, m, nxs, nys, -0.5*dx - mg*dx, ys0, 50, 1.0);
  srand(1);
  for(i = 0; i < 3*f->layer->ns; i++){
    f->layer->V3[i] = rand()/(double)RAND_MAX - 0.5;
    if(fabs(f->layer->V3[i]) > v0)
      v0 = fabs(f->layer->V3[i]);
  }
  for(n = 0; n < 300; n++)
    fd3d_step(f, NULL);
  for(i = 0; i < 3*f->layer->ns; i++)
    if(fabs(f->layer->V3[i]) > v1 || !isfinite(f->layer->V3[i]))
      v1 = isfinite(f->layer->V3[i]) ? fabs(f->layer->V3[i]) : 1e300;
  fail = !(v1 <= 1e-3*v0);
  printf("fd3d_test surfstab: surface elements dx/%d, max rate %.2e -> %.2e after 300 steps\n", m, v0, v1);
  printf("fd3d_test surfstab: %s\n", fail ? "FAIL" : "PASS");
  fd3d_free(f);
  fd3d_tab_free(t);
  free(vp);
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
  size_t pq;
  ys0 = ytop - ceil((ytop - ys0)/hs)*hs;
  nys = (int)ceil((ytop + mg*dx - ys0)/hs);
  f->nl = 1;
  f->law = 0;
  f->mu_s = f->mu_d = 0.0;
  f->d0 = 1.0;
  ly = fd3d_surf_create(f, -0.5*dx, ytop, ztop, dip, 2, m, nxs, nys, -0.5*dx - mg*dx, ys0, L, 1.0);
  if(!ly || !fo){
    fprintf(stderr, "surfsym: setup failed\n");
    return 1;
  }
  fprintf(fo, "# cS dt/hs %.4f, dt %.6e, eta %.6e, L %d; columns: i j Ks*dt W(1..L)  (wavenumber i/8, j/8 of Nyquist)\n", m*cfl, dt, 0.5*2670.0*3464.0, L);
  for(i = 0; i <= 8; i++)
    for(j = 0; j <= 4; j++){
      pq = (size_t)((i*ly->Mx2/16) % ly->Mx2)*ly->Mh2 + (j*ly->My2/8) % ly->Mh2;
      fprintf(fo, "%d %d %.10e", i, j, creal(ly->T2[0][1][pq])*dt);
      for(l = 0; l < L; l++)
        fprintf(fo, " %.10e", creal(ly->T2[0][0][(size_t)l*ly->Mx2*ly->Mh2 + pq]));
      fprintf(fo, "\n");
    }
  fclose(fo);
  f->layer = ly;
  fd3d_free(f);
  fd3d_tab_free(t);
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
  struct fd3d_tab *t;
  struct fd3d *f;
  double *rt = malloc(sizeof(double)*N);
  FILE *fr = fopen(ref, "r");
  char line[512];
  double sl;
  double tr;
  int jj;
  ys0 = ytop - ceil((ytop - ys0)/hs)*hs;
  nys = (int)ceil((ytop + mg*dx - ys0)/hs);
  L = (int)(sqrt(pow((nx + 2*mg)*dx, 2) + pow(Wd + 2*mg*dx, 2))/(3464.0*dt)) + 4;
  t = fd3d_tab_create_c(nx, nz, dx, dt, 6000.0, 3464.0, 2670.0, L, 2);
  f = fd3d_create(t);
  f->nl = 1;
  f->law = 0;
  f->mu_s = 0.677;
  f->mu_d = 0.525;
  f->d0 = 0.40;
  for(i = 0; i < N; i++){
    f->sigma[i] = 120e6;
    f->tau0[i] = ((i % nz) >= nz - 2) ? 81.6e6 : 70e6;
    rt[i] = -1.0;
  }   /* thrust: + slip up dip */
  f->layer = fd3d_surf_create(f, -0.5*dx, ytop, ztop, dip, 2, m, nxs, nys, -0.5*dx - mg*dx, ys0, L, 1.0);
  if(!f->layer || !fr){
    fprintf(stderr, "q2d: setup or reference failed\n");
    return 1;
  }
  for(i = 0; i < N; i++)
    f->layer->sig0[i] = 120e6;
  for(n = 0; n < nsteps; n++){
    fd3d_step(f, NULL);
    for(i = 0; i < N; i++)
      if(rt[i] < 0 && fabs(f->V[i]) > 0.1)
        rt[i] = f->time;
  }
  while(fgets(line, sizeof(line), fr)){
    if(line[0] == '#' || sscanf(line, "%d %lf %lf", &jj, &sl, &tr) != 3)
      continue;
    i = (nx/2)*nz + jj;
    if(fabs(fabs(f->delta[i]) - sl)/sl > emax)
      emax = fabs(fabs(f->delta[i]) - sl)/sl;
    if(fabs(rt[i] - tr) > tmax)
      tmax = fabs(rt[i] - tr);
    cnt++;
  }
  fclose(fr);
  fail = !(cnt == nz && emax <= 0.02 && tmax <= 1.01*dt);
  printf("fd3d_test q2d: central column against the 2-D prototype: max relative slip difference %.2e, max rupture-time difference %.3f s (step %.3f s)\n", emax, tmax, dt);
  printf("fd3d_test q2d: %s\n", fail ? "FAIL" : "PASS");
  fd3d_free(f);
  fd3d_tab_free(t);
  free(rt);
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
  struct fd3d_tab *t;
  struct fd3d *f;
  FILE *fo[4];
  char fn[256];
  int e;
  double sdd;
  double xs;
  double sig;
  double tau;
  int nucl;
  double sdd_2;
  double jj;
  double w1;
  int j0;
  int j1;
  int ic[2];
  int k;
  int ncol;
  double v[4];
  int e0;
  int e1;
  /* surface rows: uniform hs over the projection plus mg fault elements either side, the trace on
     a row boundary */
  ys0 = ytop - Wd - mg*dx;
  nys = (int)((ytop + mg*dx - ys0)/hs + 0.5);
  /* window: the tables vanish once the S wave has crossed the farthest pair */
  L = (int)ceil(sqrt(pow((nx + 2*mg)*dx, 2) + pow(Wd + 2*mg*dx, 2) + pow(W*s, 2))/(cT*dt)) + 4;
  if(L > nsteps + 2)
    L = nsteps + 2;
  yc = malloc(sizeof(double)*nys);
  wy = malloc(sizeof(double)*nys);
  for(q = 0; q < nys; q++){
    yc[q] = ys0 + (q + 0.5)*hs;
    wy[q] = hs;
  }
  if(getenv("TPV10_NOSURF")){
    t = fd3d_tab_create_c(nx, nz, dx, dt, cL, cT, rho, L, 2);
    f = fd3d_create(t);
  } else{
    t = fd3d_tab_create_c(nx, nz + 1, dx, dt, cL, cT, rho, L, 2);
    f = fd3d_create3(t, nz, 0, 1);
  }
  f->nl = 1;
  f->law = 0;
  f->mu_s = mus;
  f->mu_d = mud;
  f->d0 = d0;
  f->coh = malloc(sizeof(double)*N);
  f->path = 1;
  f->dpath = calloc(N, sizeof(double));
  f->sw_saturate = 1;
  if(!getenv("TPV10_NOSURF")){
    f->layer = fd3d_surf_create_rows(f, -0.5*dx, ytop, ztop, dip, 2, m, nxs, nys, -0.5*dx - mg*dx, yc, wy, L, 1.0);
    if(!f->layer){
      fprintf(stderr, "tpv10: surface setup failed\n");
      return 1;
    }
    if(getenv("TPV10_NONC"))
      f->layer->normal_coupling = 0;
  }
  fprintf(stderr, "tpv10: %d x %d elements of %g m, %d steps of %.4f s, surface %d x %d of %g m, %d lags\n", nx, nz, dx, nsteps, dt, nxs, nys, hs, L);
  for(i = 0; i < nx; i++)
    for(j = 0; j < nz; j++){
      e = i*nz + j;
      sdd = (j + 0.5)*dx;
      xs = (i + 0.5)*dx - 0.5*nx*dx;
      sig = 7378.0*sdd;
      nucl = (fabs(sdd - 12000.0) <= 1500.0 + 1e-6) && (q2d || fabs(xs) <= 1500.0 + 1e-6);
      tau = nucl ? coh + (mus + 0.0057)*sig : 0.55*sig;
      f->tau0[e] = -tau;                 /* normal faulting: the hanging wall moves down dip (negative up-dip slip here) */
      f->sigma[e] = sig;
      if(f->layer)
        f->layer->sig0[e] = sig;
      f->coh[e] = coh;
    }
  for(q = 0; q < nst; q++){
    sprintf(fn, "%s_st%05.0f.dat", out, st_dd[q]);
    fo[q] = fopen(fn, "w");
  }
  for(n = 0; n < nsteps; n++){
    fd3d_step(f, NULL);
    for(q = 0; q < nst; q++){
      /* centre column (average of the two middle columns when nx is even), interpolated along dip */
      sdd_2 = st_dd[q];
      jj = sdd_2/dx - 0.5;
      v[0] = 0;
      v[1] = 0;
      v[2] = 0;
      v[3] = 0;
      if(jj < 0)
        jj = 0;
      if(jj > nz - 1)
        jj = nz - 1;
      j0 = (int)floor(jj);
      j1 = j0 + 1 < nz ? j0 + 1 : j0;
      w1 = jj - j0;
      if(nx % 2 == 0){
        ic[0] = nx/2 - 1;
        ic[1] = nx/2;
        ncol = 2;
      } else {
        ic[0] = nx/2;
        ncol = 1;
      }
      for(k = 0; k < ncol; k++){
        e0 = ic[k]*nz + j0;
        e1 = ic[k]*nz + j1;
        v[0] += ((1 - w1)*f->delta[e0] + w1*f->delta[e1])/ncol;
        v[1] += ((1 - w1)*f->V[e0] + w1*f->V[e1])/ncol;
        v[2] += ((1 - w1)*f->tau[e0] + w1*f->tau[e1])/ncol;
        v[3] += ((1 - w1)*f->sigma[e0] + w1*f->sigma[e1])/ncol;
      }
      fprintf(fo[q], "%.6e %.6e %.6e %.6e %.6e\n", f->time, -v[0], -v[1], -v[2]/1e6, -v[3]/1e6);
    }
  }
  for(q = 0; q < nst; q++)
    fclose(fo[q]);
  if(f->layer)
    fprintf(stderr, "tpv10: done; outer/inner iterations per step %.1f / %.1f\n", (double)f->layer->n_outer/f->layer->n_steps, (double)f->layer->n_inner/f->layer->n_steps);
  fd3d_free(f);
  fd3d_tab_free(t);
  free(yc);
  free(wy);
  return 0;
}

/* two parallel planes: plane 0 slips at a prescribed rate 1 for K steps (first-order scheme,
   for which the discrete convolution of a piecewise-constant rate is exact), plane 1 is at
   rest; the shear and normal traction on plane 1 must equal the direct sum of the exact ramp
   responses (fd3d_rect_ramp_tensor at the normal offset W) over the sources of plane 0 */
static int run_pair(int mirror)
{
  double cL = 6000, cT = 3464, rho = 2670, dx = 500.0, dt = 0.5*dx/cT, W = 1000.0, es = 0, ss = 0, en = 0, sn = 0, s6[6];
  int nx = 8, nz = 4, N = nx*nz, K = 4, nsteps = 40, n, i, r, q, fail;
  struct fd3d_tab *t = fd3d_tab_create(nx, mirror ? 2*nz : nz, dx, dt, cL, cT, rho, 0);
  struct fd3d *f0 = fd3d_create2(t, nz, mirror), *f1 = fd3d_create2(t, nz, mirror);
  struct fd3d_pair *pp;
  double *v0 = malloc(sizeof(double)*N), *v1 = calloc(N, sizeof(double));
  double ex_s;
  double ex_n;
  double x1;
  double sa[2];
  int k;
  int im;
  double x2;
  double tt;
  f0->nl = f1->nl = 0;
  f0->law = f1->law = 0;
  f0->mu_s = f0->mu_d = f1->mu_s = f1->mu_d = 0.0;
  f0->d0 = f1->d0 = 1.0;
  f0->vpres = v0;
  f1->vpres = v1;
  pp = fd3d_pair_create(f0, f1, W, 1);
  if(!pp){
    fprintf(stderr, "pair setup failed\n");
    return 1;
  }
  for(n = 1; n <= nsteps; n++){
    for(i = 0; i < N; i++)
      v0[i] = (n <= K) ? 1.0 : 0.0;
    fd3d_pair_step(pp, NULL, NULL);
    for(r = 0; r < N; r++){
      ex_s = 0.0;
      ex_n = 0.0;
      for(q = 0; q < N; q++){
        x1 = (r/nz - q/nz)*dx;
        sa[0] = 0;
        sa[1] = 0;

        for(im = 0; im <= mirror; im++){            /* the source and (mirror) its image above the surface */
          x2 = im ? (r%nz + q%nz + 1)*dx : (r%nz - q%nz)*dx;
          for(k = 0; k < 2; k++){                   /* ramp from 0 minus ramp from K dt */
            tt = n*dt - (k ? K*dt : 0.0);
            if(tt <= 0.0)
              continue;
            fd3d_rect_ramp_tensor(x1, x2, W, tt, 0.5*dx, 0.5*dx, cL, cT, rho, 1, s6);
            sa[0] += (k ? -1.0 : 1.0)*s6[4];
            sa[1] += (k ? -1.0 : 1.0)*s6[2];
          }
        }
        ex_s += sa[0];
        ex_n += sa[1];
      }
      if(fabs(f1->tau[r] - ex_s) > es)
        es = fabs(f1->tau[r] - ex_s);
      if(fabs(ex_s) > ss)
        ss = fabs(ex_s);
      if(fabs(pp->tn[1][r] - ex_n) > en)
        en = fabs(pp->tn[1][r] - ex_n);
      if(fabs(ex_n) > sn)
        sn = fabs(ex_n);
    }
  }
  fail = !(es <= 1e-8*ss && en <= 1e-8*sn && ss > 0 && sn > 0);
  printf("fd3d_test pair%s: plane 1 at %g m, %d steps: max |shear diff| %.2e of %.3e, max |normal diff| %.2e of %.3e (k1 cross %.1e)\n", mirror ? " (mirror)" : "", W, nsteps, es, ss, en, sn, pp->k1c);
  printf("fd3d_test pair%s: %s\n", mirror ? " (mirror)" : "", fail ? "FAIL" : "PASS");
  fd3d_pair_free(pp);
  fd3d_free(f0);
  fd3d_free(f1);
  fd3d_tab_free(t);
  free(v0);
  free(v1);
  return fail;
}

/* 3-D stepover (fd_plan.md P1): two vertical strike-slip segments of length Lseg and depth D
   in a whole space, overlapping by `ov` along strike, plane 1 at normal offset W (sign: see
   the printed normal-stress change), slip weakening (mu_s 0.75, mu_d 0.51, d0, sigma 33.3
   MPa, tau0 20 MPa: Harris and Day's case D strength parameter S = 1.65) with normal-stress
   coupling, nucleation by an overstressed 6 x 6 km patch at the centre of segment 0; reports
   whether and when segment 1 ruptures */
static int run_stepover(double W, double dx, double tend, double d0, double mus)
{
  double cL = 6000, cT = 3464, rho = 2700, dt = 0.5*dx/cT, Lseg = 28e3, D = 12e3, ov = 5e3, mud = 0.51, sig = 33.3e6, tau0 = 20e6;
  int n1 = (int)(Lseg/dx + 0.5), nov = (int)(ov/dx + 0.5), nx = 2*n1 - nov, nz = (int)(D/dx + 0.5), N = nx*nz, i, j, n, p, nsteps = (int)(tend/dt + 0.5);
  struct fd3d_tab *t = fd3d_tab_create(nx, nz, dx, dt, cL, cT, rho, 0);
  struct fd3d *f[2];
  struct fd3d_pair *pp;
  double trig[2] = {-1, -1}, xtrig = 0, smin = 1e300, smax = -1e300;
  int e;
  int on;
  double xc;
  double zc;
  int len[2];
  double smx[2];
  int any;
  int e_2;
  for(p = 0; p < 2; p++){
    f[p] = fd3d_create(t);
    f[p]->nl = 1;
    f[p]->law = 0;
    f[p]->mu_s = mus;
    f[p]->mu_d = mud;
    f[p]->d0 = d0;
    f[p]->locked = calloc(N, sizeof(int));
    f[p]->path = 1;
    f[p]->dpath = calloc(N, sizeof(double));
    f[p]->sw_saturate = 1;
    for(i = 0; i < nx; i++)
      for(j = 0; j < nz; j++){
        e = i*nz + j;
        on = p == 0 ? (i < n1) : (i >= n1 - nov);
        xc = (i - 0.5*n1 + 0.5)*dx;
        zc = (j + 0.5)*dx - 0.5*D;    /* relative to segment 0's centre */
        f[p]->locked[e] = !on;
        f[p]->sigma[e] = sig;
        f[p]->tau0[e] = (p == 0 && fabs(xc) <= 3000.0 && fabs(zc) <= 3000.0) ? (mus + 0.005)*sig : tau0;   /* 6 x 6 km patch */
      }
  }
  pp = fd3d_pair_create(f[0], f[1], W, 1);
  if(!pp){
    fprintf(stderr, "stepover: pair setup failed (W too small for dx?)\n");
    return 1;
  }
  pp->normal_coupling = 1;
  for(p = 0; p < 2; p++)
    for(i = 0; i < N; i++)
      pp->sig0[p][i] = sig;
  fprintf(stderr, "stepover: W %+g m, dx %g m, %d x %d lattice (segments %d columns, overlap %d), %d steps, window %d lags\n", W, dx, nx, nz, n1, nov, nsteps, t->L);
  for(n = 0; n < nsteps; n++){
    fd3d_pair_step(pp, NULL, NULL);
    for(p = 0; p < 2; p++)
      if(trig[p] < 0)
        for(i = 0; i < N; i++)
          if(!f[p]->locked[i] && fabs(f[p]->V[i]) > 1e-3){
            trig[p] = f[p]->time;
            if(p == 1)
              xtrig = (i/nz - 0.5*n1 + 0.5)*dx;
            break;
          }
    for(i = 0; i < N; i++)
      if(!f[1]->locked[i]){
        if(f[1]->sigma[i] < smin)
          smin = f[1]->sigma[i];
        if(f[1]->sigma[i] > smax)
          smax = f[1]->sigma[i];
      }
  }
  {
    len[0] = 0;
    len[1] = 0;
    smx[0] = 0;
    smx[1] = 0;
    for(p = 0; p < 2; p++)
      for(i = 0; i < nx; i++){
        any = 0;
        for(j = 0; j < nz; j++){
          e_2 = i*nz + j;
          if(!f[p]->locked[e_2] && fabs(f[p]->delta[e_2]) > 0.05)
            any = 1;
          if(fabs(f[p]->delta[e_2]) > smx[p])
            smx[p] = fabs(f[p]->delta[e_2]);
        }
        len[p] += any;
      }
    printf("stepover mu_s %.3f (S %.2f) W %+6.0f m: segment 0 ruptured %.1f km (max slip %.2f m); segment 1 %s at %.2f s, x = %.1f km, ruptured %.1f km (max slip %.2f m); sigma on segment 1 %.1f .. %.1f MPa\n",
        mus, (mus*sig - tau0)/(tau0 - mud*sig), W, len[0]*dx/1e3, smx[0], trig[1] >= 0 ? "triggered" : "not triggered", trig[1], xtrig/1e3, len[1]*dx/1e3, smx[1], smin/1e6, smax/1e6);
  }
  fd3d_pair_free(pp);
  fd3d_free(f[0]);
  fd3d_free(f[1]);
  fd3d_tab_free(t);
  return 0;
}

/* several planes (fd_plan.md G1): a vertical strike-slip fault of nx x nz cells of 500 m with
   slip weakening and an overstressed patch, run (a) as one plane and (b) split into two
   planes at column nx/2 whose second half is turned by `ang` degrees about the vertical
   through the seam (0: coplanar, the split must reproduce (a) to the solver tolerance), or
   tilted (`dip` < 90: the second half dips). Prints the largest differences or the result. */
static void run_sw_setup(struct fd3d *f, int nzp, double i0, double dx, int nxall, int off)
{
  int i, N = f->N;
  int ii;
  int jj;
  double xc;
  double zc;
  (void)nxall;
  f->law = 0;
  f->mu_s = 0.677;
  f->mu_d = 0.525;
  f->d0 = 0.40;
  f->nl = 1;
  for(i = 0; i < N; i++){
    ii = i/nzp + off;
    jj = i%nzp;
    xc = (ii + 0.5)*dx - i0;
    zc = (jj + 0.5)*dx - 0.5*nzp*dx;
    f->sigma[i] = 120e6;
    f->tau0[i] = (fabs(xc) <= 1500.0 && fabs(zc) <= 1500.0) ? 81.6e6 : 70e6;
  }
}
/* TPV14 / TPV15 (SCEC branched-fault benchmarks, Barall 2011): vertical strike-slip main fault
   from x = -16 to +12 km, 15 km deep; vertical branch from the junction at x = 0, 12 km long,
   turned 30 degrees towards -y (the "near side"); half space by free-surface images (G3a);
   linear slip weakening (mu_s 0.677, mu_d 0.525, d0 0.4 m), normal stress 120 MPa with normal
   stress changes from the other fault; 3 x 3 km nucleation at x = -8 km, 7.5 km depth.
   TPV14: right-lateral, 70 MPa on both faults; TPV15: left-lateral, 70 MPa on the main fault and
   78 MPa on the branch. Frames: main e1 = x, e2 = z, e3 = -y; branch e1 = (cos 30, -sin 30, 0),
   e2 = z, e3 = e1 x e2; slip along e1 of the +e3 side, so right-lateral is negative here.
   gap: distance in m from the junction to the start of the branch (0: touching; the benchmark's
   split-node rule corresponds to one 100 m element; earlier versions took 1 for one cell).
   Writes tpv14_rt_main.dat and tpv14_rt_branch.dat (x along strike, depth, rupture time,
   final slip) and prints the ruptured fractions. */
static void tpv14_frame(double az_deg, double E[3][3])
{
  double a;
  a = az_deg*M_PI/180.0;
  E[0][0] = cos(a);
  E[0][1] = sin(a);
  E[0][2] = 0.0;
  E[1][0] = 0.0;
  E[1][1] = 0.0;
  E[1][2] = 1.0;
  E[2][0] = E[0][1]*E[1][2] - E[0][2]*E[1][1];
  E[2][1] = E[0][2]*E[1][0] - E[0][0]*E[1][2];
  E[2][2] = E[0][0]*E[1][1] - E[0][1]*E[1][0];
}
static void tpv14_friction(struct fd3d *f, int nx, int nz, double dx, double tau, double tnuc, double x0)
{
  int i, j;
  double xc, zc;
  f->law = 0;
  f->mu_s = 0.677;
  f->mu_d = 0.525;
  f->d0 = 0.40;
  f->nl = 1;
  for(i = 0; i < nx; i++)
    for(j = 0; j < nz; j++){
      xc = x0 + (i + 0.5)*dx;
      zc = (j + 0.5)*dx;
      f->sigma[i*nz + j] = 120e6;
      f->tau0[i*nz + j] = tau;
      if(tnuc != 0.0 && fabs(xc + 8000.0) <= 1500.0 && fabs(zc - 7500.0) <= 1500.0)
        f->tau0[i*nz + j] = tnuc;
    }
}
static int run_tpv14(double dx, int which, double tend, double gap)
{
  double cL, cT, rho, dt, sgn, tau_b, x0b, tnmax, vthr;
  double *rt[2], *slip[2];
  int nxm, nxb, nz, ns, n, p, i, j, nx[2], nrup[2], nrup_beyond;
  struct fd3d_tab *tm, *tb;
  struct fd3d *fa[2];
  struct fd3d_mplane g[2];
  struct fd3d_multi *m;
  FILE *out;
  char fn[64];

  cL = 6000.0;
  cT = 3464.0;
  rho = 2670.0;
  dt = 0.5*dx/cT;
  ns = (int)(tend/dt + 0.5);
  nz = (int)(15000.0/dx + 0.5);
  nxm = (int)(28000.0/dx + 0.5);
  x0b = gap;
  nxb = (int)((12000.0 - x0b)/dx + 0.5);
  nx[0] = nxm;
  nx[1] = nxb;
  sgn = (which == 15) ? 1.0 : -1.0;              /* right-lateral (TPV14) is negative slip here */
  tau_b = (which == 15) ? 78e6 : 70e6;
  vthr = 1e-3;
  tm = fd3d_tab_create_c(nxm, 2*nz, dx, dt, cL, cT, rho, 0, 1);
  tb = fd3d_tab_create_c(nxb, 2*nz, dx, dt, cL, cT, rho, 0, 1);
  fa[0] = fd3d_create2(tm, nz, 1);
  fa[1] = fd3d_create2(tb, nz, 1);
  tpv14_friction(fa[0], nxm, nz, dx, sgn*70e6, sgn*81.6e6, -16000.0);
  tpv14_friction(fa[1], nxb, nz, dx, sgn*tau_b, 0.0, 0.0);
  memset(g, 0, sizeof(g));
  for(p = 0; p < 2; p++){
    g[p].dx = dx;
    g[p].comp = 1;
    g[p].mirror = 1;
  }
  tpv14_frame(0.0, g[0].E);
  g[0].O[0] = -16000.0;
  tpv14_frame(-30.0, g[1].E);
  for(i = 0; i < 3; i++)
    g[1].O[i] = x0b*g[1].E[0][i];
  m = fd3d_multi_create(2, fa, g);
  if(!m){
    fprintf(stderr, "tpv14: setup failed\n");
    return 1;
  }
  m->normal_coupling = 1;
  for(p = 0; p < 2; p++){
    rt[p] = malloc(sizeof(double)*fa[p]->N);
    slip[p] = calloc(fa[p]->N, sizeof(double));
    for(i = 0; i < fa[p]->N; i++){
      m->sig0[p][i] = 120e6;
      rt[p][i] = -1.0;
    }
  }
  fprintf(stderr, "tpv%d: dx %g m, main %d x %d, branch %d x %d cells, %d steps of %.4f s, cross windows %d / %d lags\n",
          which, dx, nxm, nz, nxb, nz, ns, dt, m->c[0][1].L, m->c[1][0].L);
  tnmax = 0.0;
  for(n = 0; n < ns; n++){
    fd3d_multi_step(m, NULL);
    for(p = 0; p < 2; p++)
      for(i = 0; i < fa[p]->N; i++){
        slip[p][i] += fa[p]->V[i]*dt;
        if(rt[p][i] < 0.0 && fabs(fa[p]->V[i]) > vthr)
          rt[p][i] = (n + 1)*dt;
      }
    /* normal traction change on the branch's first column before it slips (unclamping > 0) */
    for(j = 0; j < nz; j++)
      if(rt[1][j] < 0.0 && m->tn[1][j] > tnmax)
        tnmax = m->tn[1][j];
  }
  for(p = 0; p < 2; p++){
    nrup[p] = 0;
    for(i = 0; i < fa[p]->N; i++)
      if(rt[p][i] >= 0.0)
        nrup[p]++;
  }
  nrup_beyond = 0;
  for(i = 0; i < fa[0]->N; i++)
    if(rt[0][i] >= 0.0 && (-16000.0 + (i/nz + 0.5)*dx) > 0.0)
      nrup_beyond++;
  printf("tpv%d: dx %g m, gap %g m, %g s: main fault ruptured %.0f %% (beyond the junction %.0f %%), branch %.0f %%; largest unclamping of the branch's first column before it slips %.2f MPa\n",
         which, dx, gap, tend, 100.0*nrup[0]/fa[0]->N, 100.0*nrup_beyond/(double)((nxm - (int)(16000.0/dx + 0.5))*nz),
         100.0*nrup[1]/fa[1]->N, tnmax/1e6);
  if(getenv("FD3D_DUMP")){                        /* exact state for byte comparisons (MPI, fd_plan.md H7) */
    fd3d_dump_state(getenv("FD3D_DUMP"), "w", 2, fa);
    for(p = 0; p < 2; p++){
      fd3d_dump_vec(getenv("FD3D_DUMP"), p ? "rt1" : "rt0", fa[p]->N, rt[p]);
      fd3d_dump_vec(getenv("FD3D_DUMP"), p ? "tn1" : "tn0", fa[p]->N, m->tn[p]);
    }
  }
  for(p = 0; p < 2; p++){
    sprintf(fn, "tpv%d_rt_%s.dat", which, p ? "branch" : "main");
    out = fopen(fn, "w");
    for(i = 0; i < nx[p]; i++)
      for(j = 0; j < nz; j++)
        fprintf(out, "%g %g %g %g\n", (p ? x0b : -16000.0) + (i + 0.5)*dx, (j + 0.5)*dx, rt[p][i*nz + j], slip[p][i*nz + j]);
    fclose(out);
    free(rt[p]);
    free(slip[p]);
  }
  fd3d_multi_free(m);
  for(p = 0; p < 2; p++)
    fd3d_free(fa[p]);
  fd3d_tab_free(tm);
  fd3d_tab_free(tb);
  return 0;
}

static int run_multi(double ang, double dip, double tend, int comp, int mir)
{
  double cL = 6000, cT = 3464, rho = 2670, dx = 500.0, dt = 0.5*dx/cT, i0 = 6.5*dx, emax = 0, smax = 0;
  int nx = 24, nz = 12, h = nx/2, n, i, ns = (int)(tend/dt + 0.5), p;
  struct fd3d_tab *t1 = fd3d_tab_create_c(nx, mir ? 2*nz : nz, dx, dt, cL, cT, rho, 0, comp), *th = fd3d_tab_create_c(h, mir ? 2*nz : nz, dx, dt, cL, cT, rho, 0, comp);
  struct fd3d *f1 = fd3d_create2(t1, nz, mir), *fa[2];
  struct fd3d_mplane g[2];
  struct fd3d_multi *m;
  double a = ang*M_PI/180.0, dp = dip*M_PI/180.0;
  double hz[3];
  double dd;
  double e1[3];
  double hz_2[3];
  double e2[3];
  int k;
  int ii;
  int jj;
  int k_2;
  double d;
  double s1;
  double s2;
  int ok;
  run_sw_setup(f1, nz, i0, dx, nx, 0);
  for(p = 0; p < 2; p++){
    fa[p] = fd3d_create2(th, nz, mir);
    run_sw_setup(fa[p], nz, i0, dx, nx, p ? h : 0);
  }
  memset(g, 0, sizeof(g));
  g[0].dx = g[1].dx = dx;
  g[0].comp = g[1].comp = comp;
  g[0].mirror = g[1].mirror = mir;
  g[0].E[0][0] = 1.0;
  g[0].E[1][2] = 1.0;
  g[0].E[2][1] = -1.0;              /* e1 = x, e2 = z up, e3 = -y */
  g[1].O[0] = h*dx;
  g[1].E[0][0] = cos(a);
  g[1].E[0][1] = sin(a);                               /* strike turned about z */
  {                                                                          /* up dip: sin(dip) z + cos(dip) (horizontal, normal to e1) */
    hz[0] = sin(a);
    hz[1] = -cos(a);
    hz[2] = 0.0;
    for(i = 0; i < 3; i++)
      g[1].E[1][i] = cos(dp)*hz[i] + sin(dp)*(i == 2);
  }
  g[1].E[2][0] = g[1].E[0][1]*g[1].E[1][2] - g[1].E[0][2]*g[1].E[1][1];
  g[1].E[2][1] = g[1].E[0][2]*g[1].E[1][0] - g[1].E[0][0]*g[1].E[1][2];
  g[1].E[2][2] = g[1].E[0][0]*g[1].E[1][1] - g[1].E[0][1]*g[1].E[1][0];
  if(dip < 0.0){
    /* rigid rotation of the whole split fault (coplanar halves, strike turned by ang, dip |dip|):
       in a whole space it must reproduce the single plane exactly */
    dd = -dip*M_PI/180.0;
    e1[0] = cos(a);
    e1[1] = sin(a);
    e1[2] = 0.0;
    hz_2[0] = sin(a);
    hz_2[1] = -cos(a);
    hz_2[2] = 0.0;
    for(i = 0; i < 3; i++)
      e2[i] = cos(dd)*hz_2[i] + sin(dd)*(i == 2);
    for(p = 0; p < 2; p++){
      for(k = 0; k < 3; k++){
        g[p].E[0][k] = e1[k];
        g[p].E[1][k] = e2[k];
      }
      g[p].E[2][0] = e1[1]*e2[2] - e1[2]*e2[1];
      g[p].E[2][1] = e1[2]*e2[0] - e1[0]*e2[2];
      g[p].E[2][2] = e1[0]*e2[1] - e1[1]*e2[0];
      for(k = 0; k < 3; k++)
        g[p].O[k] = (p ? h*dx : 0.0)*e1[k];
    }
  }
  m = fd3d_multi_create(2, fa, g);
  if(!m){
    fprintf(stderr, "multi setup failed\n");
    return 1;
  }
  fprintf(stderr, "multi: turn %g deg, dip %g deg, cross first lag %s, cross windows %d / %d lags\n", ang, dip, m->implicit ? "kept (implicit)" : "zero", m->c[0][1].L, m->c[1][0].L);
  for(n = 0; n < ns; n++){
    fd3d_step(f1, NULL);
    fd3d_multi_step(m, NULL);
    for(p = 0; p < 2; p++)
      for(i = 0; i < fa[p]->N; i++){
        ii = i/nz + (p ? h : 0);
        jj = i%nz;
        k_2 = ii*nz + jj;
        d = fabs(fa[p]->V[i] - f1->V[k_2]);
        if(d > emax)
          emax = d;
        if(fabs(f1->V[k_2]) > smax)
          smax = fabs(f1->V[k_2]);
      }
  }
  {
    s1 = 0;
    s2 = 0;
    for(i = 0; i < f1->N; i++)
      s1 += f1->delta[i];
    for(p = 0; p < 2; p++)
      for(i = 0; i < fa[p]->N; i++)
        s2 += fa[p]->delta[i];
    printf("fd3d_test multi%s: comp %d, turn %g deg, dip %g deg, %d steps: max |V(two planes) - V(one plane)| %.3e of %.3f m/s; mean slip %.5f (one plane) %.5f (two planes) m; fixed-point iterations last step %d\n",
        mir ? " (half space, images)" : "", comp, ang, dip, ns, emax, smax, s1/f1->N, s2/(2*fa[0]->N), m->last_iters);
    if(getenv("FD3D_DUMP")){                      /* exact state for byte comparisons (MPI, fd_plan.md H7) */
      fd3d_dump_state(getenv("FD3D_DUMP"), "w", 2, fa);
      for(p = 0; p < 2; p++)
        fd3d_dump_vec(getenv("FD3D_DUMP"), p ? "tn1" : "tn0", fa[p]->N, m->tn[p]);
    }
    if((ang == 0.0 && dip >= 90.0) || dip < 0.0){
      ok = emax <= (dip < 0.0 ? 1e-7 : 1e-8)*smax;   /* rotated: round-off of the rotations, amplified by the rupture */
      printf("fd3d_test multi: %s\n", ok ? "PASS" : "FAIL");
      return !ok;
    }
  }
  return 0;
}

/* dense several-plane coupling against the Toeplitz pair for two parallel planes 1 km apart:
   plane 0 slips at a prescribed rate (first order), traction on plane 1 from both objects.
   The pair's frame has e3 = e1 x (down) = -e3 of the multi frame, so the same receivers sit at
   pair offset -W, and its slip and normal-traction signs are opposite: expected tau equal,
   normal traction opposite. */
static int run_multi_pair(void)
{
  double cL = 6000, cT = 3464, rho = 2670, dx = 500.0, dt = 0.5*dx/cT, W = 1000.0, es = 0, ss = 0, en = 0, sn = 0;
  int nx = 8, nz = 4, N = nx*nz, n, i, ok;
  struct fd3d_tab *t = fd3d_tab_create(nx, nz, dx, dt, cL, cT, rho, 0);
  struct fd3d *a0 = fd3d_create(t), *a1 = fd3d_create(t), *b[2];
  struct fd3d_pair *pp;
  struct fd3d_multi *m;
  struct fd3d_mplane g[2];
  double *v0 = malloc(sizeof(double)*N), *vz = calloc(N, sizeof(double));
  b[0] = fd3d_create(t);
  b[1] = fd3d_create(t);
  a0->nl = a1->nl = b[0]->nl = b[1]->nl = 0;
  a0->vpres = v0;
  a1->vpres = vz;
  b[0]->vpres = v0;
  b[1]->vpres = vz;
  pp = fd3d_pair_create(a0, a1, -W, 1);
  memset(g, 0, sizeof(g));
  for(i = 0; i < 2; i++){
    g[i].dx = dx;
    g[i].comp = 1;
    g[i].E[0][0] = 1.0;
    g[i].E[1][2] = 1.0;
    g[i].E[2][1] = -1.0;
  }
  g[1].O[1] = -W;                                                           /* +W along e3 = -y */
  m = fd3d_multi_create(2, b, g);
  if(!pp || !m){
    fprintf(stderr, "setup failed\n");
    return 1;
  }
  for(n = 1; n <= 40; n++){
    for(i = 0; i < N; i++)
      v0[i] = (n <= 4) ? 1.0 : 0.0;
    fd3d_pair_step(pp, NULL, NULL);
    fd3d_multi_step(m, NULL);
    for(i = 0; i < N; i++){
      if(fabs(b[1]->tau[i] - a1->tau[i]) > es)
        es = fabs(b[1]->tau[i] - a1->tau[i]);
      if(fabs(a1->tau[i]) > ss)
        ss = fabs(a1->tau[i]);
      if(fabs(m->tn[1][i] + pp->tn[1][i]) > en)
        en = fabs(m->tn[1][i] + pp->tn[1][i]);
      if(fabs(pp->tn[1][i]) > sn)
        sn = fabs(pp->tn[1][i]);
    }
  }
  ok = es <= 1e-8*ss && en <= 1e-8*sn;
  printf("fd3d_test multipair: 40 steps: max |tau(multi) - tau(pair)| %.3e of %.3e, max |tn(multi) + tn(pair)| %.3e of %.3e\n", es, ss, en, sn);
  printf("fd3d_test multipair: %s\n", ok ? "PASS" : "FAIL");
  return !ok;
}

/* several faults with one virtual surface (fd_plan.md G3b): a buried dipping fault (dip 60
   degrees, 16 x 8 cells of 500 m, top 1 km deep, slip weakening, overstressed patch) with a
   graded-free uniform surface, run (a) as one fault with fd3d_surf_create and (b) split along
   dip into an upper and a lower fault sharing the surface (fd3d_msurf); (b) must reproduce (a). */
static void msurf_sw(struct fd3d *f, int nzp, int j0, double dx)
{
  int i;
  int ii;
  int jj;
  double xc;
  double zc;
  f->law = 0;
  f->mu_s = 0.677;
  f->mu_d = 0.525;
  f->d0 = 0.40;
  f->nl = 1;
  for(i = 0; i < f->N; i++){
    ii = i/nzp;
    jj = i%nzp + j0;
    xc = (ii + 0.5)*dx - 4000.0;
    zc = (jj + 0.5)*dx - 2000.0;
    f->sigma[i] = 120e6;
    f->tau0[i] = (fabs(xc) <= 1500.0 && fabs(zc) <= 1500.0) ? 81.6e6 : 70e6;
  }
}
static int run_msurf(int comp, double tend)
{
  double cL = 6000, cT = 3464, rho = 2670, dx = 500.0, dt = 0.5*dx/cT, dip = 60.0, d = dip*M_PI/180.0, emax = 0, smax = 0;
  double ytop[2], ztop[2], dips[2] = {dip, dip}, fs[2] = {1.0, 1.0};
  int nx = 16, nz = 8, h = 4, m = 2, nxs = 48, nys = 24, n, i, k, ns = (int)(tend/dt + 0.5), L;
  struct fd3d_tab *t1 = fd3d_tab_create_c(nx, nz, dx, dt, cL, cT, rho, 0, comp), *th = fd3d_tab_create_c(nx, h, dx, dt, cL, cT, rho, 0, comp);
  struct fd3d *f1 = fd3d_create(t1), *fa[2];
  struct fd3d_msurf *ms;
  double *yc = malloc(sizeof(double)*nys), *wy = malloc(sizeof(double)*nys);
  int ii;
  int jj;
  int q;
  double dd;
  double s1;
  double s2;
  int ok;
  L = t1->L;
  for(i = 0; i < nys; i++){
    yc[i] = -4000.0 + (i + 0.5)*dx/m;
    wy[i] = dx/m;
  }
  msurf_sw(f1, nz, 0, dx);
  f1->layer = fd3d_surf_create_rows(f1, 0.0, 0.0, -1000.0, dip, comp, m, nxs, nys, -2000.0, yc, wy, L, 1.0);
  if(!f1->layer){
    fprintf(stderr, "surf failed\n");
    return 1;
  }
  f1->layer->normal_coupling = 1;
  if(getenv("MSURF_ONE")){                         /* debugging: the same single fault through fd3d_msurf */
    h = nz;
    th = fd3d_tab_create_c(nx, h, dx, dt, cL, cT, rho, 0, comp);
    fa[0] = fd3d_create(th);
    msurf_sw(fa[0], h, 0, dx);
    ytop[0] = 0.0;
    ztop[0] = -1000.0;
    ms = fd3d_msurf_create(1, fa, ytop, ztop, dips, fs, comp, 0.0, m, nxs, nys, -2000.0, yc, wy, L);
  }else{
    for(k = 0; k < 2; k++){
      fa[k] = fd3d_create(th);
      msurf_sw(fa[k], h, k*h, dx);
    }
    ytop[0] = 0.0;
    ztop[0] = -1000.0;
    ytop[1] = -h*dx*cos(d);
    ztop[1] = -1000.0 - h*dx*sin(d);
    if(getenv("MSURF_DIP2"))
      dips[1] = atof(getenv("MSURF_DIP2"));     /* a kink in dip at the joint (no reference: runs) */
    ms = fd3d_msurf_create(2, fa, ytop, ztop, dips, fs, comp, 0.0, m, nxs, nys, -2000.0, yc, wy, L);
  }
  if(!ms){
    fprintf(stderr, "msurf failed\n");
    return 1;
  }
  fprintf(stderr, "msurf: comp %d, cross first lag %s\n", comp, ms->mx ? (ms->mx->implicit ? "kept" : "zero") : "none");
  for(n = 0; n < ns; n++){
    fd3d_step(f1, NULL);
    fd3d_msurf_step(ms, NULL);
    for(k = 0; k < ms->nf; k++)
      for(i = 0; i < fa[k]->N; i++){
        ii = i/h;
        jj = i%h + k*h;
        q = ii*nz + jj;
        dd = fabs(fa[k]->V[i] - f1->V[q]);
        if(dd > emax)
          emax = dd;
        if(fabs(f1->V[q]) > smax)
          smax = fabs(f1->V[q]);
      }
  }
  {
    s1 = 0;
    s2 = 0;
    ok = emax <= 1e-7*smax;
    for(i = 0; i < f1->N; i++)
      s1 += f1->delta[i];
    for(k = 0; k < ms->nf; k++)
      for(i = 0; i < fa[k]->N; i++)
        s2 += fa[k]->delta[i];
    printf("fd3d_test msurf: comp %d, dip %g, %d steps: max |V(two faults) - V(one fault)| %.3e of %.3f m/s; mean slip %.5f / %.5f m; outer iterations last step %d\n",
        comp, dip, ns, emax, smax, s1/f1->N, s2/f1->N, ms->last_iters);
    printf("fd3d_test msurf: %s\n", ok ? "PASS" : "FAIL");
    if(getenv("FD3D_DUMP")){                      /* exact state for byte comparisons (MPI, fd_plan.md H7) */
      fd3d_dump_state(getenv("FD3D_DUMP"), "w", ms->nf, fa);
      fd3d_dump_vec(getenv("FD3D_DUMP"), "surf", 3*ms->ly[0]->ns, ms->ly[0]->V3);
    }
    return !ok;
  }
}

/* kernel check: full stress tensors (slip and opening, rectangles) against tada_quad.py values */
static int run_kernels(const char *ref)
{
  FILE *fr = fopen(ref, "r");
  double v[13], s[6], emax = 0.0, smax = 0.0;
  int n = 0, k;
  double sc;
  double e;
  if(!fr){
    fprintf(stderr, "cannot open %s\n", ref);
    return 1;
  }
  while(fscanf(fr, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", v, v + 1, v + 2, v + 3, v + 4, v + 5, v + 6,
            v + 7, v + 8, v + 9, v + 10, v + 11, v + 12) == 13){
    sc = 1.0;
    fd3d_rect_ramp_tensor(v[1], v[2], v[3], v[4], v[5], v[6], 6000.0, 3464.0, 2670.0, (int)v[0], s);
    for(k = 0; k < 6; k++)
      if(fabs(v[7 + k]) > sc)
        sc = fabs(v[7 + k]);
    for(k = 0; k < 6; k++){
      e = fabs(s[k] - v[7 + k])/sc;          /* relative to the point's largest component */
      if(!(e <= emax))
        emax = isfinite(e) ? e : 1e300;
    }
    if(sc > smax)
      smax = sc;
    n++;
  }
  fclose(fr);
  printf("fd3d_test kernels: %d points (slip along x1 and x2, opening; rectangles, edge cases), max relative diff %.2e\n", n, emax);
  /* 1e-7: the closed forms cancel large terms at some points; Python and C (with fused
     multiply-add) differ at the 1e-8 level there */
  printf("fd3d_test kernels: %s\n", (n > 0 && emax <= 1e-7) ? "PASS" : "FAIL");
  return !(n > 0 && emax <= 1e-7);
}

static int main0(int argc, char **argv);
/* H6 step 4: where the largest differences between two resolutions are: per fault and traction
   component (shear, normal), the cell (along strike i, down dip j; 2 rows) and the step */
static void h6strike_where(const double *t1, const double *t2, int ns, int N)
{
  double best[2][2], d;
  int best_i[2][2], best_n[2][2];
  int n, k, i, c;
  for(k = 0; k < 2; k++)
    for(c = 0; c < 2; c++){
      best[k][c] = 0.0;
      best_i[k][c] = 0;
      best_n[k][c] = 0;
    }
  for(n = 0; n < ns; n++)
    for(k = 0; k < 2; k++)
      for(i = 0; i < 2*N; i++){
        c = (i >= N);
        d = fabs(t1[((size_t)n*2 + k)*2*N + i] - t2[((size_t)n*2 + k)*2*N + i]);
        if(d > best[k][c]){
          best[k][c] = d;
          best_i[k][c] = i % N;
          best_n[k][c] = n;
        }
      }
  for(k = 0; k < 2; k++)
    printf("h6strike: fault %d: largest shear difference %.3e at cell (%d, %d), step %d; normal %.3e at cell (%d, %d), step %d\n",
           k, best[k][0], best_i[k][0]/2, best_i[k][0]%2, best_n[k][0], best[k][1], best_i[k][1]/2, best_i[k][1]%2, best_n[k][1]);
}

int main(int argc, char **argv)
{
  int r;
#ifdef FD3D_MPI
  /* under mpirun: the tables and histories are split over the ranks (fd_plan.md S2); rank 0
     prints and writes, the other ranks work in a scratch directory with stdout discarded */
  int rank = 0, size = 1;
  char d[64];
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if(size > 1){
    fd3d_set_comm(MPI_COMM_WORLD);
    if(rank > 0){
      sprintf(d, "/tmp/fd3d_test_rank%d_XXXXXX", rank);
      if(mkdtemp(d) && chdir(d) == 0){
      }
      if(!freopen("/dev/null", "w", stdout)){
      }
    }
  }
#endif
  r = main0(argc, argv);
#ifdef FD3D_MPI
  MPI_Finalize();
#endif
  return r;
}
static int main0(int argc, char **argv)
{
  double cL = 6000, cT = 3464, rho = 2670, dx = argc > 2 ? atof(argv[2]) : 500.0, dt;
  int nx, nz, i, j, n, nsteps, s, ci[2], cj[2], ni, nj, cnt = 0, fail = 0;
  struct fd3d_tab *t;
  struct fd3d *f;
  FILE *fr = NULL;
  const char *nm[2] = {"PI", "PA"};
  double sx[2] = {7.5e3, 0.0}, sz[2] = {0.0, 6.0e3};
  int scell[2][4][2], sn[2];
  double emax[3] = {0, 0, 0}, smax[3] = {0, 0, 0};
  double *rec = NULL;
  int nrec = 0;
  int ns;
  int which;
  double e1;
  double e2;
  double e3;
  int fl;
  double ang;
  int ns_2;
  int m1;
  int m2;
  int N;
  int i_2;
  size_t nv;
  double *t1;
  double *t2;
  double es;
  double en;
  double ss;
  double sn_2;
  int isn;
  double d;
  int n_2;
  int k;
  double *t0;
  double d_2;
  int fl_2;
  struct fd3d_tab *t_2;
  double w;
  int fail_2;
  int k_2;
  double cases[6][4];
  double e;
  double x;
  double z;
  double d_3;
  double v;
  double ta;
  int k_3;
  int idx;
  FILE *fo;
  double on;
  double pk;
  double v_2;
  char name[8];
  double r[4];
  int q;
  double c;
  double e_2;
  if(argc >= 3 && !strcmp(argv[1], "tpv101"))
    return run_tpv101(argc, argv);
  if(argc >= 2 && !strcmp(argv[1], "mirror"))
    return run_mirror(argc, argv);
  if(argc >= 3 && !strcmp(argv[1], "kernels"))
    return run_kernels(argv[2]);
  if(argc >= 3 && !strcmp(argv[1], "layer"))
    return run_layer(argv[2]);
  if(argc >= 3 && !strcmp(argv[1], "surf"))
    return run_surf(argv[2]);
  if(argc >= 2 && !strcmp(argv[1], "surfstab"))
    return run_surfstab(argc >= 3 ? atoi(argv[2]) : 4);
  if(argc >= 5 && !strcmp(argv[1], "surfsym"))
    return run_surfsym(atof(argv[2]), atoi(argv[3]), argv[4]);
  if(argc >= 3 && !strcmp(argv[1], "q2d"))
    return run_q2d(argv[2]);
  if(argc >= 2 && !strcmp(argv[1], "pair"))
    return run_pair(argc > 2 ? atoi(argv[2]) : 0);
  if(argc >= 2 && !strcmp(argv[1], "multipair"))
    return run_multi_pair();
  if(argc >= 2 && !strcmp(argv[1], "h6step")){
    ns = argc > 2 ? atoi(argv[2]) : 12;
    which = argc > 3 ? atoi(argv[3]) : -1;
    e1 = (which < 0 || which == 0) ? fd3d_h6_step_check(ns, 0, 0.0, 1) : 0.0;
    e2 = (which < 0 || which == 1) ? fd3d_h6_step_check(ns, 1, 0.0, 1) : 0.0;
    e3 = (which < 0 || which == 2) ? fd3d_h6_step_check(ns, 1, 37.0, 1) : 0.0;
    fl = !(e1 >= 0.0 && e1 < 1e-8) || !(e2 >= 0.0 && e2 < 1e-8) || !(e3 >= 0.0 && e3 < 1e-6);
    printf("fd3d_test h6step: %s\n", fl ? "FAIL" : "PASS");
    return fl;
  }
  if(argc >= 4 && !strcmp(argv[1], "h6strike")){
    /* H6 step 4: tractions on both faults with strips of dx/3 and dx/4 (and dx/m2 if given) for a strike difference */
    /* usage: fd3d_test h6strike <angle> <steps> [m2 (4)] [m1 (3)]: cells dx/m1 against dx/m2 (m = 2 is
       in the surface's unstable band; dx/6 needs several GB for the dense strip tables) */
    ang = atof(argv[2]);
    ns_2 = atoi(argv[3]);
    m1 = argc > 5 ? atoi(argv[5]) : 3;
    m2 = argc > 4 ? atoi(argv[4]) : 4;
    N = 8;
    nv = (size_t)ns_2*2*2*N;
    t1 = malloc(sizeof(double)*nv);
    t2 = malloc(sizeof(double)*nv);
    es = 0;
    en = 0;
    ss = 0;
    sn_2 = 0;
    if(fd3d_h6_strike_run(ang, m1, ns_2, t1, 1) || fd3d_h6_strike_run(ang, m2, ns_2, t2, 1))
      return 1;
    for(i_2 = 0; i_2 < (int)nv; i_2++){
      isn = (i_2 % (2*N)) >= N;
      d = fabs(t1[i_2] - t2[i_2]);
      if(isn){
        if(d > en)
          en = d;
        if(fabs(t2[i_2]) > sn_2)
          sn_2 = fabs(t2[i_2]);
      } else{
        if(d > es)
          es = d;
      }
    }
    {
      t0 = malloc(sizeof(double)*2*N);
      (void)t0;
      for(n_2 = 0; n_2 < ns_2; n_2++)
        for(k = 0; k < 2; k++)
          for(i_2 = 0; i_2 < N; i_2++){
            d_2 = fabs(t2[((size_t)n_2*2 + k)*2*N + i_2] - t2[((size_t)0*2 + k)*2*N + i_2]);
            if(d_2 > ss)
              ss = d_2;
          }
      free(t0);
    }
    h6strike_where(t1, t2, ns_2, N);
    printf("fd3d_test h6strike: turn %g deg: cells dx/%d against dx/%d: shear %.3e of %.3e (largest change of shear), normal %.3e of %.3e\n", ang, m1, m2, es, ss, en, sn_2);
    free(t1);
    free(t2);
    return 0;
  }
  if(argc >= 3 && !strcmp(argv[1], "tpv14"))
    return run_tpv14(atof(argv[2]), argc > 3 ? atoi(argv[3]) : 14, argc > 4 ? atof(argv[4]) : 12.0, argc > 5 ? atof(argv[5]) : 0.0);
  if(argc >= 2 && !strcmp(argv[1], "h6strips")){
    fl_2 = fd3d_h6_strips_check(1);
    printf("fd3d_test h6strips: %s\n", fl_2 ? "FAIL" : "PASS");
    return fl_2;
  }
  if(argc >= 2 && !strcmp(argv[1], "h6frames")){
    /* H6: the generic global-frame response against the one-fault surface code, rotated */
    t_2 = fd3d_tab_create(4, 2, 1000.0, 0.5*1000.0/3464.0, 6000.0, 3464.0, 2670.0, 0);
    w = 0.0;
    fail_2 = 0;
    cases[0][0] = 60;
    cases[0][1] = 2;
    cases[0][2] = 1;
    cases[0][3] = 37;
    cases[1][0] = 30;
    cases[1][1] = 2;
    cases[1][2] = -1;
    cases[1][3] = -71;
    cases[2][0] = 90;
    cases[2][1] = 1;
    cases[2][2] = 1;
    cases[2][3] = 15;
    cases[3][0] = 45;
    cases[3][1] = 1;
    cases[3][2] = -1;
    cases[3][3] = 123;
    cases[4][0] = 60;
    cases[4][1] = 2;
    cases[4][2] = 1;
    cases[4][3] = 0;
    cases[5][0] = 75;
    cases[5][1] = 1;
    cases[5][2] = 1;
    cases[5][3] = 200;
    for(k_2 = 0; k_2 < 6; k_2++){
      e = fd3d_h6_frame_check(t_2, cases[k_2][0], (int)cases[k_2][1], cases[k_2][2], cases[k_2][3], 250.0);
      printf("fd3d_test h6frames: dip %g, comp %g, sign %g, strike %g: largest relative difference %.2e\n", cases[k_2][0], cases[k_2][1], cases[k_2][2], cases[k_2][3], e);
      if(!(e >= 0.0 && e < 1e-10))
        fail_2 = 1;
      if(e > w)
        w = e;
    }
    printf("fd3d_test h6frames: %s\n", fail_2 ? "FAIL" : "PASS");
    fd3d_tab_free(t_2);
    return fail_2;
  }
  if(argc >= 2 && !strcmp(argv[1], "msurf"))
    return run_msurf(argc > 2 ? atoi(argv[2]) : 2, argc > 3 ? atof(argv[3]) : 3.0);
  if(argc >= 2 && !strcmp(argv[1], "multi"))
    return run_multi(argc > 2 ? atof(argv[2]) : 0.0, argc > 3 ? atof(argv[3]) : 90.0, argc > 4 ? atof(argv[4]) : 4.0, argc > 5 ? atoi(argv[5]) : 1, argc > 6 ? atoi(argv[6]) : 0);
  if(argc >= 6 && !strcmp(argv[1], "stepover"))
    return run_stepover(atof(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5]), argc > 6 ? atof(argv[6]) : 0.75);
  if(argc >= 7 && !strcmp(argv[1], "tpv10"))
    return run_tpv10(atof(argv[2]), atof(argv[3]), atoi(argv[4]), atoi(argv[5]), argv[6]);
  if(argc < 3 || strcmp(argv[1], "tpv3")){
    fprintf(stderr, "usage: %s tpv3|tpv101 dx [ref] [tend]\n", argv[0]);
    return 2;
  }
  nx = (int)floor(30e3/dx + 0.5);
  nz = (int)floor(15e3/dx + 0.5);
  dt = 0.5*dx/cT;
  t = fd3d_tab_create(nx, nz, dx, dt, cL, cT, rho, 0);
  fprintf(stderr, "fd3d_test tpv3: %d x %d elements, dt %.4f s, window %d lags, FFT %d x %d, tables %.0f MB, k1 self %.6e, ks self %.6e\n",
      nx, nz, dt, t->L, t->Mx, t->Mz, 16.0*t->L*t->Mx*t->Mz/1e6, t->k1[1][1], t->ks_self);
  f = fd3d_create(t);
  f->law = 0;
  f->mu_s = 0.677;
  f->mu_d = 0.525;
  f->d0 = 0.40;
  f->nl = getenv("FD3D_NL") ? 1 : 0;
  for(i = 0; i < nx; i++)
    for(j = 0; j < nz; j++){
      x = -15e3 + (i + 0.5)*dx;
      z = -7.5e3 + (j + 0.5)*dx;
      f->tau0[i*nz + j] = (fabs(x) < 1.5e3 && fabs(z) < 1.5e3) ? 81.6e6 : 70e6;
      f->sigma[i*nz + j] = 120e6;
    }
  for(s = 0; s < 2; s++){
    ni = near2(sx[s], -15e3, dx, nx, ci);
    nj = near2(sz[s], -7.5e3, dx, nz, cj);
    sn[s] = 0;
    for(i = 0; i < ni; i++)
      for(j = 0; j < nj; j++){
        scell[s][sn[s]][0] = ci[i];
        scell[s][sn[s]][1] = cj[j];
        sn[s]++;
      }
  }
  nsteps = (int)floor(12.0/dt + 0.5);
  rec = malloc(sizeof(double)*2*4*nsteps);
  for(n = 0; n < nsteps; n++){
    fd3d_step(f, NULL);
    for(s = 0; s < 2; s++){
      d_3 = 0;
      v = 0;
      ta = 0;
      for(k_3 = 0; k_3 < sn[s]; k_3++){
        idx = scell[s][k_3][0]*nz + scell[s][k_3][1];
        d_3 += f->delta[idx];
        v += f->V[idx];
        ta += f->tau[idx];
      }
      rec[(s*nsteps + n)*4 + 0] = f->time;
      rec[(s*nsteps + n)*4 + 1] = d_3/sn[s];
      rec[(s*nsteps + n)*4 + 2] = v/sn[s];
      rec[(s*nsteps + n)*4 + 3] = ta/sn[s];
    }
  }
  nrec = nsteps;
  {
    fo = fopen("c_series.txt", "w");
    for(s = 0; s < 2; s++)
      for(n = 0; n < nrec; n++)
        fprintf(fo, "%s %.10e %.12e %.12e %.12e\n", nm[s], rec[(s*nsteps+n)*4], rec[(s*nsteps+n)*4+1], rec[(s*nsteps+n)*4+2], rec[(s*nsteps+n)*4+3]);
    fclose(fo);
  }
  for(s = 0; s < 2; s++){
    on = -1;
    pk = 0;
    for(n = 0; n < nrec; n++){
      v_2 = rec[(s*nsteps + n)*4 + 2];
      if(on < 0 && v_2 > 0.1)
        on = rec[(s*nsteps + n)*4];
      if(v_2 > pk)
        pk = v_2;
    }
    printf("%s: V > 0.1 m/s at %.3f s, peak slip rate %.3f m/s, final slip %.3f m\n", nm[s], on, pk, rec[(s*nsteps + nrec - 1)*4 + 1]);
  }
  if(argc > 3 && (fr = fopen(argv[3], "r"))){
    while(fscanf(fr, "%7s %lf %lf %lf %lf", name, r, r + 1, r + 2, r + 3) == 5){
      s = strcmp(name, "PI") ? 1 : 0;
      n = (int)floor(r[0]/dt + 0.5) - 1;
      if(n < 0 || n >= nrec)
        continue;
      for(q = 0; q < 3; q++){
        c = rec[(s*nsteps + n)*4 + 1 + q];
        e_2 = fabs(c - r[1 + q]);
        if(e_2 > emax[q])
          emax[q] = e_2;
        if(fabs(r[1 + q]) > smax[q])
          smax[q] = fabs(r[1 + q]);
      }
      cnt++;
    }
    fclose(fr);
    printf("fd3d_test tpv3 dx %g vs %s: %d samples; max |diff| slip %.2e (of %.3f), slip rate %.2e (of %.3f), stress %.2e (of %.3e)\n",
        dx, argv[3], cnt, emax[0], smax[0], emax[1], smax[1], emax[2], smax[2]);
    fail = (cnt == 0 || emax[0] > 1e-8*smax[0] || emax[1] > 1e-8*smax[1] || emax[2] > 1e-8*smax[2]);
    printf("fd3d_test tpv3: %s\n", fail ? "FAIL" : "PASS");
  }
  fd3d_free(f);
  fd3d_tab_free(t);
  free(rec);
  return fail;
}
