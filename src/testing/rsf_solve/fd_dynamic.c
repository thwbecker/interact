/*
  fd_dynamic: standalone single-rupture driver for the fully dynamic
  marcher (slip-weakening friction), for the benchmarks that the
  rate-and-state framework of rsf_solve cannot express.

  usage: fd_dynamic tpv205 [dx] [tend] [ref_file]
    TPV205-2D (SCEC; 2-D slice of TPV5): 30 km fault, dx (default 50 m),
    mu_s 0.677, mu_d 0.525, d0 0.4 m, sigma 120 MPa, tau0 70 MPa with the
    81.6 MPa nucleation patch (|x| < 1.5 km), 78 MPa (-9 < x < -6 km) and
    62 MPa (6 < x < 9 km) patches; stations at -4.5, -7.5 and +12 km.
    Writes st-045.dat, st-075.dat, st120.dat (t, slip, slip rate, shear
    stress) and tpv205_final_slip.dat. With ref_file (prototype dump:
    lines "station t slip V tau"), prints the maximum differences and
    PASS/FAIL at 1e-8 relative.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fd_march.h"

int main(int argc, char **argv)
{
  double AL = 6000.0, BE = 3464.0, RHO = 2670.0, dx = 50.0, tend = 12.0, L = 30e3, dt;
  int n, i, k, nsteps, nst = 3, ist[3], j;
  const char *stname[3] = {"st-045", "st-075", "st120"};
  double stx[3] = {-4.5e3, -7.5e3, 12.0e3};
  struct fd_segment seg; struct fd_seglist segs; struct fd_pair *pair; struct fd_sys *s; struct fd_state *st;
  struct fd_slipweak fr; char err[256]; FILE *fo[3]; int *gidx;
  const char *reffile = NULL;
  if(argc < 2 || strcmp(argv[1], "tpv205")){ fprintf(stderr, "usage: %s tpv205 [dx] [tend] [ref_file]\n", argv[0]); return 2; }
  if(argc > 2) dx = atof(argv[2]);
  if(argc > 3) tend = atof(argv[3]);
  if(argc > 4) reffile = argv[4];
  n = (int)(L / dx); dt = 0.5 * dx / BE; nsteps = (int)floor(tend / dt + 0.5);
  seg.n = n; seg.dx = dx; seg.a = dx / 2; seg.s[0] = 1.0; seg.s[1] = 0.0; seg.X0[0] = -L / 2 + dx / 2; seg.X0[1] = 0.0;
  gidx = malloc(sizeof(int) * n); for(i = 0; i < n; i++) gidx[i] = i;
  segs.nseg = 1; segs.seg = &seg; segs.gidx = &gidx;
  pair = fd_pair_build(&seg, &seg, FD_MODE_INPLANE, dt, 1.0, AL, BE, RHO, 1e-2, nsteps + 2, 1, 0, err);
  if(!pair){ fprintf(stderr, "table build failed: %s\n", err); return 1; }
  s = fd_sys_create(&segs, &pair, FD_MODE_INPLANE, dt, AL, BE, RHO, -1);
  st = fd_state_create(n);
  st->normal_coupling = 0;
  fr.N = n; fr.dugdale = 0; fr.d0 = 0.40; fr.force_time = NULL; fr.force_t0 = 0.0; fr.locked = NULL;
  fr.mu_s = malloc(sizeof(double) * n); fr.mu_d = malloc(sizeof(double) * n);
  for(i = 0; i < n; i++){
    double x = seg.X0[0] + i * dx;
    fr.mu_s[i] = 0.677; fr.mu_d[i] = 0.525;
    st->sigma0[i] = 120e6;
    st->tau0[i] = 70e6;
    if(fabs(x) < 1.5e3) st->tau0[i] = 81.6e6;
    if(x > -9e3 && x < -6e3) st->tau0[i] = 78e6;
    if(x > 6e3 && x < 9e3) st->tau0[i] = 62e6;
  }
  for(k = 0; k < nst; k++){
    double best = 1e30; ist[k] = 0;
    for(i = 0; i < n; i++){ double x = seg.X0[0] + i * dx; if(fabs(x - stx[k]) < best){ best = fabs(x - stx[k]); ist[k] = i; } }
    fo[k] = fopen(stname[k], "w"); if(!fo[k]) { char fn[64]; sprintf(fn, "%s.dat", stname[k]); fo[k] = fopen(fn, "w"); }
  }
  fprintf(stderr, "fd_dynamic tpv205: dx %g n %d dt %.3f ms steps %d window %d (%.1f s) M %d\n", dx, n, dt * 1e3, nsteps, pair->L, pair->L * dt, pair->M);
  for(k = 0; k < nsteps; k++){
    fd_step_slipweak(s, st, &fr, 1);
    for(j = 0; j < nst; j++) fprintf(fo[j], "%.10e %.10e %.10e %.10e\n", st->t, st->delta[ist[j]], st->V[ist[j]], st->tau[ist[j]]);
    if(k % 400 == 0){
      double vmax = 0; int nr = 0;
      for(i = 0; i < n; i++){ if(st->V[i] > vmax) vmax = st->V[i]; if(st->delta[i] > 0.01) nr++; }
      fprintf(stderr, "  t=%5.2f s  max V=%.3f m/s  ruptured %.1f km\n", st->t, vmax, nr * dx / 1e3);
    }
  }
  for(j = 0; j < nst; j++) fclose(fo[j]);
  {
    FILE *ff = fopen("tpv205_final_slip.dat", "w");
    for(i = 0; i < n; i++) fprintf(ff, "%.6e %.10e\n", seg.X0[0] + i * dx, st->delta[i]);
    fclose(ff);
  }
  if(reffile){
    /* reference lines: station t slip V tau; compare at matching step index */
    FILE *fr_ = fopen(reffile, "r"); char name[32]; double t, d, v, tau, md = 0, mv = 0, mt = 0, sd = 0, sv = 0, stau = 0;
    int cnt = 0, fail = 0;
    double **rec = malloc(sizeof(double *) * nst); int *nrec = calloc(nst, sizeof(int));
    for(j = 0; j < nst; j++){ rec[j] = malloc(sizeof(double) * 4 * nsteps); FILE *fi = fopen(stname[j], "r"); if(!fi){ char fn[64]; sprintf(fn, "%s.dat", stname[j]); fi = fopen(fn, "r"); }
      while(nrec[j] < nsteps && fscanf(fi, "%lf %lf %lf %lf", rec[j] + 4 * nrec[j], rec[j] + 4 * nrec[j] + 1, rec[j] + 4 * nrec[j] + 2, rec[j] + 4 * nrec[j] + 3) == 4) nrec[j]++;
      fclose(fi); }
    if(!fr_){ fprintf(stderr, "cannot open %s\n", reffile); return 1; }
    while(fscanf(fr_, "%31s %lf %lf %lf %lf", name, &t, &d, &v, &tau) == 5){
      for(j = 0; j < nst; j++) if(!strcmp(name, stname[j])) break;
      if(j == nst) continue;
      i = (int)floor(t / dt + 0.5) - 1;
      if(i < 0 || i >= nrec[j]) continue;
      if(fabs(rec[j][4 * i + 1] - d) > md) md = fabs(rec[j][4 * i + 1] - d);
      if(fabs(rec[j][4 * i + 2] - v) > mv) mv = fabs(rec[j][4 * i + 2] - v);
      if(fabs(rec[j][4 * i + 3] - tau) > mt) mt = fabs(rec[j][4 * i + 3] - tau);
      if(fabs(d) > sd) sd = fabs(d);
      if(fabs(v) > sv) sv = fabs(v);
      if(fabs(tau) > stau) stau = fabs(tau);
      cnt++;
    }
    fclose(fr_);
    printf("fd_dynamic tpv205 dx %g vs %s: %d samples; max |diff| slip %.2e (of %.2f), slip rate %.2e (of %.2f), stress %.2e (of %.2e)\n",
           dx, reffile, cnt, md, sd, mv, sv, mt, stau);
    if(cnt == 0 || md > 1e-8 * sd || mv > 1e-8 * sv || mt > 1e-8 * stau) fail = 1;
    printf("fd_dynamic tpv205: %s\n", fail ? "FAIL" : "PASS");
    return fail;
  }
  return 0;
}
