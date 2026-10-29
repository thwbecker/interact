/*
  fd_dynamic: standalone single-rupture driver for the fully dynamic
  marcher (slip-weakening friction), for the benchmarks that the
  rate-and-state framework of rsf_solve cannot express.

  usage: fd_dynamic tpv205 [dx] [tend] [ref_file]
         fd_dynamic hd W_m [dx] [tend] [ref_file]
    TPV205-2D (SCEC; 2-D slice of TPV5): 30 km fault, dx (default 50 m),
    mu_s 0.677, mu_d 0.525, d0 0.4 m, sigma 120 MPa, tau0 70 MPa with the
    81.6 MPa nucleation patch (|x| < 1.5 km), 78 MPa (-9 < x < -6 km) and
    62 MPa (6 < x < 9 km) patches; stations at -4.5, -7.5 and +12 km.
    Writes st-045.dat, st-075.dat, st120.dat (t, slip, slip rate, shear
    stress) and tpv205_final_slip.dat. With ref_file (prototype dump:
    lines "station t slip V tau"), prints the maximum differences and
    PASS/FAIL at 1e-8 relative.
  hd: Harris and Day (1993) case D (mu_s 0.75, mu_d 0.51, d0 0.1 m, tau0 20
    MPa, sigma0 33.3 MPa), two 28 km segments with 5 km overlap, step W (m;
    prototype right-lateral convention, +W compressional), forced nucleation
    over 2.5 km at 0.5 c_s, normal-stress coupling (prototype sign
    convention, nsign = -1). Writes hd_trace.dat (t, slip, slip rate and
    normal stress at four points); with ref_file (fully_dynamic/dump_hd.py)
    compares every step at 1e-8 relative. This exercises the cross-segment
    tables and the normal coupling at roundoff level.
*/
static int run_hd(int argc, char **argv);
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
  struct fd_segment seg;
  struct fd_seglist segs;
  struct fd_pair *pair;
  struct fd_sys *s;
  struct fd_state *st;
  struct fd_slipweak fr;
  char err[256];
  FILE *fo[3];
  int *gidx;
  const char *reffile = NULL;
  int e;
  double x;
  double best;
  double x_2;
  char fn[64];
  double vmax;
  int nr;
  FILE *ff;
  FILE *fr_;
  char name[32];
  double t;
  double d;
  double v;
  double tau;
  double md;
  double mv;
  double mt;
  double sd;
  double sv;
  double stau;
  int cnt;
  int fail;
  double **rec;
  int *nrec;
  FILE *fi;
  char fn_2[64];
  if(argc >= 2 && !strcmp(argv[1], "hd"))
    return run_hd(argc, argv);
  if(argc < 2 || strcmp(argv[1], "tpv205")){
    fprintf(stderr, "usage: %s tpv205 [dx] [tend] [ref_file] | hd W_m [dx] [tend] [ref_file]\n", argv[0]);
    return 2;
  }
  if(argc > 2)
    dx = atof(argv[2]);
  if(argc > 3)
    tend = atof(argv[3]);
  if(argc > 4)
    reffile = argv[4];
  n = (int)(L / dx);
  dt = 0.5 * dx / BE;
  nsteps = (int)floor(tend / dt + 0.5);
  seg.n = n;
  seg.dx = dx;
  seg.a = dx / 2;
  seg.s[0] = 1.0;
  seg.s[1] = 0.0;
  seg.X0[0] = -L / 2 + dx / 2;
  seg.X0[1] = 0.0;
  gidx = malloc(sizeof(int) * n);
  for(i = 0; i < n; i++)
    gidx[i] = i;
  segs.nseg = 1;
  segs.seg = &seg;
  segs.gidx = &gidx;
  pair = fd_pair_build(&seg, &seg, FD_MODE_INPLANE, dt, 1.0, AL, BE, RHO, 1e-2, nsteps + 2, 1, 0, err);
  if(!pair){
    fprintf(stderr, "table build failed: %s\n", err);
    return 1;
  }
  s = fd_sys_create(&segs, &pair, 1, NULL, NULL, FD_MODE_INPLANE, dt, AL, BE, RHO, -1);
  if(getenv("FD_CONV_B") && getenv("FD_CONV_B")[0]){
    e = fd_sys_set_partitioned(s, atoi(getenv("FD_CONV_B")), 1);
    fprintf(stderr, "partitioned convolution: %s, B %d, D %d\n", e ? "failed" : (s->conv ? "on" : "not needed"), s->B, s->D);
  }
  st = fd_state_create(n);
  st->normal_coupling = 0;
  st->nl = (getenv("FD_NL") && getenv("FD_NL")[0]) ? 1 : 0;       /* second-order marching */
  fr.N = n;
  fr.dugdale = 0;
  fr.d0 = 0.40;
  fr.force_time = NULL;
  fr.force_t0 = 0.0;
  fr.locked = NULL;
  fr.mu_s = malloc(sizeof(double) * n);
  fr.mu_d = malloc(sizeof(double) * n);
  for(i = 0; i < n; i++){
    x = seg.X0[0] + i * dx;
    fr.mu_s[i] = 0.677;
    fr.mu_d[i] = 0.525;
    st->sigma0[i] = 120e6;
    st->tau0[i] = 70e6;
    if(fabs(x) < 1.5e3)
      st->tau0[i] = 81.6e6;
    if(x > -9e3 && x < -6e3)
      st->tau0[i] = 78e6;
    if(x > 6e3 && x < 9e3)
      st->tau0[i] = 62e6;
  }
  for(k = 0; k < nst; k++){
    best = 1e30;
    ist[k] = 0;
    for(i = 0; i < n; i++){
      x_2 = seg.X0[0] + i * dx;
      if(fabs(x_2 - stx[k]) < best){
        best = fabs(x_2 - stx[k]);
        ist[k] = i;
      }
    }
    fo[k] = fopen(stname[k], "w");
    if(!fo[k]) {
      sprintf(fn, "%s.dat", stname[k]);
      fo[k] = fopen(fn, "w");
    }
  }
  fprintf(stderr, "fd_dynamic tpv205: dx %g n %d dt %.3f ms steps %d window %d (%.1f s) M %d\n", dx, n, dt * 1e3, nsteps, pair->L, pair->L * dt, pair->M);
  for(k = 0; k < nsteps; k++){
    fd_step_slipweak(s, st, &fr, 1);
    for(j = 0; j < nst; j++)
      fprintf(fo[j], "%.10e %.10e %.10e %.10e\n", st->t, st->delta[ist[j]], st->V[ist[j]], st->tau[ist[j]]);
    if(k % 400 == 0){
      vmax = 0;
      nr = 0;
      for(i = 0; i < n; i++){
        if(st->V[i] > vmax)
          vmax = st->V[i];
        if(st->delta[i] > 0.01)
          nr++;
      }
      fprintf(stderr, "  t=%5.2f s  max V=%.3f m/s  ruptured %.1f km\n", st->t, vmax, nr * dx / 1e3);
    }
  }
  for(j = 0; j < nst; j++)
    fclose(fo[j]);
  {
    ff = fopen("tpv205_final_slip.dat", "w");
    for(i = 0; i < n; i++)
      fprintf(ff, "%.6e %.10e\n", seg.X0[0] + i * dx, st->delta[i]);
    fclose(ff);
  }
  if(reffile){
    /* reference lines: station t slip V tau; compare at matching step index */
    fr_ = fopen(reffile, "r");
    md = 0;
    mv = 0;
    mt = 0;
    sd = 0;
    sv = 0;
    stau = 0;
    cnt = 0;
    fail = 0;
    rec = malloc(sizeof(double *) * nst);
    nrec = calloc(nst, sizeof(int));
    for(j = 0; j < nst; j++){
      rec[j] = malloc(sizeof(double) * 4 * nsteps);
      fi = fopen(stname[j], "r");
      if(!fi){
        sprintf(fn_2, "%s.dat", stname[j]);
        fi = fopen(fn_2, "r");
      }
      while(nrec[j] < nsteps && fscanf(fi, "%lf %lf %lf %lf", rec[j] + 4 * nrec[j], rec[j] + 4 * nrec[j] + 1, rec[j] + 4 * nrec[j] + 2, rec[j] + 4 * nrec[j] + 3) == 4)
        nrec[j]++;
      fclose(fi);
    }
    if(!fr_){
      fprintf(stderr, "cannot open %s\n", reffile);
      return 1;
    }
    while(fscanf(fr_, "%31s %lf %lf %lf %lf", name, &t, &d, &v, &tau) == 5){
      for(j = 0; j < nst; j++)
        if(!strcmp(name, stname[j]))
          break;
      if(j == nst)
        continue;
      i = (int)floor(t / dt + 0.5) - 1;
      if(i < 0 || i >= nrec[j])
        continue;
      if(fabs(rec[j][4 * i + 1] - d) > md)
        md = fabs(rec[j][4 * i + 1] - d);
      if(fabs(rec[j][4 * i + 2] - v) > mv)
        mv = fabs(rec[j][4 * i + 2] - v);
      if(fabs(rec[j][4 * i + 3] - tau) > mt)
        mt = fabs(rec[j][4 * i + 3] - tau);
      if(fabs(d) > sd)
        sd = fabs(d);
      if(fabs(v) > sv)
        sv = fabs(v);
      if(fabs(tau) > stau)
        stau = fabs(tau);
      cnt++;
    }
    fclose(fr_);
    printf("fd_dynamic tpv205 dx %g vs %s: %d samples; max |diff| slip %.2e (of %.2f), slip rate %.2e (of %.2f), stress %.2e (of %.2e)\n",
        dx, reffile, cnt, md, sd, mv, sv, mt, stau);
    if(cnt == 0 || md > 1e-8 * sd || mv > 1e-8 * sv || mt > 1e-8 * stau)
      fail = 1;
    printf("fd_dynamic tpv205: %s\n", fail ? "FAIL" : "PASS");
    return fail;
  }
  return 0;
}

static int run_hd(int argc, char **argv)
{
  double AL = 6000.0, BE = 3464.0, RHO = 2700.0, W, dx = 100.0, tend = 8.0, L = 28e3, ov = 5e3, dt;
  int n, N, i, k, nsteps, pts[4], j, fail = 0, cnt = 0;
  struct fd_segment seg[2];
  struct fd_seglist segs;
  struct fd_pair *pairs[4];
  struct fd_sys *s;
  struct fd_state *st;
  struct fd_slipweak fr;
  char err[256];
  int *gidx[2];
  FILE *fo, *fr_;
  const char *reffile = NULL;
  double emax[3] = {0, 0, 0}, smax[3] = {0, 0, 0};
  int e;
  double x;
  char line[4096];
  double r[13];
  double c[13];
  int q;
  double d;
  double m2;
  if(argc < 3){
    fprintf(stderr, "usage: fd_dynamic hd W_m [dx] [tend] [ref_file]\n");
    return 2;
  }
  W = atof(argv[2]);
  if(argc > 3)
    dx = atof(argv[3]);
  if(argc > 4)
    tend = atof(argv[4]);
  if(argc > 5)
    reffile = argv[5];
  n = (int)(L / dx);
  N = 2 * n;
  dt = 0.5 * dx / BE;
  nsteps = (int)floor(tend / dt + 0.5);
  for(j = 0; j < 2; j++){
    seg[j].n = n;
    seg[j].dx = dx;
    seg[j].a = dx / 2;
    seg[j].s[0] = 1.0;
    seg[j].s[1] = 0.0;
    gidx[j] = malloc(sizeof(int) * (size_t)n);
    for(i = 0; i < n; i++)
      gidx[j][i] = j * n + i;
  }
  seg[0].X0[0] = -L / 2 + dx / 2;
  seg[0].X0[1] = 0.0;
  seg[1].X0[0] = L / 2 - ov + dx / 2;
  seg[1].X0[1] = W;
  segs.nseg = 2;
  segs.seg = seg;
  segs.gidx = gidx;
  for(k = 0; k < 2; k++)
    for(j = 0; j < 2; j++){
      pairs[k * 2 + j] = fd_pair_build(seg + k, seg + j, FD_MODE_INPLANE, dt, 1.0, AL, BE, RHO, 1e-2, nsteps + 2, 1, 0, err);
      if(!pairs[k * 2 + j]){
        fprintf(stderr, "table build failed: %s\n", err);
        return 1;
      }
    }
  s = fd_sys_create(&segs, pairs, 2, NULL, NULL, FD_MODE_INPLANE, dt, AL, BE, RHO, -1);
  if(getenv("FD_CONV_B") && getenv("FD_CONV_B")[0]){
    e = fd_sys_set_partitioned(s, atoi(getenv("FD_CONV_B")), 1);
    fprintf(stderr, "partitioned convolution: %s, B %d, D %d\n", e ? "failed" : (s->conv ? "on" : "not needed"), s->B, s->D);
  }
  st = fd_state_create(N);
  st->normal_coupling = getenv("FD_NOCOUP") ? 0 : 1;
  st->nl = (getenv("FD_NL") && getenv("FD_NL")[0]) ? 1 : 0;
  fr.N = N;
  fr.dugdale = 0;
  fr.d0 = 0.10;
  fr.force_t0 = 0.2;
  fr.locked = NULL;
  fr.mu_s = malloc(sizeof(double) * (size_t)N);
  fr.mu_d = malloc(sizeof(double) * (size_t)N);
  fr.force_time = malloc(sizeof(double) * (size_t)N);
  for(i = 0; i < N; i++){
    x = (i < n) ? seg[0].X0[0] + i * dx : seg[1].X0[0] + (i - n) * dx;
    fr.mu_s[i] = 0.75;
    fr.mu_d[i] = 0.51;
    st->sigma0[i] = 33.3e6;
    st->tau0[i] = 20e6;
    fr.force_time[i] = (i < n && fabs(x) < 2.5e3) ? fabs(x) / (0.5 * BE) : INFINITY;
  }
  pts[0] = n - 1;
  pts[1] = n / 2 + 50;
  pts[2] = n + (int)(5.5e3 / dx);
  pts[3] = n + (int)(1e3 / dx);
  fo = fopen("hd_trace.dat", "w");
  fr_ = reffile ? fopen(reffile, "r") : NULL;
  if(reffile && !fr_){
    fprintf(stderr, "cannot open %s\n", reffile);
    return 1;
  }
  if(fr_){
    if(!fgets(line, sizeof(line), fr_))
      return 1;
  }   /* header */
  fprintf(stderr, "fd_dynamic hd: W %g m, dx %g, N %d, steps %d, window %d\n", W, dx, N, nsteps, s->L);
  for(k = 0; k < nsteps; k++){
    fd_step_slipweak(s, st, &fr, 1);
    c[0] = st->t;
    for(j = 0; j < 4; j++){
      c[1 + j] = st->delta[pts[j]];
      c[5 + j] = st->V[pts[j]];
      c[9 + j] = st->sigma[pts[j]];
    }
    for(j = 0; j < 13; j++)
      fprintf(fo, "%.10e%c", c[j], j == 12 ? '\n' : ' ');
    if(fr_){
      for(j = 0; j < 13; j++)
        if(fscanf(fr_, "%lf", r + j) != 1){
          fprintf(stderr, "reference too short\n");
          return 1;
        }
      for(j = 0; j < 4; j++){
        for(q = 0; q < 3; q++){
          d = fabs(c[1 + 4 * q + j] - r[1 + 4 * q + j]);
          if(d > emax[q])
            emax[q] = d;
          if(fabs(r[1 + 4 * q + j]) > smax[q])
            smax[q] = fabs(r[1 + 4 * q + j]);
        }
      }
      cnt++;
    }
  }
  fclose(fo);
  {
    m2 = 0;
    for(i = n; i < N; i++)
      if(st->delta[i] > m2)
        m2 = st->delta[i];
    fprintf(stderr, "fd_dynamic hd: final max slip on segment 2 %.6f m\n", m2);
  }
  if(fr_){
    fclose(fr_);
    printf("fd_dynamic hd W %g dx %g vs %s: %d steps; max |diff| slip %.2e (of %.3f), slip rate %.2e (of %.3f), normal stress %.2e (of %.3e)\n",
        W, dx, reffile, cnt, emax[0], smax[0], emax[1], smax[1], emax[2], smax[2]);
    if(cnt == 0 || emax[0] > 1e-8 * smax[0] || emax[1] > 1e-8 * smax[1] || emax[2] > 1e-8 * smax[2])
      fail = 1;
    printf("fd_dynamic hd: %s\n", fail ? "FAIL" : "PASS");
  }
  return fail;
}
