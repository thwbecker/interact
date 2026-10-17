/*
  test of fd_kernels/fd_tables against table dumps of the Python prototype
  (fully_dynamic/dump_tables.py). usage: test_fd_tables ref_file [tol]
  compares the window length, static rows and every W entry; exit 0 on PASS
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "fd_tables.h"

int main(int argc, char **argv)
{
  FILE *f; int mode, ns, nr, lmax, M, Lref, i, l, nm, mm, fail = 0, idx;
  double dx, dt, et, al, be, rho, tol, ks, kn, v, wmax = 0, emax = 0, kmax = 0, kerr = 0, rtol = argc > 2 ? atof(argv[2]) : 1e-10;
  struct fd_segment s, r; struct fd_pair *p; char err[256];
  double _Complex *c;
  if(argc < 2){ fprintf(stderr, "usage: %s ref_file [tol]\n", argv[0]); return 2; }
  f = fopen(argv[1], "r"); if(!f){ perror(argv[1]); return 2; }
  if(fscanf(f, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf", &mode, &ns, &nr, &dx, &dt, &et, &al, &be, &rho,
            &tol, &lmax, &M, &Lref, &s.X0[0], &s.X0[1], &s.s[0], &s.s[1], &r.X0[0], &r.X0[1], &r.s[0], &r.s[1]) != 21){ fprintf(stderr, "bad header\n"); return 2; }
  /* the prototype's X0 is the segment start; element 0 centre is half an element along */
  s.n = ns; s.dx = dx; s.a = dx / 2; s.X0[0] += 0.5 * dx * s.s[0]; s.X0[1] += 0.5 * dx * s.s[1];
  r.n = nr; r.dx = dx; r.a = dx / 2; r.X0[0] += 0.5 * dx * r.s[0]; r.X0[1] += 0.5 * dx * r.s[1];
  p = fd_pair_build(&s, &r, mode, dt, et, al, be, rho, tol, lmax, 1, M, err);
  if(!p){ fprintf(stderr, "build failed: %s\n", err); return 1; }
  nm = ns + nr - 1;
  printf("%s: mode %d ns %d nr %d: L %d (ref %d), M %d, nm %d\n", argv[1], mode, ns, nr, p->L, Lref, p->M, nm);
  if(p->L != Lref) fail = 1;
  for(i = 0; i < nm; i++){
    if(fscanf(f, "%d %lf %lf", &mm, &ks, &kn) != 3){ fprintf(stderr, "short static rows\n"); return 2; }
    if(mm != p->m[i]) fail = 1;
    kmax = fmax(kmax, fmax(fabs(ks), fabs(kn)));
    kerr = fmax(kerr, fmax(fabs(ks - p->ks_row[i]), fabs(kn - p->kn_row[i])));
  }
  printf("  static rows: max abs err %.2e relative to max %.2e -> %.2e\n", kerr, kmax, kerr / kmax);
  if(kerr > rtol * kmax) fail = 1;
  c = malloc(sizeof(*c) * p->M);
  for(l = 0; l < Lref && l < p->L; l++){
    fd_fft_exec(p->pb, p->Ws_hat + (size_t)l * p->M, c);
    for(i = 0; i < nm; i++){
      if(fscanf(f, "%lf", &v) != 1){ fprintf(stderr, "short W shear\n"); return 2; }
      idx = p->m[i] >= 0 ? p->m[i] : p->M + p->m[i];
      wmax = fmax(wmax, fabs(v)); emax = fmax(emax, fabs(v - creal(c[idx]) / p->M));
    }
  }
  printf("  W shear: max abs err %.2e relative to max %.2e -> %.2e\n", emax, wmax, emax / wmax);
  if(emax > rtol * wmax) fail = 1;
  if(p->has_normal){
    wmax = emax = 0;
    for(l = 0; l < Lref && l < p->L; l++){
      fd_fft_exec(p->pb, p->Wn_hat + (size_t)l * p->M, c);
      for(i = 0; i < nm; i++){
        if(fscanf(f, "%lf", &v) != 1){ fprintf(stderr, "short W normal\n"); return 2; }
        idx = p->m[i] >= 0 ? p->m[i] : p->M + p->m[i];
        wmax = fmax(wmax, fabs(v)); emax = fmax(emax, fabs(v - creal(c[idx]) / p->M));
      }
    }
    printf("  W normal: max abs err %.2e relative to max %.2e -> %.2e\n", emax, wmax, emax / wmax);
    if(emax > rtol * wmax) fail = 1;
  }
  free(c); fd_pair_free(p); fclose(f);
  printf("test_fd_tables %s: %s\n", argv[1], fail ? "FAIL" : "PASS");
  return fail;
}
