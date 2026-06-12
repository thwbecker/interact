/*
  hmmvp_c_shim.cpp: minimal C interface to hmmvp (A.M. Bradley,
  https://github.com/ambrad/hmmvp, EPL-1.0) for use from
  compress_interaction_matrix.c (-use_hmatrix 4)

  - the GreensFn subclass calls interact's index based ckernel_func
    (the same kernel that fills the dense matrix and drives the
    HACApK construction), converting hmmvp's 1-based block index sets
    to 0-based entries
  - compression uses the MREM whole-matrix Frobenius tolerance
    (tm_mrem_fro): ||B - A||_F <= tol ||B||_F, i.e. tol bounds exactly
    the global relative error that the test tool measures
  - build with -std=c++14 (hmmvp uses pre-C++17 dynamic exception
    specifications) and the hmmvp include paths, link
    -lhmmvp_omp -lstdc++ -fopenmp (see makefile.petsc)

  compile example:
    mpicxx -std=c++14 -O2 -I$(HMMVP_DIR)/hmmvp/include \
       -I$(HMMVP_DIR)/util/include -c hmmvp_c_shim.cpp
*/
#include <vector>
#include <cstdio>
#include <mutex>
#include "Compress.hpp"
#include "Hmat.hpp"

extern "C" {
  /* interact's per-entry kernel (0-based), defined in
     compress_interaction_matrix.c */
  double ckernel_func(int, int, void *);

  void *chmmvp_compress_in_memory(int, double *, double *, double *,
				  double, double, int, void *);
  void chmmvp_mvp(void *, double *, double *);
  void chmmvp_get_info(void *, int *, int *, long *);
  void chmmvp_delete(void *);
}

namespace {
class InteractGF : public hmmvp::GreensFn {
public:
  InteractGF(void *ctx) : kctx(ctx) {}
  virtual bool Call(const hmmvp::CompressBlockInfo &cbi,
		    const std::vector<hmmvp::UInt> &rs,
		    const std::vector<hmmvp::UInt> &cs, double *B) const {
    /* 
       B is column-major with the fast index on rs; hmmvp indices are
       1-based. kernel calls MUST be serialized: interact's Green's
       function chain is not thread safe (the Okada dc3d routines use
       COMMON blocks) and hmmvp's threaded compression calls this
       concurrently. the mutex costs construction the kernel-call
       parallelism (ACA algebra still threads); the MVP never calls
       the kernel and threads freely. NOTE: dc3d.F now carries
       THREADPRIVATE COMMON blocks (compile with -fopenmp), making
       the rectangular Okada chain thread safe - the mutex is
       therefore DISABLED; re-enable it if unaudited element types
       (e.g. triangles) are used with threading:
       std::lock_guard<std::mutex> lock(kmutex);
    */
    for (std::size_t ic = 0; ic < cs.size(); ic++)
      for (std::size_t ir = 0; ir < rs.size(); ir++)
	*B++ = ckernel_func((int)rs[ir] - 1, (int)cs[ic] - 1, kctx);
    return true;
  }
private:
  void *kctx;
  mutable std::mutex kmutex;
};
}

/*
  build the H matrix in memory for n points with coordinates x,y,z,
  whole-matrix relative Frobenius tolerance tol, admissibility eta
  (hmmvp default 3), nthreads OpenMP threads, and the interact kernel
  context kctx; returns an opaque hmmvp::Hmat handle (NULL on error)
*/
void *chmmvp_compress_in_memory(int n, double *x, double *y, double *z,
				double tol, double eta, int nthreads,
				void *kctx)
{
  try {
    hmmvp::Matrix<double> D(3, n);
    for (int i = 0; i < n; i++) {
      D(1, i + 1) = x[i];
      D(2, i + 1) = y[i];
      D(3, i + 1) = z[i];
    }
    hmmvp::Hd *hd = hmmvp::NewHd(D, NULL, eta);
    InteractGF gf(kctx);
    hmmvp::Compressor *c = hmmvp::NewCompressor(hd, &gf);
    c->SetTolMethod(hmmvp::Compressor::tm_mrem_fro);
    c->SetBfroEstimate(c->EstimateBfro()); /* operator is diagonally
					      dominant, so this
					      estimate is reliable */
    c->SetTol(tol);
    c->SetOmpNthreads(nthreads);
    c->AvoidRedundantGfCalls(true);
    hmmvp::Hmat *hm = c->CompressInMemory(1, nthreads);
    hmmvp::DeleteCompressor(c);
    hmmvp::DeleteHd(hd);
    return (void *)hm;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "chmmvp_compress_in_memory: exception: %s\n",
		 e.what());
    return NULL;
  } catch (...) {
    std::fprintf(stderr, "chmmvp_compress_in_memory: unknown exception\n");
    return NULL;
  }
}

/* y = A x with full (global) vectors */
void chmmvp_mvp(void *hm, double *x, double *y)
{
  ((hmmvp::Hmat *)hm)->Mvp(x, y, 1);
}

/* matrix dimensions and number of stored scalars (for the
   compression ratio) */
void chmmvp_get_info(void *hm, int *m, int *n, long *nnz)
{
  hmmvp::Hmat *h = (hmmvp::Hmat *)hm;
  *m = (int)h->GetM();
  *n = (int)h->GetN();
  *nnz = (long)h->GetNnz();
}

void chmmvp_delete(void *hm)
{
  hmmvp::DeleteHmat((hmmvp::Hmat *)hm);
}
