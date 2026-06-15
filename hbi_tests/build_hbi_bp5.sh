#!/usr/bin/env bash
#
# build_hbi_bp5.sh -- reproducibly build HBI (sozawa94/hbi) for the BP5 comparison,
#   starting from a pristine clone. Clones HBI (if absent), applies the one-line
#   HACApK bounds fix needed to run BP5, and compiles -> ./hbi/lhbiem.
#
# The fix is embedded below (no external .patch needed). It is the ONLY change to
# HBI's sources; everything else is stock. The fix guards a rank that owns zero
# near-field rows in HACApK's param(61)==3 diagonal-scaling path:
#     m_HACApK_use.f90:  do il=ndnr_s,ndnr_e  ->  do il=max(ndnr_s,1),ndnr_e
#
# PARAMETERS ARE POSITIONAL (no environment-variable overrides, so nothing in your
# shell can silently change the build). Edit the CONFIG block for paths/flags.
#   $1 = eps_h  (optional) sets HBI's *default* ACA tolerance at main_LH.f90:145.
#
# NOTE (recommended): you do NOT need to set eps_h at build time. read_inputfile()
# runs (main_LH.f90:171) before matrix generation (493), so an `eps_h <val>` line in
# the .in file overrides the default with no recompile -- and hbi_bp5_scaling_test.sh
# injects exactly that from its $3 argument. So just build once (no $1) and set the
# tolerance per run. For reference at BP5 1 km / N=4000: stock 1d-4 is near-dense
# (~25.3 MB), eps_h=1d-1 -> ~12.3 MB (matches rsf_solve HACApK ztol 1e-1 ~12.1 MB).
# The $1 argument is retained only as a convenience for changing the compiled default.
#
# The ONLY environment variable consulted is $MPIF90 (your MPI Fortran wrapper, e.g.
# PETSc's build/bin/mpif90); it falls back to mpifort if unset. Set it explicitly:
#   MPIF90=/path/to/mpif90 ./build_hbi_bp5.sh          # build once; set eps_h per run
# ---------------------------------------------------------------------------
set -euo pipefail

# ============================ CONFIG =======================================
HBI_REPO="https://github.com/sozawa94/hbi.git"
HBI_DIR="./hbi"
# $MPIF90 is the one intentional environment input (your toolchain wrapper);
# fall back to mpifort only if it is unset. We do NOT read $F90/$LDFLAGS - those
# are commonly exported for other builds and would silently hijack this one.
: "${MPIF90:=mpifort}"
FFLAGS="-fallow-argument-mismatch -ffree-form -ffree-line-length-none -O3 -march=native -fopenmp"
LDLIBS="-llapack -lblas"          # HACApK needs BLAS/LAPACK; edit to -lopenblas / -qmkl as needed
EPS_H="${1-}"                      # $1: empty -> stock 1d-4 ; e.g. 1d-1 -> match rsf_solve band
# ===========================================================================

# --- 1. pristine clone (skip if a tree is already there) ---
if [ -f "$HBI_DIR/main_LH.f90" ]; then
  echo "[1/5] using existing HBI tree at $HBI_DIR"
else
  echo "[1/5] cloning $HBI_REPO -> $HBI_DIR"
  git clone --depth 1 "$HBI_REPO" "$HBI_DIR"
fi
cd "$HBI_DIR"

# --- 2. write the embedded fix ---
echo "[2/5] writing hbi_bp5_fix.patch"
cat > hbi_bp5_fix.patch <<'PATCH_EOF'
--- m_HACApK_use.f90.orig
+++ m_HACApK_use.f90
@@ -191,7 +191,7 @@
  elseif(st_ctl%param(61)==3)then
    ndnr_s=st_ctl%lpmd(6); ndnr_e=st_ctl%lpmd(7); ndnr=st_ctl%lpmd(5)
    allocate(st_bemv%ao(nd)); st_bemv%ao(:)=0.0d0; zsqnd=sqrt(real(nd))
-   do il=ndnr_s,ndnr_e
+   do il=max(ndnr_s,1),ndnr_e
      zad=HACApK_entry_ij(il,il,st_bemv)
      st_bemv%ao(il)=1.0d0/dsqrt(zad/zsqnd)
    enddo
PATCH_EOF

# --- 3. apply it idempotently (-p0 from the repo root) ---
echo "[3/5] applying fix"
if patch -p0 -N --dry-run < hbi_bp5_fix.patch >/dev/null 2>&1; then
  patch -p0 -N < hbi_bp5_fix.patch
elif patch -p0 -R --dry-run < hbi_bp5_fix.patch >/dev/null 2>&1; then
  echo "      already applied, skipping"
else
  echo "ERROR: fix does not apply cleanly to this HBI version - inspect m_HACApK_use.f90 near the param(61)==3 block" >&2
  exit 1
fi
grep -q "max(ndnr_s,1)" m_HACApK_use.f90 || { echo "ERROR: fix not present after apply" >&2; exit 1; }

# --- 4. optional eps_h override (matched-accuracy comparison) ---
if [ -n "$EPS_H" ]; then
  echo "[4/5] setting eps_h=$EPS_H in main_LH.f90 (matched-accuracy build)"
  sed -i -E "s/(; *eps_h=)[0-9.dDeE+-]+/\1${EPS_H}/" main_LH.f90
  grep -nE "eps_r=.*eps_h=" main_LH.f90 | head -1
else
  echo "[4/5] keeping stock eps_h (1d-4, near-dense for the smooth BP5 kernel)"
fi

# --- 5. build (serial: OBJS order encodes Fortran module deps; do NOT use -j) ---
# make's command-line F90= / LDFLAGS= override BOTH the Makefile and the environment,
# so an exported F90=gfortran etc. cannot leak in here.
echo "[5/5] building with F90=\"$MPIF90 $FFLAGS\"  LDFLAGS=\"$LDLIBS\""
make clean >/dev/null 2>&1 || true
make F90="$MPIF90 $FFLAGS" LDFLAGS="$LDLIBS"

if [ -x ./lhbiem ]; then
  echo ""
  echo "BUILD OK -> $(pwd)/lhbiem"
  echo "point the scaling script at it:   HBI=$(pwd)/lhbiem  (edit its CONFIG)"
  [ -n "$EPS_H" ] && echo "built at eps_h=$EPS_H  -> set HBI_EPS_H=$EPS_H label in the scaling script"
else
  echo "ERROR: lhbiem not produced" >&2; exit 1
fi
