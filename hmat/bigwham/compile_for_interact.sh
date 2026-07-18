#!/bin/bash
#
# compile_for_interact.sh
#
# Build BigWham (https://github.com/GeoEnergyLab-EPFL/BigWham) as a static
# library for linking into interact, in a similar way that HMMVP/HACApK are built and
# linked. Run this from inside the BigWham checkout that lives at
# interact/bigwham/ :
#
#   cd interact
#   mkdir bigwham
#   cd bigwham
#   git clone https://github.com/GeoEnergyLab-EPFL/BigWham.git
#   ./compile_for_interact.sh
#
# (or `git pull` to update an existing checkout, then re-run). The result is
#   interact/bigwham/BigWham/build/libBigWham.a
# which config/makefile.petsc references through the BIGWHAM_* variables when
# interact is built with USE_BIGWHAM.
#
# NOTE on licensing: BigWham is GPL v3. Linking it into interact makes the
# combined, distributed work subject to GPL v3. Keep the BigWham backend
# optional (compile-time guarded via USE_BIGWHAM) unless that is acceptable.
#
# Dependencies:
#   - CMake >= 3.x and a C++17 compiler
#   - OpenMP
#   - a BLAS that provides the CBLAS header (cblas.h) and LAPACKE (lapacke.h)
#       Ubuntu/Debian:  apt-get install libopenblas-dev liblapacke-dev
#       (these provide /usr/include/.../cblas.h, /usr/include/lapacke.h,
#        and libopenblas with the LAPACKE entry points)
#   - on clusters with MKL (e.g. TACC), set IL_MATH_VENDOR=MKL_sequential and
#       load the MKL module instead of installing OpenBLAS:
#       IL_MATH_VENDOR=MKL_sequential ./compile_for_interact.sh
#
# Environment overrides:
#   IL_MATH_VENDOR   BLAS vendor for the bundled IL library (default OpenBLAS)
#   NP               parallel build jobs (default 4)
#
set -e

IL_MATH_VENDOR=${IL_MATH_VENDOR:-OpenBLAS}
NP=${NP:-4}

echo "compile_for_interact.sh: configuring BigWham (IL_MATH_VENDOR=$IL_MATH_VENDOR)"

#git clone https://github.com/GeoEnergyLab-EPFL/BigWham.git
# sub dir 
cd BigWham
cmake -B build -S . \
  -DBIGWHAM_PYTHON_INTERFACE=OFF \
  -DBIGWHAM_JULIA_INTERFACE=OFF \
  -DBIGWHAM_MATHEMATICA_INTERFACE=OFF \
  -DBIGWHAM_TESTING=OFF \
  -DBIGWHAM_OPENMP=ON \
  -DCMAKE_BUILD_TYPE=Release \
  -DIL_MATH_VENDOR="$IL_MATH_VENDOR"

echo "compile_for_interact.sh: building static library (-j$NP)"
cmake --build build --target BigWhamStatic -j"$NP"

if [ -f build/libBigWham.a ]; then
    echo "compile_for_interact.sh: OK -> $(pwd)/build/libBigWham.a"
    echo "  now build interact with USE_BIGWHAM (see config/makefile.petsc)"
else
    echo "compile_for_interact.sh: ERROR: build/libBigWham.a was not produced" >&2
    exit 1
fi
