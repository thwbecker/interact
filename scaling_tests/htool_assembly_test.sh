#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
#
#makefault -n 200 -m 80 -l 50 -w 20 -z -20 -strike 0 -dip 90 > /tmp/g.in   # ~16000 patches, 0.5 km cells
makefault -n 400 -m 160 -l 50 -w 20 -z -20 -strike 0 -dip 90 > /tmp/g.in   # 

for P in 1 2 4 8 16 24 48; do
  asm=$(mpirun --bind-to core --map-by core -np $P compress_interaction_matrix \
        -geom_file /tmp/g.in -use_hmatrix 1 -mat_htool_epsilon 3e-5 -mat_htool_eta 100 \
        -nrandom 0 -skip_dense  2>&1 | awk '/H matrix assembly took/{print $(NF-1)}')
  echo "np=$P  htool_assembly_s=$asm"
done

