#!/usr/bin/env bash
#
# quad2trimesh.sh - mesh a planar quad (4 corners in 3-D) with triangles using
# Shewchuk's Triangle, then print the triangle corners in 3-D.
#
# Triangle is purely 2-D, so the script:
#   1. builds an orthonormal in-plane basis (u_hat, v_hat) from the 4 points,
#   2. projects the corners to 2-D (u, v),
#   3. meshes in 2-D (see modes below),
#   4. lifts every output vertex back with  P = O + u*u_hat + v*v_hat.
#
# MESH MODE - selected by the SIGN of the quality-angle argument:
#   qangle < 0   EQUILATERAL LATTICE. Seed the interior with a triangular
#                lattice at the chosen spacing plus matching boundary points,
#                then take the constrained Delaunay with no refinement. The
#                interior is equilateral up to floating point; a thin rim of
#                boundary triangles carries the 60-vs-90 mismatch. Magnitude
#                of qangle is ignored here.
#   qangle >= 0  QUALITY (irregular). Hand Triangle just the 4 corners and let
#                it refine with -q<qangle> -a<area>. Unstructured mesh.
#
# SNAP - selected by the SIGN of dx (lattice mode only):
#   dx > 0       use |dx| as the lattice spacing as given.
#   dx < 0       snap: treat |dx| as a target, pick an integer number of
#                columns and rows near it, and use the resulting spacing so the
#                rectangle is tiled EXACTLY (no ragged partial cells; left/right
#                edges become clean half-triangles). A rectangle cannot be tiled
#                by perfectly equilateral triangles, so the snapped triangles
#                are slightly isoceles; the chosen counts and the deviation are
#                printed on stderr. (Snap assumes the quad is a rectangle.)
#
# Usage:
#   ./quad2trimesh.sh dx [qangle]  < vertices.txt
#   printf '%s\n' "0 0 0" "4 0 0" "4 3 0" "0 3 0" | ./quad2trimesh.sh 1.0       # lattice, spacing as given
#   printf '%s\n' "0 0 0" "4 0 0" "4 3 0" "0 3 0" | ./quad2trimesh.sh -1.0      # lattice, snapped to exact fit
#   printf '%s\n' "0 0 0" "4 0 0" "4 3 0" "0 3 0" | ./quad2trimesh.sh 1.0 30    # old quality mesh
#
#     dx      target equilateral edge length; negative => snap (lattice mode)
#     qangle  quality angle; default -1 (lattice). Positive values above ~33
#             may prevent Triangle from terminating in quality mode.
#     stdin   four vertices, one "x y z" per line, CCW:
#               bottom-left, bottom-right, top-right, top-left
#             (blank lines and lines starting with # are ignored)
#
# Output (stdout): one triangle per line, nine numbers -
#   ax ay az  bx by bz  cx cy cz
#
set -euo pipefail

command -v triangle >/dev/null 2>&1 || {
  echo "error: 'triangle' (Shewchuk's Triangle) not found in PATH" >&2
  exit 1
}

DX=${1-0.1}
QANGLE=${2--1}

# Lattice mode only: how far (as a fraction of the column spacing) an interior
# lattice point must sit from every edge to be kept. Stay below 0.5 so the
# offset rows' edge-most points (which lie half a column from a vertical edge)
# survive.
MARGIN_FRAC=0.30

# snap from the sign of dx; mesh mode from the sign of qangle (string tests
# avoid float comparison)
case "$DX"     in -*) SNAP=1; DXABS=${DX#-} ;; *) SNAP=0; DXABS=$DX ;; esac
case "$QANGLE" in -*) MODE=lattice          ;; *) MODE=quality       ;; esac

# --- read four vertices from stdin: BL, BR, TR, TL ------------------------
verts=(); n=0
while IFS= read -r line || [ -n "$line" ]; do
  case "$line" in ''|'#'*) continue ;; esac
  read -r a b c _ <<<"$line"
  if [ -z "$a" ] || [ -z "$b" ] || [ -z "$c" ]; then
    echo "error: each vertex line needs three numbers 'x y z'" >&2; exit 1
  fi
  verts+=("$a" "$b" "$c"); n=$((n+1))
  [ "$n" -eq 4 ] && break
done
if [ "$n" -ne 4 ]; then
  echo "error: expected 4 vertices on stdin, got $n" >&2; exit 1
fi
X0=${verts[0]}; Y0=${verts[1]}; Z0=${verts[2]}
X1=${verts[3]}; Y1=${verts[4]}; Z1=${verts[5]}
X2=${verts[6]}; Y2=${verts[7]}; Z2=${verts[8]}
X3=${verts[9]}; Y3=${verts[10]}; Z3=${verts[11]}

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
POLY="$WORK/mesh.poly"

# --- steps 1-3: basis, projection, write .poly for the chosen mode --------
# awk writes the .poly file and prints: "area ox oy oz ux uy uz vx vy vz"
read -r AREA OX OY OZ UX UY UZ VX VY VZ < <(
  awk -v x0="$X0" -v y0="$Y0" -v z0="$Z0" \
      -v x1="$X1" -v y1="$Y1" -v z1="$Z1" \
      -v x2="$X2" -v y2="$Y2" -v z2="$Z2" \
      -v x3="$X3" -v y3="$Y3" -v z3="$Z3" \
      -v dx="$DXABS" -v snap="$SNAP" -v mfrac="$MARGIN_FRAC" \
      -v mode="$MODE" -v poly="$POLY" '
    function dot(ax,ay,az,bx,by,bz){ return ax*bx+ay*by+az*bz }
    # inside the CCW quad by at least margin m (cu[],cv[] are the 4 corners)
    function inside(pu,pv,m,   k,k2,ex,ey,len,cr){
      for(k=0;k<4;k++){
        k2=(k+1)%4
        ex=cu[k2]-cu[k]; ey=cv[k2]-cv[k]
        len=sqrt(ex*ex+ey*ey)
        if(len==0) continue
        cr=ex*(pv-cv[k])-ey*(pu-cu[k])      # >0 means left of edge (inside)
        if(cr/len <= m) return 0
      }
      return 1
    }
    BEGIN{
      # two adjacent edges leaving P0 (bottom edge and left edge)
      ax=x1-x0; ay=y1-y0; az=z1-z0
      bx=x3-x0; by=y3-y0; bz=z3-z0
      nx=ay*bz-az*by; ny=az*bx-ax*bz; nz=ax*by-ay*bx   # normal a x b
      la=sqrt(ax*ax+ay*ay+az*az); ln=sqrt(nx*nx+ny*ny+nz*nz)
      if(la==0||ln==0){ print "error: degenerate / collinear points" > "/dev/stderr"; exit 1 }
      ux=ax/la; uy=ay/la; uz=az/la                     # u_hat along bottom edge
      nhx=nx/ln; nhy=ny/ln; nhz=nz/ln                  # n_hat unit normal
      vx=nhy*uz-nhz*uy; vy=nhz*ux-nhx*uz; vz=nhx*uy-nhy*ux   # v_hat = n_hat x u_hat

      cu[0]=0; cv[0]=0
      cu[1]=dot(x1-x0,y1-y0,z1-z0,ux,uy,uz); cv[1]=dot(x1-x0,y1-y0,z1-z0,vx,vy,vz)
      cu[2]=dot(x2-x0,y2-y0,z2-z0,ux,uy,uz); cv[2]=dot(x2-x0,y2-y0,z2-z0,vx,vy,vz)
      cu[3]=dot(x3-x0,y3-y0,z3-z0,ux,uy,uz); cv[3]=dot(x3-x0,y3-y0,z3-z0,vx,vy,vz)

      area = sqrt(3)/4 * dx*dx     # used by -a in quality mode

      if(mode=="quality"){
        # original PSLG: just the 4 corners and 4 boundary segments
        s="4 2 0 0\n"
        for(i=0;i<4;i++) s=s sprintf("%d %.12g %.12g\n", i+1, cu[i], cv[i])
        s=s "4 0\n1 1 2\n2 2 3\n3 3 4\n4 4 1\n0\n"
        printf "%s", s > poly
      } else {
        # ---- lattice mode ----
        r=sqrt(3)/2                 # equilateral row-height / spacing ratio
        if(snap){
          # rectangle extents from the bottom and left edges, in (u,v)
          W=sqrt(cu[1]*cu[1]+cv[1]*cv[1])
          H=sqrt(cu[3]*cu[3]+cv[3]*cv[3])
          ncol=int(W/dx+0.5);     if(ncol<1)ncol=1
          nrow=int(H/(dx*r)+0.5); if(nrow<1)nrow=1
          su=W/ncol; sv=H/nrow      # column / row spacing for an exact fit
          slant=sqrt(su*su/4+sv*sv)
          dev=(slant/su-1)*100; if(dev<0)dev=-dev
          printf "snap: %d cols x %d rows; su=%.6g sv=%.6g; edge ~%.6g (%.2f%% from equilateral)\n", \
                 ncol,nrow,su,sv,slant,dev > "/dev/stderr"
        } else {
          su=dx; sv=dx*r            # exact equilateral lattice at spacing dx
        }
        margin = mfrac*su

        # boundary vertices around the perimeter. Roughly-horizontal edges
        # split at su (columns); roughly-vertical edges at sv (rows), so the
        # rim lines up with the interior lattice.
        B=0
        for(k=0;k<4;k++){
          k2=(k+1)%4
          ex=cu[k2]-cu[k]; ey=cv[k2]-cv[k]
          len=sqrt(ex*ex+ey*ey)
          aex=(ex<0?-ex:ex); aey=(ey<0?-ey:ey)
          target=(aey>aex? sv : su)
          N=int(len/target+0.5); if(N<1)N=1
          B++; bu[B]=cu[k]; bv[B]=cv[k]                # corner k (added once)
          for(i=1;i<N;i++){ t=i/N; B++; bu[B]=cu[k]+t*ex; bv[B]=cv[k]+t*ey }
        }

        # interior lattice over the (u,v) bounding box, clipped to the quad
        # with an inward margin so it never meets the boundary points
        umin=cu[0];umax=cu[0];vmin=cv[0];vmax=cv[0]
        for(k=1;k<4;k++){
          if(cu[k]<umin)umin=cu[k]; if(cu[k]>umax)umax=cu[k]
          if(cv[k]<vmin)vmin=cv[k]; if(cv[k]>vmax)vmax=cv[k]
        }
        jmin=int(vmin/sv)-1; jmax=int(vmax/sv)+2
        I=0
        for(j=jmin;j<=jmax;j++){
          vv=j*sv
          off=((j%2)!=0 ? su/2 : 0)                    # offset odd rows by su/2
          iimin=int((umin-off)/su)-1; iimax=int((umax-off)/su)+2
          for(ii=iimin;ii<=iimax;ii++){
            uu=off+ii*su
            if(inside(uu,vv,margin)){ I++; iu[I]=uu; iv[I]=vv }
          }
        }

        V=B+I
        printf "%d 2 0 0\n", V > poly
        for(i=1;i<=B;i++) printf "%d %.12g %.12g\n", i,   bu[i], bv[i] > poly
        for(i=1;i<=I;i++) printf "%d %.12g %.12g\n", B+i, iu[i], iv[i] > poly
        printf "%d 0\n", B > poly
        for(i=1;i<=B;i++){ jn=(i%B)+1; printf "%d %d %d\n", i, i, jn > poly }
        printf "0\n" > poly
      }

      printf "%.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g\n", \
             area, x0,y0,z0, ux,uy,uz, vx,vy,vz
    }'
)

[ -n "${AREA:-}" ] || { echo "error: projection failed" >&2; exit 1; }

# --- step 3b: run Triangle for the chosen mode ----------------------------
if [ "$MODE" = lattice ]; then
  # -p triangulate the PSLG, -Y no Steiner points on the boundary, -Q quiet.
  # No -q/-a/-D: any refinement would insert points and destroy the lattice.
  triangle -p -Y -Q "$POLY" >/dev/null
else
  # original quality refinement: -D conforming Delaunay, -q min angle, -a area
  triangle -D -p -q"$QANGLE" -a"$AREA" -Q "$POLY" >/dev/null
fi

NODE="$WORK/mesh.1.node"
ELE="$WORK/mesh.1.ele"

# --- step 4: lift each triangle's three corners back to 3-D ---------------
awk -v ox="$OX" -v oy="$OY" -v oz="$OZ" \
    -v ux="$UX" -v uy="$UY" -v uz="$UZ" \
    -v vx="$VX" -v vy="$VY" -v vz="$VZ" '
  function emit(i,   u,v){
    u=U[i]; v=V[i]
    printf "%.6f %.6f %.6f", ox+u*ux+v*vx, oy+u*uy+v*vy, oz+u*uz+v*vz
  }
  /^[[:space:]]*#/ { next }
  /^[[:space:]]*$/ { next }
  FNR==NR {
    if (nh++ == 0) next
    U[$1]=$2; V[$1]=$3; next
  }
  {
    if (eh++ == 0) next
    emit($2); printf "  "; emit($3); printf "  "; emit($4); printf "\n"
  }
' "$NODE" "$ELE"
