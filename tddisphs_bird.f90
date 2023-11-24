#define ANGDISDISP  Angdisdisp_Bird
#define ANGDISDISPFSC  Angdisdispfsc_Bird
!#define ANGDISDISP  Angdisdisp
!#define ANGDISDISPFSC  Angdisdispfsc


!THIS IS THE MAIN Nikkhoo & Walter [2015] SUBPROGRAM WHICH WILL BE CALLED FROM OUTSIDE THE MODULE:
!function [ue,un,uv]=TDdispHS(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,nu)
SUBROUTINE tddisphs_bird(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, nu, & ! inputs
     & ue, un, uv)                            ! outputs


  ! is a translation from MatLab code to Fortran 90 code of the
  ! "TDdispHS" analytic solution for displacements ("disp") due to a
  ! triangular dislocation ("TD") patch, of uniform Burgers vector,
  ! in an elastic halfspace ("HS").

  ! The original solution is by Nikkhoo & Walter [2015, Geophys. J. Int.].
  ! The translation is by Peter Bird, UCLA, 2016.10, for use in
  ! versions 5+ of NeoKinema (etc.).

  ! Note that in the MatLab original, the coordinates X, Y, Z of the
  ! test point can be arrays.  However, in this Fortran90 version, the
  ! test point is a single point and thus X, Y, Z are scalar REAL*8 variables.

  ! In this file, the original MatLab code lines are commented-out,
  ! and each is followed by its Fortran90 equivalent.


  ! However, it is better to avoid any test points which touch the triangle,
  ! so that these special flags (or numerical abends) are avoided.
  !% TDdispHS 
  !% Calculates displacements associated with a triangular dislocation in an 
  !% elastic half-space.
  !%
  !% TD: Triangular Dislocation
  !% EFCS: Earth-Fixed Coordinate System
  !% TDCS: Triangular Dislocation Coordinate System
  !% ADCS: Angular Dislocation Coordinate System
  !% 
  !% INPUTS
  !% X, Y and Z: 
  !% Coordinates of calculation points in EFCS (East, North, Up). X, Y and Z 
  !% must have the same size. 

  !Note: The comment above that X, Y, and Z "must have the same size" means that
  !      in the MatLab version, each of these is allowed to be an array of dimension 0, 1, 2, or 3,
  !      thus permitting the simultaneous computation of displacements at many test points
  !      with a single function call to TDdispHS.
  !HOWEVER, in this Fortran version X, Y, and Z are only allowed to be simple scalars (real numbers)
  !      giving the 3 coordinates of a single observation point.  To obtain results for other
  !      observation/test points, CALL TDdispHS from within a Fortran DO loop (or nested DO loops). 

  IMPLICIT NONE
  REAL*8, INTENT(IN) :: X, Y, Z

  !% P1,P2 and P3:
  !% Coordinates of TD vertices in EFCS. 
  REAL*8, DIMENSION(3), INTENT(IN) :: P1, P2, P3

  !% Ss, Ds and Ts:
  !% TD slip vector components (Strike-slip, Dip-slip, Tensile-slip).
  REAL*8, INTENT(IN) :: Ss, Ds, Ts

  !% nu:
  !% Poisson's ratio.
  REAL*8, INTENT(IN) :: nu

  !% OUTPUTS
  !% ue, un and uv:
  !% Calculated displacement vector components in EFCS. ue, un and uv have
  !% the same units (e.g., m, or m/s) as Ss, Ds and Ts in the inputs.
  REAL*8, INTENT(OUT) :: ue, un, uv

  REAL*8, DIMENSION(3) :: P1_image, P2_image, P3_image
  REAL*8 :: ueMS,  unMS,  uvMS
  REAL*8 :: ueFSC, unFSC, uvFSC
  REAL*8 :: ueIS,  unIS,  uvIS

  !% Example: Calculate and plot the first component of displacement vector 
  !% on a regular grid.
  !% 
  !% [X,Y,Z] = meshgrid(-3:.02:3,-3:.02:3,-5);
  !% [ue,un,uv] = TDdispHS(X,Y,Z,[-1 0 0],[1 -1 -1],[0 1.5 -2],-1,2,3,.25);
  !% h = surf(X,Y,reshape(ue,size(X)),'edgecolor','none');
  !% view(2)
  !% axis equal
  !% axis tight
  !% set(gcf,'renderer','painters')
  !
  !% Reference journal article: 
  !% Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An analytical, 
  !% artefact-free solution. 
  !% Submitted to Geophysical Journal International 
  !
  !% Copyright (c) 2014 Mehdi Nikkhoo
  !% 
  !% Permission is hereby granted, free of charge, to any person obtaining a 
  !% copy of this software and associated documentation files 
  !% (the "Software"), to deal in the Software without restriction, including 
  !% without limitation the rights to use, copy, modify, merge, publish, 
  !% distribute, sublicense, and/or sell copies of the Software, and to permit
  !% persons to whom the Software is furnished to do so, subject to the 
  !% following conditions:
  !% 
  !% The above copyright notice and this permission notice shall be included 
  !% in all copies or substantial portions of the Software.
  !% 
  !% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
  !% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
  !% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
  !% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
  !% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
  !% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
  !% USE OR OTHER DEALINGS IN THE SOFTWARE.
  !
  !% I appreciate any comments or bug reports.
  !
  !% Mehdi Nikkhoo
  !% created: 2013.1.24
  !% Last modified: 2014.7.30
  !% 
  !% VolcanoTectonics Research Group
  !% Section 2.1, Physics of Earthquakes and Volcanoes
  !% Department 2, Physics of the Earth
  !% Helmholtz Centre Potsdam
  !% German Research Centre for Geosciences (GFZ)
  !% 
  !% email: 
  !% mehdi.nikkhoo@gfz-potsdam.de 
  !% mehdi.nikkhoo@gmail.com
  !
  !if any(Z>0 | P1(3)>0 | P2(3)>0 | P3(3)>0)
  !    error('Half-space solution: Z coordinates must be negative!')
  !end
  IF ((Z > 0.0D0).OR.(P1(3) > 0.0D0).OR.(P2(3) > 0.0D0).OR.(P3(3) > 0.0D0)) THEN
     WRITE (*, "(' ERROR: Entering TDdispHS, all z-coordinates (Z, P1(3), P2(3), P3(3)) must be non-positive.')") 
     !CALL Pause()
     STOP
  END IF

  !X = X(:);
  !Y = Y(:);
  !Z = Z(:);
  !These MatLab statements reshape X, Y, Z into column vectors, regardless of their original size and shape.
  !(Each could be a single point, a row of points, and planar grid of points, or a 3-D grid of points.)

  !P1 = P1(:);
  !P2 = P2(:);
  !P3 = P3(:);
  !N.B. These MatLab statements reshape input vectors P1, P2, P3 to be column vectors,
  !     in case they were accidentally input as row vectors.
  !
  !% Calculate main dislocation contribution to displacements
  ![ueMS,unMS,uvMS] = Tddispfs_Bird(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,nu);
  !print *,p1,p2,p3
  !print *,x,y,z
  !print *,Ss,Ds,Ts
  CALL Tddispfs_Bird(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, nu, & ! inputs
       & ueMS, unMS, uvMS)                      ! outputs
  !print *,ueMS,unMS,uvMS
  !stop
  !% Calculate harmonic function contribution to displacements
  ![ueFSC,unFSC,uvFSC] = Tddisp_Harfunc_Bird(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,nu);
  CALL Tddisp_Harfunc_Bird(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, nu, & ! inputs
       & ueFSC, unFSC, uvFSC)                 ! outputs
  !print *,ueFSC,unFSC,uvFSC
  !% Calculate image dislocation contribution to displacements
  !P1(3) = -P1(3);
  !P2(3) = -P2(3);
  !P3(3) = -P3(3);
  !N.B. I do not approve of changing input parameters
  !    (and with Fortran 90 "INTENT(IN)" this will not be allowed).
  !     Therefore, I take copies, and operate on the copies:
  P1_image(1) = P1(1); P1_image(2) = P1(2); P1_image(3) = -P1(3)
  P2_image(1) = P2(1); P2_image(2) = P2(2); P2_image(3) = -P2(3)
  P3_image(1) = P3(1); P3_image(2) = P3(2); P3_image(3) = -P3(3)

  ![ueIS,unIS,uvIS] = Tddispfs_Bird(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,nu);
  CALL Tddispfs_Bird(X, Y, Z, P1_image, P2_image, P3_image, Ss, Ds, Ts, nu, & ! inputs
       & ueIS, unIS, uvIS)                                        ! outputs
  !print *,ueIS,unIS,uvIS
  !stop
  !if P1(3)==0 && P2(3)==0 && P3(3)==0
  !    uvIS = -uvIS;
  !end
  IF ((P1(3) == 0.0D0).AND.(P2(3) == 0.0D0).AND.(P3(3) == 0.0D0)) THEN
     uvIS = -uvIS
  END IF

  !% Calculate the complete displacement vector components in EFCS
  !ue = ueMS+ueIS+ueFSC;
  !un = unMS+unIS+unFSC;
  !uv = uvMS+uvIS+uvFSC;
  ue = ueMS + ueIS + ueFSC
  un = unMS + unIS + unFSC
  uv = uvMS + uvIS + uvFSC

  !if P1(3)==0 && P2(3)==0 && P3(3)==0
  !    ue = -ue;
  !    un = -un;
  !    uv = -uv;
  IF ((P1(3) == 0.0D0).AND.(P2(3) == 0.0D0).AND.(P3(3) == 0.0D0)) THEN
     ue = -ue
     un = -un
     uv = -uv
  END IF

  !end
END SUBROUTINE TDdispHS_Bird

!----------------------------------------------------------------------------

!function [ue,un,uv]=Tddispfs_Bird(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,nu)
SUBROUTINE Tddispfs_Bird(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, nu, & ! inputs
     & ue, un, uv)                            ! outputs
  !% Tddispfs_Bird 
  !% Calculates displacements associated with a triangular dislocation in an
  !% elastic full-space.

  USE, INTRINSIC :: IEEE_ARITHMETIC ! ONLY used to support setting (and reading) NaN flags in REAL*8 variables.
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: X, Y, Z ! <==== In my version, these are scalars,
  !       and (X, Y, Z) is the single observation point.
  !       However, in the MatLab version, these are (potentially very long) column vectors,
  !       with values rearranged from 0-D, 1-D, 2-D, or 3-D arrays!
  REAL*8, DIMENSION(3), INTENT(IN) :: P1, P2, P3 ! The three vertices of the triangular dislocation.
  REAL*8, INTENT(IN) :: Ss, Ds, Ts ! 3 components of Burger's vector: Strike-slip, Dip-slip, Tensile-slip
  REAL*8, INTENT(IN) :: nu
  REAL*8, INTENT(OUT) :: ue, un, uv ! 3 components of displacement at the observation point (East, North, Vertical/Up).

  INTEGER :: Trimode
  LOGICAL :: casepLog, casenLog, casezLog
  REAL*8 :: A, B, C, bx, by, bz, Fi, na, nb, nc, &
       & u, u1Tn, u2Tn, u3Tn, u1Tp, u2Tp, u3Tp, v, &
       v1Tn, v2Tn, v3Tn, v1Tp, v2Tp, v3Tp, w, w1Tn, w2Tn, w3Tn, w1Tp, w2Tp, w3Tp, &
       & x_, xn, xp, y_, yn, yp, z_, zn, zp
  REAL*8, DIMENSION(3) :: a_, b_, c_  ! N.B. "_" = "lower-case" (because Fortran is NOT case-sensitive).
  REAL*8, DIMENSION(3) :: e12, e13, e23, eY, eZ, P2mP1, P3mP1, p_1, p_2, p_3, Vnorm, Vstrike, Vdip
  REAL*8, DIMENSION(3, 3) :: At, matrix3x3

  !bx = Ts; % Tensile-slip
  !by = Ss; % Strike-slip
  !bz = Ds; % Dip-slip
  bx = Ts
  by = Ss
  bz = Ds

  !% Calculate unit strike, dip and normal to TD vectors: For a horizontal TD 
  !% as an exception, if the normal vector points upward, the strike and dip 
  !% vectors point Northward and Westward, whereas if the normal vector points
  !% downward, the strike and dip vectors point Southward and Westward, 
  !% respectively.
  !Vnorm = cross(P2-P1,P3-P1);
  !Vnorm = Vnorm/norm(Vnorm);
  P2mP1 = P2 - P1 ! for components (1:3)
  P3mP1 = P3 - P1
  CALL Dcross_Bird(P2mP1, P3mP1, Vnorm) ! N.B. Vnorm is actually not normalized yet; see next line:
  CALL DMake_Uvec(Vnorm, Vnorm)

  !eY = [0 1 0]';
  !eZ = [0 0 1]';
  !Vstrike = cross(eZ,Vnorm);
  eY = (/ 0.0D0, 1.0D0, 0.0D0 /)
  eZ = (/ 0.0D0, 0.0D0, 1.0D0 /)
  CALL Dcross_Bird(eZ, Vnorm, Vstrike)

  !if norm(Vstrike)==0
  !    Vstrike = eY*Vnorm(3);
  !    % For horizontal elements in case of half-space calculation!!!
  !    % Correct the strike vector of image dislocation only
  !    if P1(3)>0
  !        Vstrike = -Vstrike;
  !    end
  !end
  !Vstrike = Vstrike/norm(Vstrike);
  !Vdip = cross(Vnorm,Vstrike);
  IF (norm2(Vstrike) == 0.0D0) THEN
     Vstrike = eY * Vnorm(3)
     IF (P1(3) > 0.0D0) THEN
        Vstrike = -Vstrike
     END IF
  END IF
  CALL DMake_Uvec(Vstrike, Vstrike)
  CALL Dcross_Bird(Vnorm, Vstrike, Vdip)

  !% Transform coordinates and slip vector components from EFCS into TDCS
  !p1 = zeros(3,1);
  !p2 = zeros(3,1);
  !p3 = zeros(3,1);
  p_1 = 0.0D0 ! "_" = "lower-case"
  p_2 = 0.0D0 ! N.B. p_1..3 previously specified as REAL*8, DIMENSION(3).  (However, Fortran does not distinguish between row-vectors and column-vectors.)
  p_3 = 0.0D0 ! N.B. Two of these 3 statements are not actually necessary, due to code below.  However, "p_2 = 0.0D0" is necessary because this will not be set below.

  !At = [Vnorm Vstrike Vdip]';
  At(1, 1:3) = Vnorm(1:3)
  At(2, 1:3) = Vstrike(1:3)
  At(3, 1:3) = Vdip(1:3)
  !N.B. The MatLab line has an apostrophe after the brackets, indicating a transpose.    

  ![x,y,z] = Coordtrans_Bird(X'-P2(1),Y'-P2(2),Z'-P2(3),At);
  !N.B. The apostrophes in the MatLab line above convert (very long?) column vectors into (very long?) row vectors for each of X, Y, Z.
  CALL Coordtrans_Bird(X-P2(1), Y-P2(2), Z-P2(3), At, & ! inputs
       & x_, y_, z_)                      ! outputs
  !However, in my Fortran version, each of X, Y, Z is a simple REAL*8 scalar; together, they represent a single observation point.
  !Therefore x_, y_, and z_ are also simple scalars.

  ![p1(1),p1(2),p1(3)] = Coordtrans_Bird(P1(1)-P2(1),P1(2)-P2(2),P1(3)-P2(3),At);
  ![p3(1),p3(2),p3(3)] = Coordtrans_Bird(P3(1)-P2(1),P3(2)-P2(2),P3(3)-P2(3),At);
  CALL Coordtrans_Bird(P1(1)-P2(1), P1(2)-P2(2), P1(3)-P2(3), At, & ! inputs
       & p_1(1), p_1(2), p_1(3))                      ! outputs
  CALL Coordtrans_Bird(P3(1)-P2(1), P3(2)-P2(2), P3(3)-P2(3), At, & ! inputs
       & p_3(1), p_3(2), p_3(3))                      ! outputs

  !% Calculate the unit vectors along TD sides in TDCS
  !e12 = (p2-p1)/norm(p2-p1);
  !e13 = (p3-p1)/norm(p3-p1);
  !e23 = (p3-p2)/norm(p3-p2);
  e12 = p_2 - p_1; CALL DMake_Uvec(e12, e12)
  e13 = p_3 - p_1; CALL DMake_Uvec(e13, e13)
  e23 = p_3 - p_2; CALL DMake_Uvec(e23, e23)

  !% Calculate the TD angles
  !A = acos(e12'*e13);
  !B = acos(-e12'*e23);
  !C = acos(e23'*e13);
  A = ACOS(DOT_PRODUCT(e12, e13))
  B = ACOS(DOT_PRODUCT(-e12, e23))
  C = ACOS(DOT_PRODUCT(e23, e13))
  !print *,'abc',A,B,C
  !% Determine the best artefact-free configuration for each calculation point
  !Trimode = trimodefinder_bird(y,z,x,p1(2:3),p2(2:3),p3(2:3));
  CALL trimodefinder_bird(y_, z_, x_, p_1(2:3), p_2(2:3), p_3(2:3), & ! inputs
       & Trimode)                                    ! output
  !N.B. In the MatLab original, Trimode is an array of numbers the same size/shape as x_ (or y_ or z_);
  !     however, in my Fortran version, it is a single scalar integer.

  !casepLog = Trimode==1;
  !casenLog = Trimode==-1;
  !casezLog = Trimode==0;
  casepLog = (Trimode == 1) ! Note that Fortran .TRUE. and .FALSE. are not necessarily represented the same way as MatLab 1 and 0 !
  casenLog = (Trimode == -1)
  casezLog = (Trimode == 0)
  !N.B. In the MatLab original, casepLog, casenLog, & casexLog all have the same size/shape as Trimode, which 
  !     is the same size/shape as x_ (or y_ or z_).
  !     However, in my Fortran version there is only a single observation point, and so each of these
  !     logical variables is a single T/F value.

  !xp = x(casepLog);
  !yp = y(casepLog);
  !zp = z(casepLog);
  !xn = x(casenLog);
  !yn = y(casenLog);
  !zn = z(casenLog);
  !N.B. The statements above use "logical indexing" in MatLab.  The result is to select all the elements of the
  !     array for which the matching logical element is true, and arrange them into a smaller column vector!
  !     Therefore, these statements divide the lists of observation points into 2 smaller column vectors.
  !     In my Fortran code, there is only a single observation point, so only ONE of (xp, yp, zp) OR (xn, yn, zn)
  !     will be a meaningful vector that will be used in computation!
  IF (casepLog) THEN
     xp = x_
     yp = y_
     zp = z_
  END IF
  IF (casenLog) THEN
     xn = x_
     yn = y_
     zn = z_
  END IF

  !% Configuration I
  !if nnz(casepLog)~=0
  !    % Calculate first angular dislocation contribution
  !    [u1Tp,v1Tp,w1Tp] = Tdsetupd_Bird(xp,yp,zp,A,bx,by,bz,nu,p1,-e13);
  !    % Calculate second angular dislocation contribution
  !    [u2Tp,v2Tp,w2Tp] = Tdsetupd_Bird(xp,yp,zp,B,bx,by,bz,nu,p2,e12);
  !    % Calculate third angular dislocation contribution
  !    [u3Tp,v3Tp,w3Tp] = Tdsetupd_Bird(xp,yp,zp,C,bx,by,bz,nu,p3,e23);
  !end
  IF (casepLog) THEN
     !print *,xp, yp, zp, A, bx, by, bz, nu, p_1
     CALL Tdsetupd_Bird(xp, yp, zp, A, bx, by, bz, nu, p_1, -e13, & ! inputs
          & u1Tp, v1Tp, w1Tp)
     !print *,-e13,u1Tp, v1Tp, w1Tp
     CALL Tdsetupd_Bird(xp, yp, zp, B, bx, by, bz, nu, p_2, e12, & ! inputs
          & u2Tp, v2Tp, w2Tp)
     CALL Tdsetupd_Bird(xp, yp, zp, C, bx, by, bz, nu, p_3, e23, & ! inputs
          & u3Tp, v3Tp, w3Tp)
  END IF

  !% Configuration II
  !if nnz(casenLog)~=0
  !    % Calculate first angular dislocation contribution
  !    [u1Tn,v1Tn,w1Tn] = Tdsetupd_Bird(xn,yn,zn,A,bx,by,bz,nu,p1,e13);
  !    % Calculate second angular dislocation contribution
  !    [u2Tn,v2Tn,w2Tn] = Tdsetupd_Bird(xn,yn,zn,B,bx,by,bz,nu,p2,-e12);
  !    % Calculate third angular dislocation contribution
  !    [u3Tn,v3Tn,w3Tn] = Tdsetupd_Bird(xn,yn,zn,C,bx,by,bz,nu,p3,-e23);
  !end
  IF (caseNLog) THEN
     CALL Tdsetupd_Bird(xn, yn, zn, A, bx, by, bz, nu, p_1, e13, & ! inputs
          & u1Tn, v1Tn, w1Tn)
     CALL Tdsetupd_Bird(xn, yn, zn, B, bx, by, bz, nu, p_2, -e12, & ! inputs
          & u2Tn, v2Tn, w2Tn)
     CALL Tdsetupd_Bird(xn, yn, zn, C, bx, by, bz, nu, p_3, -e23, & ! inputs
          & u3Tn, v3Tn, w3Tn)
    
  END IF

  !% Calculate the "incomplete" displacement vector components in TDCS
  !if nnz(casepLog)~=0
  !    u(casepLog,1) = u1Tp+u2Tp+u3Tp;
  !    v(casepLog,1) = v1Tp+v2Tp+v3Tp;
  !    w(casepLog,1) = w1Tp+w2Tp+w3Tp;
  !end
  !if nnz(casenLog)~=0
  !    u(casenLog,1) = u1Tn+u2Tn+u3Tn;
  !    v(casenLog,1) = v1Tn+v2Tn+v3Tn;
  !    w(casenLog,1) = w1Tn+w2Tn+w3Tn;
  !end
  !if nnz(casezLog)~=0
  !    u(casezLog,1) = nan;
  !    v(casezLog,1) = nan;
  !    w(casezLog,1) = nan;
  !end
  IF (casepLog) THEN
     u = u1Tp + u2Tp + u3Tp;
     v = v1Tp + v2Tp + v3Tp
     w = w1Tp + w2Tp + w3Tp
  END IF
  IF (casenLog) THEN
     u = u1Tn + u2Tn + u3Tn
     v = v1Tn + v2Tn + v3Tn
     w = w1Tn + w2Tn + w3Tn
  END IF
  IF (casezLog) THEN
     u = IEEE_VALUE(u, IEEE_QUIET_NAN)
     v = IEEE_VALUE(v, IEEE_QUIET_NAN)
     w = IEEE_VALUE(w, IEEE_QUIET_NAN)
  END IF

  !% Calculate the Burgers' function contribution corresponding to the TD
  !a = [-x p1(2)-y p1(3)-z];
  !b = [-x -y -z];
  !c = [-x p3(2)-y p3(3)-z];
  !na = sqrt(sum(a.^2,2));
  !nb = sqrt(sum(b.^2,2));
  !nc = sqrt(sum(c.^2,2));
  a_ = (/ -x_, p_1(2)-y_, p_1(3)-z_ /)
  b_ = (/ -x_, -y_, -z_ /)
  c_ = (/ -x_, p_3(2)-y_, p_3(3)-z_ /)
  !MatLab manual: "If A is a matrix, then sum(A,2) is a column vector containing the sum of each row."
  !However, in this Fortran version, a_, b_, c_ have only one row each, so every sum( ,2) is a scalar.
  na = norm2(a_)
  nb = norm2(b_)
  nc = norm2(c_)

  !Fi = -2*atan2((a(:,1).*(b(:,2).*c(:,3)-b(:,3).*c(:,2))-...
  !    a(:,2).*(b(:,1).*c(:,3)-b(:,3).*c(:,1))+...
  !    a(:,3).*(b(:,1).*c(:,2)-b(:,2).*c(:,1))),...
  !    (na.*nb.*nc+sum(a.*b,2).*nc+sum(a.*c,2).*nb+sum(b.*c,2).*na))/4/pi;
  !N.B. The MatLab function atan2() works just like Fortran ATAN2; the first argument is y and the second is x.
  !     Also note that a_, b_, c_ have only one row each in this Fortran version.
  Fi = -2.0D0 * ATAN2((a_(1)*(b_(2)*c_(3) - b_(3)*c_(2)) - &
       &  a_(2)*(b_(1)*c_(3) - b_(3)*c_(1)) + &
       &  a_(3)*(b_(1)*c_(2) - b_(2)*c_(1))), & ! This concludes the "y" argument.
       &  (na*nb*nc + DOT_PRODUCT(a_, b_)*nc + DOT_PRODUCT(a_, c_)*nb + DOT_PRODUCT(b_, c_)*na)) & ! This concludes the "x" argument.
       &  / 4.0D0 / 3.14159265358979D0;

  !% Calculate the complete displacement vector components in TDCS
  !u = bx.*Fi+u;
  !v = by.*Fi+v;
  !w = bz.*Fi+w;
  u = bx * Fi + u
  v = by * Fi + v
  w = bz * Fi + w
  !print *,u,v,w
  !where all variables in the 3 lines above are simple scalars.
  !print *,'uvw',u,v,w
  !% Transform the complete displacement vector components from TDCS into EFCS
  ![ue,un,uv] = Coordtrans_Bird(u,v,w,[Vnorm Vstrike Vdip]);
  matrix3x3(1:3, 1) = Vnorm(1:3)
  matrix3x3(1:3, 2) = Vstrike(1:3)
  matrix3x3(1:3, 3) = Vdip(1:3)
  CALL Coordtrans_Bird(u, v, w, matrix3x3, & ! inputs
       & ue, un, uv)           ! outputs
  !print *,ue,un,uv
END SUBROUTINE Tddispfs_Bird

!----------------------------------------------------------------------------------

!function [ue,un,uv]=Tddisp_Harfunc_Bird(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,nu)
SUBROUTINE Tddisp_Harfunc_Bird(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, nu, & ! inputs
     & ue, un, uv)                            ! outputs

  !% Tddisp_Harfunc_Bird calculates the harmonic function contribution to the
  !% displacements associated with a triangular dislocation in a half-space.
  !% The function cancels the surface normal tractions induced by the main and
  !% image dislocations.

  IMPLICIT NONE
  REAL*8, INTENT(IN) :: X, Y, Z, Ss, Ds, Ts, nu
  REAL*8, DIMENSION(3), INTENT(IN) :: P1, P2, P3
  REAL*8, INTENT(OUT) :: ue, un, uv

  REAL*8 :: bx, by, bz, bXcap, bYcap, bZcap, & ! "cap" = "capital" (because Fortran is not case-sensitive)
       & u1, u2, u3, v1, v2, v3, w1, w2, w3
  REAL*8, DIMENSION(3) ::  eY, eZ, P2mP1, P3mP1, Vdip, Vnorm, Vstrike ! "cap" = "capital" (because Fortran is not case-sensitive)
  REAL*8, DIMENSION(3, 3) :: A

  !bx = Ts; % Tensile-slip
  !by = Ss; % Strike-slip
  !bz = Ds; % Dip-slip
  bx = Ts
  by = Ss
  bz = Ds

  !% Calculate unit strike, dip, and normal-to-TD vectors: For a horizontal TD,
  !% as an exception, if the normal vector points upward, the strike and dip 
  !% vectors point Northward and Westward, whereas if the normal vector points
  !% downward, the strike and dip vectors point Southward and Westward, 
  !% respectively.
  !Vnorm = cross(P2-P1,P3-P1);
  !Vnorm = Vnorm/norm(Vnorm);
  P2mP1 = P2 - P1 ! all 3 components
  P3mP1 = P3 - P1
  CALL Dcross_Bird(P2mP1, P3mP1, Vnorm) ! but, this still needs to be normalized:
  CALL DMake_Uvec(Vnorm, Vnorm)

  !eY = [0 1 0]';
  !eZ = [0 0 1]';
  !Vstrike = cross(eZ,Vnorm);
  eY = (/ 0.0D0, 1.0D0, 0.0D0 /)
  eZ = (/ 0.0D0, 0.0D0, 1.0D0 /)
  CALL Dcross_Bird(eZ, Vnorm, Vstrike)

  !if norm(Vstrike)==0
  !    Vstrike = eY*Vnorm(3);
  !end
  !Vstrike = Vstrike/norm(Vstrike);
  !Vdip = cross(Vnorm,Vstrike);
  IF (norm2(Vstrike) == 0.0D0) THEN
     Vstrike = eY * Vnorm(3)
  END IF
  CALL DMake_Uvec(Vstrike, Vstrike)
  CALL Dcross_Bird(Vnorm, Vstrike, Vdip) ! (apparently, no normalization is needed here)

  !% Transform slip vector components from TDCS into EFCS
  !A = [Vnorm Vstrike Vdip];
  ![bX,bY,bZ] = Coordtrans_Bird(bx,by,bz,A);
  A(1:3, 1) = Vnorm(1:3)
  A(1:3, 2) = Vstrike(1:3)
  A(1:3, 3) = Vdip(1:3)
  CALL Coordtrans_Bird(bx, by, bz, A, &     ! inputs
       & bXcap, bYcap, bZcap) ! outputs

  !% Calculate contribution of angular dislocation pair on each TD side 
  ![u1,v1,w1] = Angsetupfsc_Bird(X,Y,Z,bX,bY,bZ,P1,P2,nu); % Side P1P2
  ![u2,v2,w2] = Angsetupfsc_Bird(X,Y,Z,bX,bY,bZ,P2,P3,nu); % Side P2P3
  ![u3,v3,w3] = Angsetupfsc_Bird(X,Y,Z,bX,bY,bZ,P3,P1,nu); % Side P3P1
  CALL Angsetupfsc_Bird(X, Y, Z, bXcap, bYcap, bZcap, P1, P2, nu, & ! inputs
       & u1, v1, w1)                                 ! outputs
  CALL Angsetupfsc_Bird(X, Y, Z, bXcap, bYcap, bZcap, P2, P3, nu, & ! inputs
       & u2, v2, w2)                                 ! outputs
  CALL Angsetupfsc_Bird(X, Y, Z, bXcap, bYcap, bZcap, P3, P1, nu, & ! inputs
       & u3, v3, w3)                                 ! outputs

  !% Calculate total harmonic function contribution to displacements
  !ue = u1+u2+u3;
  !un = v1+v2+v3;
  !uv = w1+w2+w3;
  ue = u1 + u2 + u3
  un = v1 + v2 + v3
  uv = w1 + w2 + w3

END SUBROUTINE Tddisp_Harfunc_Bird

!-----------------------------------------------------------------------

!function [X1,X2,X3]=Coordtrans_Bird(x1,x2,x3,A)
SUBROUTINE Coordtrans_Bird(x1, x2, x3, A, & ! inputs
     & Xcap1, Xcap2, Xcap3)      ! outputs
  !% Coordtrans_Bird transforms the coordinates of the vectors, from
  !% x1x2x3 coordinate system to X1X2X3 coordinate system. "A" is the
  !% transformation matrix, whose columns e1,e2 and e3 are the unit base 
  !% vectors of the x1x2x3. The coordinates of e1,e2 and e3 in A must be given 
  !% in X1X2X3. The transpose of A (i.e., A') will transform the coordinates 
  !% from X1X2X3 into x1x2x3.

  IMPLICIT NONE
  REAL*8, INTENT(IN) :: x1, x2, x3
  REAL*8, DIMENSION(3, 3) :: A
  REAL*8, INTENT(OUT) :: Xcap1, Xcap2, Xcap3 ! N.B. "cap" = "capital" added because Fortran does not distinguish upper- from lower-case.

  !x1 = x1(:);
  !x2 = x2(:);
  !x3 = x3(:);
  !N.B. The lines above re-arrange the input vectors (lists of points to be transformed) into column vectors.
  !However, in my Fortran version there will only be one point transformed per CALL, so x1, x2, x3 are simple scalars.

  !r = A*[x1';x2';x3'];
  !X1 = r(1,:)';
  !X2 = r(2,:)';
  !X3 = r(3,:)';
  Xcap1 = A(1, 1) * x1 + A(1, 2) * x2 + A(1, 3) * x3
  Xcap2 = A(2, 1) * x1 + A(2, 2) * x2 + A(2, 3) * x3
  Xcap3 = A(3, 1) * x1 + A(3, 2) * x2 + A(3, 3) * x3

END SUBROUTINE Coordtrans_Bird

!-------------------------------------------------------------------------

!function [trimode]=trimodefinder_bird(x,y,z,p1,p2,p3)
SUBROUTINE trimodefinder_bird(x, y, z, p1, p2, p3, & ! inputs
     & trimode)               ! output

  !% trimodefinder_bird calculates the normalized barycentric coordinates of 
  !% the points with respect to the TD vertices and specifies the appropriate
  !% artefact-free configuration of the angular dislocations for the 
  !% calculations. The input matrices x, y and z share the same size and
  !% correspond to the y, z and x coordinates in the TDCS, respectively. p1,
  !% p2 and p3 are two-component matrices representing the y and z coordinates
  !% of the TD vertices in the TDCS, respectively.
  !% The components of the output (trimode) corresponding to each calculation 
  !% points, are 1 for the first configuration, -1 for the second 
  !% configuration and 0 for the calculation point that lie on the TD sides.

  !An important distinction is that this Fortran version of the code only
  !operates on one single test point; thus, x, y, z are simple scalars, never arrays.

  IMPLICIT NONE
  REAL*8, INTENT(IN) :: x, y, z
  REAL*8, DIMENSION(2), INTENT(IN) :: p1, p2, p3
  INTEGER, INTENT(OUT) :: trimode

  REAL*8 :: a, b, c

  !x = x(:);
  !y = y(:);
  !z = z(:);
  !The MatLab lines above convert the x, y, and z arrays to (long?) column vectors.
  !However, in this Fortran version, x, y, and z are simple scalar coordinates of a single test point.
  !Therefore, the output "trimode" will also be a simple integer scalar, not an array.

  !a = ((p2(2)-p3(2)).*(x-p3(1))+(p3(1)-p2(1)).*(y-p3(2)))./...
  !    ((p2(2)-p3(2)).*(p1(1)-p3(1))+(p3(1)-p2(1)).*(p1(2)-p3(2)));
  !b = ((p3(2)-p1(2)).*(x-p3(1))+(p1(1)-p3(1)).*(y-p3(2)))./...
  !    ((p2(2)-p3(2)).*(p1(1)-p3(1))+(p3(1)-p2(1)).*(p1(2)-p3(2)));
  !c = 1-a-b;
  a = ((p2(2)-p3(2))*(x-p3(1))+(p3(1)-p2(1))*(y-p3(2)))/  &
       & ((p2(2)-p3(2))*(p1(1)-p3(1))+(p3(1)-p2(1))*(p1(2)-p3(2))) 
  b = ((p3(2)-p1(2))*(x-p3(1))+(p1(1)-p3(1))*(y-p3(2)))/  &
       & ((p2(2)-p3(2))*(p1(1)-p3(1))+(p3(1)-p2(1))*(p1(2)-p3(2)))
  c = 1.0D0 - a - b;

  !trimode = ones(length(x),1);
  !trimode(a<=0 & b>c & c>a) = -1;
  !trimode(b<=0 & c>a & a>b) = -1;
  !trimode(c<=0 & a>b & b>c) = -1;
  !trimode(a==0 & b>=0 & c>=0) = 0;
  !trimode(a>=0 & b==0 & c>=0) = 0;
  !trimode(a>=0 & b>=0 & c==0) = 0;
  !trimode(trimode==0 & z~=0) = 1;
  trimode = 1 
  IF ((a <= 0.0D0) .AND. (b > c) .AND. (c > a)) trimode = -1
  IF ((b <= 0.0D0) .AND. (c > a) .AND. (a > b)) trimode = -1
  IF ((c <= 0.0D0) .AND. (a > b) .AND. (b > c)) trimode = -1
  IF ((a == 0.0D0) .AND. (b >= 0.0D0) .AND. (c >= 0.0D0)) trimode = 0
  IF ((a >= 0.0D0) .AND. (b == 0.0D0) .AND. (c >= 0.0D0)) trimode = 0
  IF ((a >= 0.0D0) .AND. (b >= 0.0D0) .AND. (c == 0.0D0)) trimode = 0
  IF ((trimode == 0) .AND. (z /= 0.0D0)) trimode = 1

END SUBROUTINE trimodefinder_bird
!--------------------------------------------------------------------

!function [u,v,w]=Tdsetupd_Bird(x,y,z,alpha,bx,by,bz,nu,TriVertex,SideVec)
SUBROUTINE Tdsetupd_Bird(x, y, z, alpha, bx, by, bz, nu, TriVertex, SideVec, & ! inputs
     & u, v, w)                                              ! outputs

  !% Tdsetupd_Bird transforms coordinates of the calculation points as well as 
  !% slip vector components from ADCS into TDCS. It then calculates the 
  !% displacements in ADCS and transforms them into TDCS.

  !In this Fortran version, only a single test point (x, y, z) is transformed;
  !we do not allow x, y, or z to be an array (as in the MatLab version).
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: alpha, bx, by, bz, nu, x, y, z
  REAL*8, DIMENSION(3), INTENT(IN) :: SideVec, TriVertex
  REAL*8, INTENT(OUT) :: u, v, w

  REAL*8 :: by1, bz1, v0, w0, y1, z1
  REAL*8, DIMENSION(2) :: r1, r2, r3
  REAL*8, DIMENSION(2, 2) :: A

  !% Transformation matrix
  !A = [[SideVec(3);-SideVec(2)] SideVec(2:3)]';
  A(1, 1) = SideVec(3)
  A(1, 2) = -SideVec(2)
  A(2, 1) = SideVec(2)
  A(2, 2) = SideVec(3)

  !% Transform coordinates of the calculation points from TDCS into ADCS
  !r1 = A*[y'-TriVertex(2);z'-TriVertex(3)];
  !y1 = r1(1,:)';
  !z1 = r1(2,:)';
  r1(1) = A(1, 1) * (y - TriVertex(2)) + A(1, 2) * (z -TriVertex(3))
  r1(2) = A(2, 1) * (y - TriVertex(2)) + A(2, 2) * (z -TriVertex(3))
  y1 = r1(1)
  z1 = r1(2)

  !% Transform the in-plane slip vector components from TDCS into ADCS
  !r2 = A*[by;bz];
  !by1 = r2(1,:)';
  !bz1 = r2(2,:)';
  r2(1) = A(1, 1) * by + A(1, 2) * bz
  r2(2) = A(2, 1) * by + A(2, 2) * bz
  by1 = r2(1)
  bz1 = r2(2)

  !% Calculate displacements associated with an angular dislocation in ADCS
  ![u,v0,w0] = Angdisdisp_Bird(x,y1,z1,-pi+alpha,bx,by1,bz1,nu);
  !print *,'bin',x,y1,z1,-3.14159265358979D0+alpha,bx,by1,bz1,nu
  CALL ANGDISDISP(x, y1, z1, -3.14159265358979D0+alpha, bx, by1, bz1, nu, & ! inputs
       & u, v0, w0)                                ! outputs
  !print *,'bout',v0,w0
  !% Transform displacements from ADCS into TDCS
  !r3 = A'*[v0';w0'];
  !v = r3(1,:)';
  !w = r3(2,:)';
  r3(1) = A(1, 1) * v0 + A(2, 1) * w0
  r3(2) = A(1, 2) * v0 + A(2, 2) * w0
  v = r3(1)
  w = r3(2)

END SUBROUTINE Tdsetupd_Bird
!---------------------------------------------------------------

!function [ue,un,uv]=Angsetupfsc_Bird(X,Y,Z,bX,bY,bZ,PA,PB,nu)
SUBROUTINE Angsetupfsc_Bird(X, Y, Z, bX, bY, bZ, PA, PB, nu, & ! inputs
     & ue, un, uv)                        ! outputs

  !% Angsetupfsc_Bird calculates the Free Surface Correction to displacements 
  !% associated with angular dislocation pair on each TD side.

  !N.B. This difference between the MatLab original and this Fortran version is
  !     that here X, Y, and Z are scalars (not arrays), so only one observation point
  !     is processed by this subprogram, and only a single vector (ue, un, uv)
  !     is returned.
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: bX, bY, bZ, nu, X, Y, Z
  REAL*8, DIMENSION(3), INTENT(IN) :: PA, PB
  REAL*8, INTENT(OUT) :: ue, un, uv

  REAL*8 :: b1, b2, b3, beta, v1, v1A, v1B, v2, v2A, v2B, v3, v3A, v3B, y1A, y1AB, y1B, y2A, y2AB, y2B, y3A, y3AB, y3B
  REAL*8, DIMENSION(3) :: ey1, ey2, ey3, eZ, SideVec, SideUvec
  REAL*8, DIMENSION(3, 3) :: A, Atranspose
  REAL*8, PARAMETER :: eps = 5.0D-15 ! a crude approximation of the MatLab "eps" constant.
  LOGICAL :: I

  !% Calculate TD side vector and the angle of the angular dislocation pair
  !SideVec = PB-PA;
  !eZ = [0 0 1]';
  !beta = acos(-SideVec'*eZ/norm(SideVec));
  SideVec = PB - PA ! all 3 components
  CALL DMake_Uvec(SideVec, SideUvec)
  eZ = (/ 0.0D0, 0.0D0, 1.0D0 /)
  beta = ACOS(DOT_PRODUCT(-SideUvec, eZ))
  !GPBhere     

  !if abs(beta)<eps || abs(pi-beta)<eps  ! N.B. "||" is the "short-circuit OR" in MatLab.
  !    ue = zeros(length(X),1);
  !    un = zeros(length(X),1);
  !    uv = zeros(length(X),1);
  !else
  !    ey1 = [SideVec(1:2);0];
  !    ey1 = ey1/norm(ey1);
  !    ey3 = -eZ;
  !    ey2 = cross(ey3,ey1);
  !    A = [ey1,ey2,ey3]; % Transformation matrix
  IF ((ABS(beta) < eps).OR.(ABS(3.14159265358979D0 - beta) < eps)) THEN
     ue = 0.0D0
     un = 0.0D0
     uv = 0.0D0
  ELSE
     ey1 = (/ SideVec(1), SideVec(2), 0.0D0 /)
     CALL DMake_Uvec(ey1, ey1)
     ey3 = -eZ ! all 3 components
     CALL Dcross_Bird(ey3, ey1, ey2)
     A(1:3, 1) = ey1(1:3)
     A(1:3, 2) = ey2(1:3)
     A(1:3, 3) = ey3(1:3)
!  END IF THIS WAS WRONG

  !    % Transform coordinates from EFCS to the first ADCS
  !    [y1A,y2A,y3A] = Coordtrans_Bird(X-PA(1),Y-PA(2),Z-PA(3),A);
  !    % Transform coordinates from EFCS to the second ADCS
  !    [y1AB,y2AB,y3AB] = Coordtrans_Bird(SideVec(1),SideVec(2),SideVec(3),A);
  !    y1B = y1A-y1AB;
  !    y2B = y2A-y2AB;
  !    y3B = y3A-y3AB;
  CALL Coordtrans_Bird(X - PA(1), Y - PA(2), Z - PA(3), A, & ! inputs
       & y1A, y2A, y3A)                        ! outputs
  CALL Coordtrans_Bird(SideVec(1), SideVec(2), SideVec(3), A, & ! inputs
       & y1AB, y2AB, y3AB)                        ! outputs
  y1B = y1A - y1AB
  y2B = y2A - y2AB
  y3B = y3A - y3AB

  !    % Transform slip vector components from EFCS to ADCS
  !    [b1,b2,b3] = Coordtrans_Bird(bX,bY,bZ,A);
  CALL Coordtrans_Bird(bX, bY, bZ, A, & ! inputs
       & b1, b2, b3)      ! outputs

  !    % Determine the best artefact-free configuration for the calculation
  !    % points near the free furface
  !    I = (beta*y1A)>=0;
  I = ((beta * y1A) >= 0.0D0)

  !    % Configuration I
  !    [v1A(I),v2A(I),v3A(I)] = Angdisdispfsc_Bird(y1A(I),y2A(I),y3A(I),...
  !        -pi+beta,b1,b2,b3,nu,-PA(3));
  !    [v1B(I),v2B(I),v3B(I)] = Angdisdispfsc_Bird(y1B(I),y2B(I),y3B(I),...
  !        -pi+beta,b1,b2,b3,nu,-PB(3));
  !    % Configuration II
  !    [v1A(~I),v2A(~I),v3A(~I)] = Angdisdispfsc_Bird(y1A(~I),y2A(~I),y3A(~I),...
  !        beta,b1,b2,b3,nu,-PA(3));
  !    [v1B(~I),v2B(~I),v3B(~I)] = Angdisdispfsc_Bird(y1B(~I),y2B(~I),y3B(~I),...
  !        beta,b1,b2,b3,nu,-PB(3));
  IF (I) THEN
     CALL ANGDISDISPFSC(y1A, y2A, y3A, -3.14159265358979D0+beta, b1, b2, b3, nu, -PA(3), & ! inputs
          & v1A, v2A, v3A)                                                     ! outputs
     CALL ANGDISDISPFSC(y1B, y2B, y3B, -3.14159265358979D0+beta, b1, b2, b3, nu, -PB(3), & ! inputs
          & v1B, v2B, v3B)                                                     ! outputs
  ELSE ! (~I), or .NOT.I
     CALL ANGDISDISPFSC(y1A, y2A, y3A, beta, b1, b2, b3, nu, -PA(3), & ! inputs
          & v1A, v2A, v3A)                                 ! outputs
     CALL ANGDISDISPFSC(y1B, y2B, y3B, beta, b1, b2, b3, nu, -PB(3), & ! inputs
          & v1B, v2B, v3B)                                 ! outputs            
  END IF

  !    % Calculate total Free Surface Correction to displacements in ADCS
  !    v1 = v1B-v1A;
  !    v2 = v2B-v2A;
  !    v3 = v3B-v3A;
  v1 = v1B - v1A
  v2 = v2B - v2A
  v3 = v3B - v3A

  !    % Transform total Free Surface Correction to displacements from ADCS 
  !    % to EFCS
  !    [ue,un,uv] = Coordtrans_Bird(v1,v2,v3,A');
  Atranspose = TRANSPOSE(A)
  CALL Coordtrans_Bird(v1, v2, v3, Atranspose, & ! inputs
       & ue, un, uv)               ! outputs
end IF
END SUBROUTINE Angsetupfsc_Bird

!----------------------------------------------------------------------------

!function [u,v,w]=Angdisdisp_Bird(x,y,z,alpha,bx,by,bz,nu)
SUBROUTINE Angdisdisp_Bird(x, y, z, alpha, bx, by, bz, nu, & ! inputs
     & u, v, w) ! outputs

  !% Angdisdisp_Bird calculates the "incomplete" displacements (without the 
  !% Burgers' function contribution) associated with an angular dislocation in
  !% an elastic full-space.

  !Note that the orginal MatLab version allows x, y, and z to be arrays of test points
  !(in 0-D, 1-D, 2-D, or 3-D); however, in this Fortran version there is only a single
  !test point at (x, y, z), and each of these is a simple REAL*8 scalar number.

  IMPLICIT NONE
  REAL*8, PARAMETER :: pi = 3.14159265358979D0

  REAL*8, INTENT(IN) :: x, y, z, alpha, bx, by, bz, nu
  REAL*8, INTENT(OUT) :: u, v, w

  REAL*8 :: cosA, eta, r, sinA, ux, uy, uz, vx, vy, vz, wx, wy, wz, zz, zeta

  !cosA = cos(alpha);
  !sinA = sin(alpha);
  !eta = y*cosA-z*sinA;
  !zeta = y*sinA+z*cosA;
  !r = sqrt(x.^2+y.^2+z.^2);
  cosA = COS(alpha)
  sinA = SIN(alpha)
  eta = y * cosA - z * sinA
  zeta = y * sinA + z * cosA
  r = SQRT(x**2 + y**2 + z**2)

  !% Avoid complex results for the logarithmic terms
  !zeta(zeta>r) = r(zeta>r);
  !z(z>r) = r(z>r);
  IF (zeta > r) zeta = r
  !In Fortran, we have to avoid changing an input argument of type INTENT(IN), so I substitute a copy: zz.
  IF (z > r) THEN
     zz = r ! apply the limit
  ELSE
     zz = z ! just copy the value, without any change
  END IF

  !ux = bx/8/pi/(1-nu)*(x.*y./r./(r-z)-x.*eta./r./(r-zeta));
  !vx = bx/8/pi/(1-nu)*(eta*sinA./(r-zeta)-y.*eta./r./(r-zeta)+...
  !    y.^2./r./(r-z)+(1-2*nu)*(cosA*log(r-zeta)-log(r-z)));
  !wx = bx/8/pi/(1-nu)*(eta*cosA./(r-zeta)-y./r-eta.*z./r./(r-zeta)-...
  !    (1-2*nu)*sinA*log(r-zeta));
  ux = bx/8.0D0/pi/(1.0D0-nu) * (x*y/r/(r-zz)-x*eta/r/(r-zeta))
  vx = bx/8.0D0/pi/(1.0D0-nu) * (eta*sinA/(r-zeta)-y*eta/r/(r-zeta) + &
       & y**2/r/(r-zz) + (1.0D0-2.0D0*nu) * (cosA*log(r-zeta)-log(r-zz)))
  wx = bx/8.0D0/pi/(1.0D0-nu) * (eta*cosA/(r-zeta)-y/r-eta*zz/r/(r-zeta)- &
       & (1.0D0-2.0D0*nu) * sinA * log(r-zeta));

  !uy = by/8/pi/(1-nu)*(x.^2*cosA./r./(r-zeta)-x.^2./r./(r-z)-...
  !    (1-2*nu)*(cosA*log(r-zeta)-log(r-z)));
  !vy = by*x/8/pi/(1-nu).*(y.*cosA./r./(r-zeta)-...
  !    sinA*cosA./(r-zeta)-y./r./(r-z));
  !wy = by*x/8/pi/(1-nu).*(z*cosA./r./(r-zeta)-cosA^2./(r-zeta)+1./r);
  uy = by/8.0D0/pi/(1.0D0-nu) * (x**2*cosA/r/(r-zeta) - x**2/r/(r-zz) - &
       & (1.0D0-2.0D0*nu) * (cosA*log(r-zeta) - log(r-zz)))
  vy = by*x/8.0D0/pi/(1.0D0-nu) * (y*cosA/r/(r-zeta) - &
       & sinA*cosA/(r-zeta) - y/r/(r-zz))
  wy = by*x/8.0D0/pi/(1.0D0-nu) * (zz*cosA/r/(r-zeta) - cosA**2/(r-zeta)+1.0D0/r)

  !uz = bz*sinA/8/pi/(1-nu).*((1-2*nu)*log(r-zeta)-x.^2./r./(r-zeta));
  !vz = bz*x*sinA/8/pi/(1-nu).*(sinA./(r-zeta)-y./r./(r-zeta));
  !wz = bz*x*sinA/8/pi/(1-nu).*(cosA./(r-zeta)-z./r./(r-zeta));
  uz = bz*sinA/8.0D0/pi/(1.0D0-nu) * ((1.0D0-2.0D0*nu)*log(r-zeta)-x**2/r/(r-zeta))
  vz = bz*x*sinA/8.0D0/pi/(1.0D0-nu) * (sinA/(r-zeta)-y/r/(r-zeta))
  wz = bz*x*sinA/8.0D0/pi/(1.0D0-nu)*(cosA/(r-zeta)-zz/r/(r-zeta))

  !u = ux+uy+uz;
  !v = vx+vy+vz;
  !w = wx+wy+wz;
  u = ux + uy + uz
  v = vx + vy + vz
  w = wx + wy + wz

END SUBROUTINE Angdisdisp_Bird

!----------------------------------------------------------------------------

!function [v1 v2 v3] = Angdisdispfsc_Bird(y1,y2,y3,beta,b1,b2,b3,nu,a)
SUBROUTINE Angdisdispfsc_Bird(y1, y2, y3, beta, b1, b2, b3, nu, a, & ! inputs
     & v1, v2, v3)                            ! outputs

  !% Angdisdispfsc_Bird calculates the harmonic function contribution to the 
  !% displacements associated with an angular dislocation in an elastic 
  !% half-space.

  !Note that the orginal MatLab version allows y1, y2, and y3 to be arrays of test points
  !(in 0-D, 1-D, 2-D, or 3-D); however, in this Fortran version there is only a single
  !test point at (y1, y2, y3), and each of these is a simple REAL*8 scalar number.

  IMPLICIT NONE
  REAL*8, PARAMETER :: pi = 3.14159265358979D0
  REAL*8, INTENT(IN) :: y1, y2, y3, beta, b1, b2, b3, nu, a
  REAL*8, INTENT(OUT) :: v1, v2, v3

  REAL*8 :: cosB, cotB, Fib, r2b, rb, sinB, &
       & v1cb1, v1cb2, v1cb3, v2cb1, v2cb2, v2cb3, v3cb1, v3cb2, v3cb3, &
       & y3b, z1b, z3b

  !sinB = sin(beta);
  !cosB = cos(beta);
  !cotB = cot(beta);
  !y3b = y3+2*a;
  !z1b = y1*cosB+y3b*sinB;
  !z3b = -y1*sinB+y3b*cosB;
  !r2b = y1.^2+y2.^2+y3b.^2;
  !rb = sqrt(r2b);
  sinB = SIN(beta)
  cosB = COS(beta)
  cotB = 1.00D0 / TAN(beta)
  y3b = y3 + 2.0D0 * a
  z1b = y1 * cosB + y3b * sinB
  z3b = -y1 * sinB + y3b * cosB
  r2b = y1**2 + y2**2 + y3b**2
  rb = SQRT(r2b)

  !Fib = 2*atan(-y2./(-(rb+y3b)*cot(beta/2)+y1)); % The Burgers' function
  Fib = 2.0D0 * ATAN(-y2 / (-(rb + y3b) / TAN(beta / 2.0D0) + y1)) ! The Burgers' function

  !v1cb1 = b1/4/pi/(1-nu)*(-2*(1-nu)*(1-2*nu)*Fib*cotB.^2+(1-2*nu)*y2./...
  !    (rb+y3b).*((1-2*nu-a./rb)*cotB-y1./(rb+y3b).*(nu+a./rb))+(1-2*nu).*...
  !    y2.*cosB*cotB./(rb+z3b).*(cosB+a./rb)+a*y2.*(y3b-a)*cotB./rb.^3+y2.*...
  !    (y3b-a)./(rb.*(rb+y3b)).*(-(1-2*nu)*cotB+y1./(rb+y3b).*(2*nu+a./rb)+...
  !    a*y1./rb.^2)+y2.*(y3b-a)./(rb.*(rb+z3b)).*(cosB./(rb+z3b).*((rb*...
  !    cosB+y3b).*((1-2*nu)*cosB-a./rb).*cotB+2*(1-nu)*(rb*sinB-y1)*cosB)-...
  !    a.*y3b*cosB*cotB./rb.^2));
  v1cb1 = b1/4.0D0/pi/(1.0D0-nu)*(-2.0D0*(1.0D0-nu)*(1.0D0-2.0D0*nu)*Fib*cotB**2+(1.0D0-2.0D0*nu)*y2/ &
       & (rb+y3b)*((1.0D0-2.0D0*nu-a/rb)*cotB-y1/(rb+y3b)*(nu+a/rb))+(1.0D0-2.0D0*nu)*               &
       & y2*cosB*cotB/(rb+z3b)*(cosB+a/rb)+a*y2*(y3b-a)*cotB/rb**3+y2*                               &
       & (y3b-a)/(rb*(rb+y3b))*(-(1.0D0-2.0D0*nu)*cotB+y1/(rb+y3b)*(2.0D0*nu+a/rb)+                  &
       & a*y1/rb**2)+y2*(y3b-a)/(rb*(rb+z3b))*(cosB/(rb+z3b)*((rb*                                   &
       & cosB+y3b)*((1.0D0-2.0D0*nu)*cosB-a/rb)*cotB+2.0D0*(1.0D0-nu)*(rb*sinB-y1)*cosB)-            &
       & a*y3b*cosB*cotB/rb**2));

  !v2cb1 = b1/4/pi/(1-nu)*((1-2*nu)*((2*(1-nu)*cotB^2-nu)*log(rb+y3b)-(2*...
  !    (1-nu)*cotB^2+1-2*nu)*cosB*log(rb+z3b))-(1-2*nu)./(rb+y3b).*(y1*...
  !    cotB.*(1-2*nu-a./rb)+nu*y3b-a+y2.^2./(rb+y3b).*(nu+a./rb))-(1-2*...
  !    nu).*z1b*cotB./(rb+z3b).*(cosB+a./rb)-a*y1.*(y3b-a)*cotB./rb.^3+...
  !    (y3b-a)./(rb+y3b).*(-2*nu+1./rb.*((1-2*nu).*y1*cotB-a)+y2.^2./(rb.*...
  !    (rb+y3b)).*(2*nu+a./rb)+a*y2.^2./rb.^3)+(y3b-a)./(rb+z3b).*(cosB^2-...
  !    1./rb.*((1-2*nu).*z1b*cotB+a*cosB)+a*y3b.*z1b*cotB./rb.^3-1./(rb.*...
  !    (rb+z3b)).*(y2.^2*cosB^2-a*z1b*cotB./rb.*(rb*cosB+y3b))));
  v2cb1 = b1/4.0D0/pi/(1.0D0-nu)*((1.0D0-2.0D0*nu)*((2.0D0*(1.0D0-nu)*cotB**2-nu)*log(rb+y3b)-(2.0D0* &
       & (1.0D0-nu)*cotB**2+1.0D0-2.0D0*nu)*cosB*log(rb+z3b))-(1.0D0-2.0D0*nu)/(rb+y3b)*(y1*         &
       & cotB*(1.0D0-2.0D0*nu-a/rb)+nu*y3b-a+y2**2/(rb+y3b)*(nu+a/rb))-(1.0D0-2.0D0*                 &
       & nu)*z1b*cotB/(rb+z3b)*(cosB+a/rb)-a*y1*(y3b-a)*cotB/rb**3+                                  &
       & (y3b-a)/(rb+y3b)*(-2.0D0*nu+1.0D0/rb*((1.0D0-2.0D0*nu)*y1*cotB-a)+y2**2/(rb*                &
       & (rb+y3b))*(2.0D0*nu+a/rb)+a*y2**2/rb**3)+(y3b-a)/(rb+z3b)*(cosB**2-                         &
       & 1.0D0/rb*((1.0D0-2.0D0*nu)*z1b*cotB+a*cosB)+a*y3b*z1b*cotB/rb**3-1.0D0/(rb*                 &
       & (rb+z3b))*(y2**2*cosB**2-a*z1b*cotB/rb*(rb*cosB+y3b))))

  !v3cb1 = b1/4/pi/(1-nu)*(2*(1-nu)*(((1-2*nu)*Fib*cotB)+(y2./(rb+y3b).*(2*...
  !    nu+a./rb))-(y2*cosB./(rb+z3b).*(cosB+a./rb)))+y2.*(y3b-a)./rb.*(2*...
  !    nu./(rb+y3b)+a./rb.^2)+y2.*(y3b-a)*cosB./(rb.*(rb+z3b)).*(1-2*nu-...
  !    (rb*cosB+y3b)./(rb+z3b).*(cosB+a./rb)-a*y3b./rb.^2));
  v3cb1 = b1/4.0D0/pi/(1.0D0-nu)*(2.0D0*(1.0D0-nu)*(((1.0D0-2.0D0*nu)*Fib*cotB)+(y2/(rb+y3b)*(2*  &
       & nu+a/rb))-(y2*cosB/(rb+z3b)*(cosB+a/rb)))+y2*(y3b-a)/rb*(2*                             &
       & nu/(rb+y3b)+a/rb**2)+y2*(y3b-a)*cosB/(rb*(rb+z3b))*(1.0D0-2.0D0*nu-                     &
       & (rb*cosB+y3b)/(rb+z3b)*(cosB+a/rb)-a*y3b/rb**2))

  !v1cb2 = b2/4/pi/(1-nu)*((1-2*nu)*((2*(1-nu)*cotB^2+nu)*log(rb+y3b)-(2*...
  !    (1-nu)*cotB^2+1)*cosB*log(rb+z3b))+(1-2*nu)./(rb+y3b).*(-(1-2*nu).*...
  !    y1*cotB+nu*y3b-a+a*y1*cotB./rb+y1.^2./(rb+y3b).*(nu+a./rb))-(1-2*...
  !    nu)*cotB./(rb+z3b).*(z1b*cosB-a*(rb*sinB-y1)./(rb*cosB))-a*y1.*...
  !    (y3b-a)*cotB./rb.^3+(y3b-a)./(rb+y3b).*(2*nu+1./rb.*((1-2*nu).*y1*...
  !    cotB+a)-y1.^2./(rb.*(rb+y3b)).*(2*nu+a./rb)-a*y1.^2./rb.^3)+(y3b-a)*...
  !    cotB./(rb+z3b).*(-cosB*sinB+a*y1.*y3b./(rb.^3*cosB)+(rb*sinB-y1)./...
  !    rb.*(2*(1-nu)*cosB-(rb*cosB+y3b)./(rb+z3b).*(1+a./(rb*cosB)))));
  v1cb2 = b2/4.0D0/pi/(1.0D0-nu)*((1.0D0-2.0D0*nu)*((2.0D0*(1.0D0-nu)*cotB**2+nu)*log(rb+y3b)-(2* &
       & (1.0D0-nu)*cotB**2+1)*cosB*log(rb+z3b))+(1.0D0-2.0D0*nu)/(rb+y3b)*(-(1.0D0-2.0D0*nu)*   &
       & y1*cotB+nu*y3b-a+a*y1*cotB/rb+y1**2/(rb+y3b)*(nu+a/rb))-(1.0D0-2.0D0*                   &
       & nu)*cotB/(rb+z3b)*(z1b*cosB-a*(rb*sinB-y1)/(rb*cosB))-a*y1*                             &
       & (y3b-a)*cotB/rb**3+(y3b-a)/(rb+y3b)*(2.0D0*nu+1.0D0/rb*((1.0D0-2.0D0*nu)*y1*            &
       & cotB+a)-y1**2/(rb*(rb+y3b))*(2.0D0*nu+a/rb)-a*y1**2/rb**3)+(y3b-a)*                     &
       & cotB/(rb+z3b)*(-cosB*sinB+a*y1*y3b/(rb**3*cosB)+(rb*sinB-y1)/                           &
       & rb*(2.0D0*(1.0D0-nu)*cosB-(rb*cosB+y3b)/(rb+z3b)*(1.0D0+a/(rb*cosB)))))

  !v2cb2 = b2/4/pi/(1-nu)*(2*(1-nu)*(1-2*nu)*Fib*cotB.^2+(1-2*nu)*y2./...
  !    (rb+y3b).*(-(1-2*nu-a./rb)*cotB+y1./(rb+y3b).*(nu+a./rb))-(1-2*nu)*...
  !    y2*cotB./(rb+z3b).*(1+a./(rb*cosB))-a*y2.*(y3b-a)*cotB./rb.^3+y2.*...
  !    (y3b-a)./(rb.*(rb+y3b)).*((1-2*nu)*cotB-2*nu*y1./(rb+y3b)-a*y1./rb.*...
  !    (1./rb+1./(rb+y3b)))+y2.*(y3b-a)*cotB./(rb.*(rb+z3b)).*(-2*(1-nu)*...
  !    cosB+(rb*cosB+y3b)./(rb+z3b).*(1+a./(rb*cosB))+a*y3b./(rb.^2*cosB)));
  v2cb2 = b2/4.0D0/pi/(1.0D0-nu)*(2.0D0*(1.0D0-nu)*(1.0D0-2.0D0*nu)*Fib*cotB**2+(1.0D0-2.0D0*nu)*y2/  &
       & (rb+y3b)*(-(1.0D0-2.0D0*nu-a/rb)*cotB+y1/(rb+y3b)*(nu+a/rb))-(1.0D0-2.0D0*nu)*              &
       & y2*cotB/(rb+z3b)*(1.0D0+a/(rb*cosB))-a*y2*(y3b-a)*cotB/rb**3+y2*                            &
       & (y3b-a)/(rb*(rb+y3b))*((1.0D0-2.0D0*nu)*cotB-2.0D0*nu*y1/(rb+y3b)-a*y1/rb*                  &
       & (1.0D0/rb+1.0D0/(rb+y3b)))+y2*(y3b-a)*cotB/(rb*(rb+z3b))*(-2.0D0*(1.0D0-nu)*                &
       & cosB+(rb*cosB+y3b)/(rb+z3b)*(1.0D0+a/(rb*cosB))+a*y3b/(rb**2*cosB)))

  !v3cb2 = b2/4/pi/(1-nu)*(-2*(1-nu)*(1-2*nu)*cotB*(log(rb+y3b)-cosB*...
  !    log(rb+z3b))-2*(1-nu)*y1./(rb+y3b).*(2*nu+a./rb)+2*(1-nu)*z1b./(rb+...
  !    z3b).*(cosB+a./rb)+(y3b-a)./rb.*((1-2*nu)*cotB-2*nu*y1./(rb+y3b)-a*...
  !    y1./rb.^2)-(y3b-a)./(rb+z3b).*(cosB*sinB+(rb*cosB+y3b)*cotB./rb.*...
  !    (2*(1-nu)*cosB-(rb*cosB+y3b)./(rb+z3b))+a./rb.*(sinB-y3b.*z1b./...
  !    rb.^2-z1b.*(rb*cosB+y3b)./(rb.*(rb+z3b)))));
  v3cb2 = b2/4.0D0/pi/(1.0D0-nu)*(-2.0D0*(1.0D0-nu)*(1.0D0-2.0D0*nu)*cotB*(log(rb+y3b)-cosB*  &
       &  log(rb+z3b))-2.0D0*(1.0D0-nu)*y1/(rb+y3b)*(2.0D0*nu+a/rb)+2.0D0*(1.0D0-nu)*z1b/(rb+ &
       &  z3b)*(cosB+a/rb)+(y3b-a)/rb*((1.0D0-2.0D0*nu)*cotB-2.0D0*nu*y1/(rb+y3b)-a*          &
       &  y1/rb**2)-(y3b-a)/(rb+z3b)*(cosB*sinB+(rb*cosB+y3b)*cotB/rb*                        &
       &  (2.0D0*(1.0D0-nu)*cosB-(rb*cosB+y3b)/(rb+z3b))+a/rb*(sinB-y3b*z1b/                  &
       &  rb**2-z1b*(rb*cosB+y3b)/(rb*(rb+z3b)))))

  !v1cb3 = b3/4/pi/(1-nu)*((1-2*nu)*(y2./(rb+y3b).*(1+a./rb)-y2*cosB./(rb+...
  !    z3b).*(cosB+a./rb))-y2.*(y3b-a)./rb.*(a./rb.^2+1./(rb+y3b))+y2.*...
  !    (y3b-a)*cosB./(rb.*(rb+z3b)).*((rb*cosB+y3b)./(rb+z3b).*(cosB+a./...
  !    rb)+a.*y3b./rb.^2));
  v1cb3 = b3/4.0D0/pi/(1.0D0-nu)*((1.0D0-2.0D0*nu)*(y2/(rb+y3b)*(1+a/rb)-y2*cosB/(rb+  &
       & z3b)*(cosB+a/rb))-y2*(y3b-a)/rb*(a/rb**2+1.0D0/(rb+y3b))+y2*                 &
       & (y3b-a)*cosB/(rb*(rb+z3b))*((rb*cosB+y3b)/(rb+z3b)*(cosB+a/                  &
       & rb)+a*y3b/rb**2))

  !v2cb3 = b3/4/pi/(1-nu)*((1-2*nu)*(-sinB*log(rb+z3b)-y1./(rb+y3b).*(1+a./...
  !    rb)+z1b./(rb+z3b).*(cosB+a./rb))+y1.*(y3b-a)./rb.*(a./rb.^2+1./(rb+...
  !    y3b))-(y3b-a)./(rb+z3b).*(sinB*(cosB-a./rb)+z1b./rb.*(1+a.*y3b./...
  !    rb.^2)-1./(rb.*(rb+z3b)).*(y2.^2*cosB*sinB-a*z1b./rb.*(rb*cosB+y3b))));
  v2cb3 = b3/4.0D0/pi/(1.0D0-nu)*((1.0D0-2.0D0*nu)*(-sinB*log(rb+z3b)-y1/(rb+y3b)*(1+a/  &
       & rb)+z1b/(rb+z3b)*(cosB+a/rb))+y1*(y3b-a)/rb*(a/rb**2+1.0D0/(rb+                &
       & y3b))-(y3b-a)/(rb+z3b)*(sinB*(cosB-a/rb)+z1b/rb*(1.0D0+a*y3b/                  &
       & rb**2)-1.0D0/(rb*(rb+z3b))*(y2**2*cosB*sinB-a*z1b/rb*(rb*cosB+y3b))))

  !v3cb3 = b3/4/pi/(1-nu)*(2*(1-nu)*Fib+2*(1-nu)*(y2*sinB./(rb+z3b).*(cosB+...
  !    a./rb))+y2.*(y3b-a)*sinB./(rb.*(rb+z3b)).*(1+(rb*cosB+y3b)./(rb+...
  !    z3b).*(cosB+a./rb)+a.*y3b./rb.^2));
  v3cb3 = b3/4.0D0/pi/(1.0D0-nu)*(2.0D0*(1.0D0-nu)*Fib+2.0D0*(1.0D0-nu)*(y2*sinB/(rb+z3b)*(cosB+  &
       & a/rb))+y2*(y3b-a)*sinB/(rb*(rb+z3b))*(1.0D0+(rb*cosB+y3b)/(rb+                          &
       & z3b)*(cosB+a/rb)+a*y3b/rb**2))

  !v1 = v1cb1+v1cb2+v1cb3;
  !v2 = v2cb1+v2cb2+v2cb3;
  !v3 = v3cb1+v3cb2+v3cb3;
  v1 = v1cb1 + v1cb2 + v1cb3
  v2 = v2cb1 + v2cb2 + v2cb3
  v3 = v3cb1 + v3cb2 + v3cb3

END SUBROUTINE Angdisdispfsc_Bird


subroutine dmake_uvec(x,out)
  implicit none
  real*8, dimension(3), intent(in) :: x
  real*8, dimension(3), intent(out) :: out
  real*8 norm
  norm = norm2(x)
  out(:) = x(:)/norm
end subroutine dmake_uvec

SUBROUTINE dcross_bird (a_vec, b_vec, c_vec)
  ! vector cross product: a x b = c
  IMPLICIT NONE
  REAL*8, DIMENSION(3), INTENT(IN)  :: a_vec, b_vec
  REAL*8, DIMENSION(3), INTENT(OUT) :: c_vec
  c_vec(1) = a_vec(2)*b_vec(3) - a_vec(3)*b_vec(2)
  c_vec(2) = a_vec(3)*b_vec(1) - a_vec(1)*b_vec(3)
  c_vec(3) = a_vec(1)*b_vec(2) - a_vec(2)*b_vec(1)
END SUBROUTINE Dcross_Bird

#undef ANGDISDISP
#undef ANGDISDISP_FSC

