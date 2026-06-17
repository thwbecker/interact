  ! TDdispHS 
! calculates displacements associated with a triangular dislocation in an 
! elastic half-space.
!
! TD: Triangular Dislocation
! EFCS: Earth-Fixed Coordinate System
! TDCS: Triangular Dislocation Coordinate System
! ADCS: Angular Dislocation Coordinate System
! 
! INPUTS
! X, Y and Z: 
! Coordinates of calculation points in EFCS (East, North, Up). X, Y and Z 
! must have the same size. 
!
! P1,P2 and P3:
! Coordinates of TD vertices in EFCS. 
! 
! Ss, Ds and Ts:
! TD slip vector components (Strike-slip, Dip-slip, Tensile-slip).
!
! nu:
! Poisson's ratio.
!
! OUTPUTS
! ue, un and uv:
! Calculated displacement vector components in EFCS. ue, un and uv have
! the same unit as Ss, Ds and Ts in the inputs.
! 
! 
! Example: Calculate and plot the first component of displacement vector 
! on a regular grid.
! 
! [X,Y,Z] = meshgrid(-3:.02:3,-3:.02:3,-5);
! [ue,un,uv] = TDdispHS(X,Y,Z,[-1 0 0],[1 -1 -1],[0 1.5 -2],-1,2,3,.25);
! h = surf(X,Y,reshape(ue,size(X)),'edgecolor','none');
! view(2)
! axis equal
! axis tight
! set(gcf,'renderer','painters')

! Reference journal article: 
! Nikkhoo, M., Walter, T. R. (2015): Triangular dislocation: an analytical,
! artefact-free solution. - Geophysical Journal International, 201, 
! 1117-1139. doi: 10.1093/gji/ggv035

! Copyright (c) 2014 Mehdi Nikkhoo
! 
! Permission is hereby granted, free of charge, to any person obtaining a 
! copy of this software and associated documentation files 
! (the "Software"), to deal in the Software without restriction, including 
! without limitation the rights to use, copy, modify, merge, publish, 
! distribute, sublicense, and/or sell copies of the Software, and to permit
! persons to whom the Software is furnished to do so, subject to the 
! following conditions:
! 
! The above copyright notice and this permission notice shall be included 
! in all copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
! NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
! DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
! OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
! USE OR OTHER DEALINGS IN THE SOFTWARE.

! I appreciate any comments or bug reports.

! Update report No. 1:
! 1) The solid angle calculation was modified to account for the "sign bit"
! of the absolute zero values in the numerator of the "atan2" function. 
! These zeros correspond to the points that are exactly on the TD plane 
! (see Lines 240-245).
! 
! 2) The three vertices of a TD must not be located on the free surface 
! simultaneously. The displacements corresponding to such a case are set to
! zero (see Lines 99-103 here). Accordingly, the two "if" expressions at
! Lines 107-109 and 116-120 in the previous version were removed.

! Mehdi Nikkhoo
! created: 2013.1.24
! Last modified: 2016.10.18
! 
! Section 2.1, Physics of Earthquakes and Volcanoes
! Department 2, Geophysics
! Helmholtz Centre Potsdam
! German Research Centre for Geosciences (GFZ)
! 
! email: 
! mehdi.nikkhoo@gfz-potsdam.de 
! mehdi.nikkhoo@gmail.com
! 
! website:
! http://www.volcanodeformation.com
!
!
! F90 CONVERSION - November 2023
!
! downloaded 11,2023 by TWB - converted to F90, with some
! cross-referencing/inspiration by Peter Bird's half space
! displacement conversion
!
! loc is location vector, X,Y,Z
! u is displacement vector, East, North, Vertical 
subroutine  TDdispHS(loc,P1,P2,P3,Ss,Ds,Ts,nu,u)
  implicit none
  REAL*8, PARAMETER :: pi = 3.14159265358979D0
  real*8,intent(in) :: ss,ds,ts,nu
  real*8,intent(in),dimension(3) :: p1,p2,p3,loc
  real*8,intent(out),dimension(3) :: u
  real*8 :: ueMS,unMS,uvMS,ueFSC,unFSC,uvFSC,ueIS,unIS,uvIS,&
       u1Tn,v1Tn,w1Tn,u2Tn,v2Tn,w2Tn,u3Tn,v3Tn,w3Tn,xp,yp,zp,xn,yn,zn,&
       u1Tp,v1Tp,w1Tp,u3Tp,v3Tp,w3Tp
  real*8,dimension(3,3) :: At
  real*8,dimension(3) :: p1n,p2n,p3n,p_1,p_2,p_3,a,b,c
  logical :: caseplog,caselog,casezlog
  integer :: trimode
  if ((loc(3) .gt.0.d0).or.(P1(3)>0.d0).or.(P2(3)>0.d0).or.(P3(3)>0))then
     print *,'Half-space solution: Z coordinates must be negative!'
     stop
  else if((P1(3)==0.d0).and.(P2(3)==0.d0).and.(P3(3)==0.d0))then
     u = 0.d0
     return
  end  if


  ! Calculate main dislocation contribution to displacements
  call TDdispFS(loc(1),loc(2),loc(3),P1,P2,P3,Ss,Ds,Ts,nu,ueMS,unMS,uvMS);
  !print *,ueMS,unMS,uvMS
  !stop
  ! Calculate harmonic function contribution to displacements
  call TDdisp_HarFunc(loc(1),loc(2),loc(3),P1,P2,P3,Ss,Ds,Ts,nu,ueFSC,unFSC,uvFSC);
  !print *,ueFSC,unFSC,uvFSC
  ! Calculate image dislocation contribution to displacements
  p1n(1:2) = p1(1:2);p2n(1:2)=p2(1:2);p3n(1:2) =  p3(1:2);
  p1n(3)  = -p1(3);  p2n(3)  = -p2(3);p3n(3)   = -p3(3);
  call TDdispFS(loc(1),loc(2),loc(3),P1n,P2n,P3n,Ss,Ds,Ts,nu,ueIS,unIS,uvIS);
  !print *,ueIS,unIS,uvIS
  ! Calculate the complete displacement vector components in EFCS
  u(1) = ueMS+ueIS+ueFSC;
  u(2) = unMS+unIS+unFSC;
  u(3) = uvMS+uvIS+uvFSC;

end subroutine  TDdispHS

subroutine TDdispFS(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,nu,ue,un,uv)
  USE, INTRINSIC :: IEEE_ARITHMETIC
  implicit none
  real*8,intent(in) :: x,y,z,ss,ds,ts,nu
  real*8,intent(in),dimension(3) :: p1,p2,p3
  real*8,intent(out) :: ue,un,uv
  REAL*8, PARAMETER :: pi = 3.14159265358979D0
  
  real*8,dimension(3) :: p1_,p2_,p3_,e12,e13,e23,a,b,c
  real*8 :: bx,by,bz,x_,y_,z_,fin,fid,fi,aA,aB,aC,na,nb,nc,u,v,w,&
       u1,u2,u3,v1,v2,v3,w1,w2,w3
  real*8,dimension(3,3) :: Ar
  logical :: casepLog,casenlog,casezlog
  integer trimode
  ! TDdispFS 
  ! Calculates displacements associated with a triangular dislocation in an
  ! elastic full-space.

  
  call setup_geometry(x,y,z,Ts,Ss,Ds,p1,p2,p3,bx,by,bz,x_,y_,z_,p1_,p2_,p3_,e12,e13,e23,aA,aB,aC,Ar,trimode)
  casepLog  = (Trimode==1);
  casenLog =  (Trimode==-1);
  casezLog =  (Trimode==0);
 
  !print *,trimode,caseplog,casenlog,casezlog
 

  ! Configuration I
  if(casepLog)then
     ! Calculate first angular dislocation contribution
     call TDSetupD(x_,y_,z_,aA,bx,by,bz, nu,p1_,-e13,u1,v1,w1);
     ! Calculate second angular dislocation contribution
     call TDSetupD(x_,y_,z_,aB,bx,by,bz, nu,p2_, e12,u2,v2,w2);
     ! Calculate third angular dislocation contribution
     call TDSetupD(x_,y_,z_,aC,bx,by,bz, nu,p3_, e23,u3,v3,w3)
  end if

  ! Configuration II
 if(casenLog)then 
    ! Calculate first angular dislocation contribution
    call TDSetupD(x_,y_,z_,aA,bx,by,bz,nu,p1_, e13,u1,v1,w1);
    ! Calculate second angular dislocation contribution
    call TDSetupD(x_,y_,z_,aB,bx,by,bz,nu,p2_,-e12,u2,v2,w2);
    ! Calculate third angular dislocation contribution
    call TDSetupD(x_,y_,z_,aC,bx,by,bz,nu,p3_,-e23,u3,v3,w3);
   
 end  if
 if(caseplog.or.casenlog)then
    u = u1+u2+u3;
    v = v1+v2+v3;
    w = w1+w2+w3;
 else if (casezLog)then
    u = IEEE_VALUE(u, IEEE_QUIET_NAN);
    v = IEEE_VALUE(v, IEEE_QUIET_NAN);
    w = IEEE_VALUE(w, IEEE_QUIET_NAN);
 else
    print *,'logic error'
    stop
 end if
    

 ! Calculate the Burgers' function contribution corresponding to the TD
 a = (/-x_,  p1_(2) - y_,  p1_(3) - z_ /);
 b = (/-x_, -y_         , -z_ /);
 c = (/-x_,  p3_(2) - y_,  p3_(3) - z_ /);
 na = norm2(a)
 nb = norm2(b)
 nc = norm2(c)
 
 FiN = (a(1)*(b(2)*c(3)-b(3)*c(2)) - &
      a(2)*(b(1)*c(3)-b(3)*c(1)) + &
      a(3)*(b(1)*c(2)-b(2)*c(1)) ); 
 FiD = (na*nb*nc+dot_product(a,b)*nc+dot_product(a,c)*nb+dot_product(b,c)*na);
 if(FiN == 0.d0)then
    fin = -0.d0; ! Fix the "sign bit" of FiN for x = 0
 end if
 !print *,'fi',fin,fid
 !print *,'b',bx,by,bz
 !print *,'u',u,v,w
 Fi = -2.d0*atan2(FiN,FiD)/4.d0/pi;

 ! Calculate the complete displacement vector components in TDCS
 u = bx*Fi+u;
 v = by*Fi+v;
 w = bz*Fi+w;
 !print *,'uvw',u,v,w
 ! Transform the complete displacement vector components from TDCS into EFCS
 !call matrix_from_3vec(vnorm,vstrike,vdip,Ar)
 call CoordTrans(u,v,w,Ar,ue,un,uv);
end  subroutine TDdispFS

subroutine TDdisp_HarFunc(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,nu,ue,un,uv)
  ! TDdisp_HarFunc calculates the harmonic function contribution to the
  ! displacements associated with a triangular dislocation in a half-space.
  ! The function cancels the surface normal tractions induced by the main and
  ! image dislocations.
  USE, INTRINSIC :: IEEE_ARITHMETIC
  implicit none
  real*8,intent(in) :: x,y,z,ss,ds,ts,nu
  real*8,intent(in),dimension(3) :: p1,p2,p3
  real*8,intent(out) :: ue,un,uv
  real*8,dimension(3,3) :: Ar
  real*8,dimension(3) :: p_1,p_2,p_3,vnorm,vdip,vstrike
       
  real*8 bx,by,bz,b_x,b_y,b_z,u3,v3,w3,u1,v1,w1,u2,v2,w2
  
  bx = Ts; ! Tensile-slip
  by = Ss; ! Strike-slip
  bz = Ds; ! Dip-slip
  
  call get_base_vectors(p1,p2,p3,vstrike,vdip,vnorm)

  ! Transform slip vector components from TDCS into EFCS
  call matrix_from_3vec(vnorm,vstrike,vdip,Ar)

  call CoordTrans(bx,by,bz,Ar,b_x,b_y,b_z);

! Calculate contribution of angular dislocation pair on each TD side 
  call AngSetupFSC(X,Y,Z,b_X,b_Y,b_Z,P1,P2,nu,u1,v1,w1); ! Side P1P2
  call AngSetupFSC(X,Y,Z,b_X,b_Y,b_Z,P2,P3,nu,u2,v2,w2); ! Side P2P3
  call AngSetupFSC(X,Y,Z,b_X,b_Y,b_Z,P3,P1,nu,u3,v3,w3); ! Side P3P1

  ! Calculate total harmonic function contribution to displacements
  ue = u1+u2+u3;
  un = v1+v2+v3;
  uv = w1+w2+w3;
end subroutine TDdisp_HarFunc

subroutine CoordTrans(x1,x2,x3,A,X_1,X_2,X_3)
  implicit none
  real*8,dimension(3,3) :: A
  real*8,intent(in) :: x1,x2,x3
  real*8,intent(out) :: x_1,x_2,x_3
  ! CoordTrans transforms the coordinates of the vectors, from
  ! x1x2x3 coordinate system to X1X2X3 coordinate system. "A" is the
  ! transformation matrix, whose columns e1,e2 and e3 are the unit base 
  ! vectors of the x1x2x3. The coordinates of e1,e2 and e3 in A must be given 
  ! in X1X2X3. The transpose of A (i.e., A') will transform the coordinates 
  ! from X1X2X3 into x1x2x3.
  X_1 = A(1, 1) * x1 + A(1, 2) * x2 + A(1, 3) * x3
  X_2 = A(2, 1) * x1 + A(2, 2) * x2 + A(2, 3) * x3
  X_3 = A(3, 1) * x1 + A(3, 2) * x2 + A(3, 3) * x3

end  subroutine CoordTrans

subroutine trimodefinder(x,y,z,p1,p2,p3,trimode)
  integer,intent(out) :: trimode
  real*8,intent(in) :: x,y,z
  real*8,intent(in),dimension(2) :: p1,p2,p3
  real*8 :: a,b,c
  ! trimodefinder calculates the normalized barycentric coordinates of 
  ! the points with respect to the TD vertices and specifies the appropriate
  ! artefact-free configuration of the angular dislocations for the 
  ! calculations. The input matrices x, y and z share the same size and
  ! correspond to the y, z and x coordinates in the TDCS, respectively. p1,
  ! p2 and p3 are two-component matrices representing the y and z coordinates
  ! of the TD vertices in the TDCS, respectively.
  ! The components of the output (trimode) corresponding to each calculation 
  ! points, are 1 for the first configuration, -1 for the second 
  ! configuration and 0 for the calculation point that lie on the TD sides.
  

  a = ((p2(2)-p3(2))*(x-p3(1))+(p3(1)-p2(1))*(y-p3(2)))/&
       ((p2(2)-p3(2))*(p1(1)-p3(1))+(p3(1)-p2(1))*(p1(2)-p3(2)));
  b = ((p3(2)-p1(2))*(x-p3(1))+(p1(1)-p3(1))*(y-p3(2)))/&
       ((p2(2)-p3(2))*(p1(1)-p3(1))+(p3(1)-p2(1))*(p1(2)-p3(2)));
  c = 1.d0-a-b;

  trimode = 1
  IF ((a <= 0.0D0) .AND. (b > c) .AND. (c > a)) trimode = -1
  IF ((b <= 0.0D0) .AND. (c > a) .AND. (a > b)) trimode = -1
  IF ((c <= 0.0D0) .AND. (a > b) .AND. (b > c)) trimode = -1
  IF ((a == 0.0D0) .AND. (b >= 0.0D0) .AND. (c >= 0.0D0)) trimode = 0
  IF ((a >= 0.0D0) .AND. (b == 0.0D0) .AND. (c >= 0.0D0)) trimode = 0
  IF ((a >= 0.0D0) .AND. (b >= 0.0D0) .AND. (c == 0.0D0)) trimode = 0
  IF ((trimode == 0) .AND. (z /= 0.0D0)) trimode = 1

end  subroutine trimodefinder

subroutine TDSetupD(x,y,z,alpha,bx,by,bz,nu,TriVertex,SideVec,u,v,w)
  ! TDSetupD transforms coordinates of the calculation points as well as 
  ! slip vector components from ADCS into TDCS. It then calculates the 
  ! displacements in ADCS and transforms them into TDCS.
  implicit none
  real*8,intent(in) :: x,y,z,alpha,bx,by,bz,nu
  real*8,intent(out) :: u,v,w
  real*8,intent(in),dimension(3) :: sidevec,trivertex
  REAL*8 :: by1, bz1, v0, w0,y1,z1
  REAL*8, DIMENSION(2, 2) :: A
  REAL*8, PARAMETER :: pi = 3.14159265358979D0

  call get_tdcd_adcs(sidevec,trivertex,y,z,by,bz,A,y1,z1,by1,bz1)

  !print *,x,y1,z1,-pi+alpha,bx,by1,bz1,nu
  ! Calculate displacements associated with an angular dislocation in ADCS
  call AngDisDisp(x,y1,z1,-pi+alpha,bx,by1,bz1,nu,u,v0,w0);

  ! Transform displacements from ADCS into TDCS
  !r3 = A'*[v0';w0'];
  !v = r3(1,:)';
  !w = r3(2,:)';
  v = A(1, 1) * v0 + A(2, 1) * w0
  w = A(1, 2) * v0 + A(2, 2) * w0

  !print *,'uvw',u,v,w
end  subroutine TDSetupD

subroutine AngSetupFSC(X,Y,Z,bX,bY,bZ,PA,PB,nu,ue,un,uv)
  ! AngSetupFSC calculates the Free Surface Correction to displacements 
  ! associated with angular dislocation pair on each TD side.
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: bX, bY, bZ, nu, X, Y, Z
  REAL*8, DIMENSION(3), INTENT(IN) :: PA, PB
  REAL*8, INTENT(OUT) :: ue, un, uv
  REAL*8, PARAMETER :: pi = 3.14159265358979D0
  REAL*8 :: b1, b2, b3, beta, v1, v1A, v1B, v2, v2A, v2B, &
       v3, v3A, v3B, y1A, y1AB, y1B, y2A, y2AB, y2B, y3A, y3AB, y3B
  REAL*8, DIMENSION(3) :: ey1, ey2, ey3, eZ, SideVec, SideUvec
  REAL*8, DIMENSION(3, 3) :: A, At
  REAL*8, PARAMETER :: eps = 1.0D-15 ! a crude approximation of the MatLab "eps" constant.
  LOGICAL :: I
  
  ! Calculate TD side vector and the angle of the angular dislocation pair
  SideVec = PB-PA;
  eZ = (/0.d0, 0.d0, 1.0d0/);
  beta = acos(dot_product(-SideVec,eZ)/norm2(SideVec));

  if ((abs(beta)<eps).or.(abs(pi-beta)<eps))then 
     ue = 0.0d0
     un = 0.d0
     uv = 0.d0
  else
     ey1 = (/ SideVec(1),Sidevec(2),0.d0/);
     call normalize_vec(ey1,ey1)
     ey3 = -eZ;
     call dcross(ey3,ey1,ey2);
     call matrix_from_3vec(ey1,ey2,ey3,A)

  
    ! Transform coordinates from EFCS to the first ADCS
    call CoordTrans(X-PA(1),Y-PA(2),Z-PA(3),A,y1A,y2A,y3A);
    ! Transform coordinates from EFCS to the second ADCS
    call CoordTrans(SideVec(1),SideVec(2),SideVec(3),A,y1AB,y2AB,y3AB);
    y1B = y1A-y1AB;
    y2B = y2A-y2AB;
    y3B = y3A-y3AB;
    
    ! Transform slip vector components from EFCS to ADCS
    call CoordTrans(bX,bY,bZ,A,b1,b2,b3);
    
    ! Determine the best artefact-free configuration for the calculation
    ! points near the free surface
    I = ((beta*y1A)>=0);
    if(I)then 
       ! Configuration I
       call AngDisDispFSC(y1A,y2A,y3A,&
            -pi+beta,b1,b2,b3,nu,-PA(3),v1A,v2A,v3A);
       call AngDisDispFSC(y1B,y2B,y3B,&
            -pi+beta,b1,b2,b3,nu,-PB(3),v1B,v2B,v3B);
    else
       ! Configuration II
       call AngDisDispFSC(y1A,y2A,y3A,&
            beta,b1,b2,b3,nu,-PA(3),v1A,v2A,v3A);
       call AngDisDispFSC(y1B,y2B,y3B,&
            beta,b1,b2,b3,nu,-PB(3),v1B,v2B,v3B);
    end if
    
    ! Calculate total Free Surface Correction to displacements in ADCS
    v1 = v1B-v1A;
    v2 = v2B-v2A;
    v3 = v3B-v3A;
    At = TRANSPOSE(A)
    ! Transform total Free Surface Correction to displacements from ADCS 
    ! to EFCS
    call CoordTrans(v1,v2,v3,At,ue,un,uv);
 end if
end subroutine AngSetupFSC

subroutine AngDisDisp(x,y,z,alpha,bx,by,bz,nu,u,v,w)
  ! AngDisDisp calculates the "incomplete" displacements (without the 
  ! Burgers' function contribution) associated with an angular dislocation in
  ! an elastic full-space.
  implicit none
  real*8,intent(in) :: x,y,z,alpha,bx,by,bz,nu
  real*8,intent(out) :: u,v,w
  real*8 :: cosA,sinA,eta,zeta,r,zz,ux,uy,uz,vx,vy,vz,wx,wy,wz;

  REAL*8, PARAMETER :: pi = 3.14159265358979D0
  
  cosA = cos(alpha);
  sinA = sin(alpha);
  eta =  y*cosA - z*sinA;
  zeta = y*sinA + z*cosA;
  r = sqrt(x**2 + y**2 + z**2);
  
  ! Avoid complex results for the logarithmic terms
  if(zeta>r) zeta =  r
  if(z>r)then
     zz = r
  else
     zz = z
  endif
  ! NEED TO OPTIMIZE
  ux = bx/8.0d0/pi/(1.0d0-nu)*(x*y/r/(r-zz)-x*eta/r/(r-zeta));
  vx = bx/8.0d0/pi/(1.0d0-nu)*(eta*sinA/(r-zeta)-y*eta/r/(r-zeta)+&
       y**2/r/(r-zz)+(1.0d0-2.d0*nu)*(cosA*log(r-zeta)-log(r-zz)));
  wx = bx/8.0d0/pi/(1.0d0-nu)*(eta*cosA/(r-zeta)-y/r-eta*zz/r/(r-zeta)-&
       (1.0d0-2.0d0*nu)*sinA*log(r-zeta));
  
  uy = by/8.0d0/pi/(1.0d0-nu)*(x**2*cosA/r/(r-zeta)-x**2/r/(r-zz)-&
       (1.0d0-2.0d0*nu)*(cosA*log(r-zeta)-log(r-zz)));
  vy = by*x/8.0d0/pi/(1.0d0-nu)*(y*cosA/r/(r-zeta)-&
       sinA*cosA/(r-zeta)-y/r/(r-zz));
  wy = by*x/8.0d0/pi/(1.0d0-nu)*(zz*cosA/r/(r-zeta)-cosA**2/(r-zeta)+1.0d0/r);
  
  uz = bz*sinA/8.0d0/pi/(1.0d0-nu)*((1.0d0-2.0d0*nu)*log(r-zeta)-x**2/r/(r-zeta));
  vz = bz*x*sinA/8.0d0/pi/(1.0d0-nu)*(sinA/(r-zeta)-y/r/(r-zeta));
  wz = bz*x*sinA/8.0d0/pi/(1.0d0-nu)*(cosA/(r-zeta)-zz/r/(r-zeta));
  
  u = ux+uy+uz;
  v = vx+vy+vz;
  w = wx+wy+wz;
end  subroutine AngDisDisp

subroutine AngDisDispFSC(y1,y2,y3,beta,b1,b2,b3,nu,a,v1,v2,v3)
  ! AngDisDispFSC calculates the harmonic function contribution to the 
  ! displacements associated with an angular dislocation in an elastic 
  ! half-space.
  implicit none 
  real*8,intent(in) :: y1,y2,y3,beta,b1,b2,b3,nu,a
  real*8,intent(out) ::v1,v2,v3
  REAL*8, PARAMETER :: pi = 3.14159265358979D0
  REAL*8 :: cosB, cotB, Fib, r2b, rb, sinB, &
       & v1cb1, v1cb2, v1cb3, v2cb1, v2cb2, v2cb3, v3cb1, v3cb2, v3cb3, &
       & y3b, z1b, z3b

 
  sinB = sin(beta);
  cosB = cos(beta);
  cotB = 1.00D0 / TAN(beta);
  y3b = y3+2*a;
  z1b = y1*cosB+y3b*sinB;
  z3b = -y1*sinB+y3b*cosB;
  r2b = y1**2+y2**2+y3b**2;
  rb = sqrt(r2b);
  ! replace with atan2?
  Fib = 2.d0*atan(-y2/(-(rb+y3b)/tan(beta/2.d0)+y1)); ! The Burgers' function

  v1cb1 = b1/4.0d0/pi/(1.0d0-nu)*(-2.0d0*(1.0d0-nu)*(1.0d0-2.0d0*nu)*Fib*cotB**2+(1.0d0-2.0d0*nu)*y2/&
       (rb+y3b)*((1.0d0-2.0d0*nu-a/rb)*cotB-y1/(rb+y3b)*(nu+a/rb))+(1.0d0-2.0d0*nu)*&
       y2*cosB*cotB/(rb+z3b)*(cosB+a/rb)+a*y2*(y3b-a)*cotB/rb**3+y2*&
       (y3b-a)/(rb*(rb+y3b))*(-(1.0d0-2.0d0*nu)*cotB+y1/(rb+y3b)*(2.0d0*nu+a/rb)+&
       a*y1/rb**2)+y2*(y3b-a)/(rb*(rb+z3b))*(cosB/(rb+z3b)*((rb*&
       cosB+y3b)*((1.0d0-2.0d0*nu)*cosB-a/rb)*cotB+2.0d0*(1.0d0-nu)*(rb*sinB-y1)*cosB)-&
       a*y3b*cosB*cotB/rb**2));
  
  v2cb1 = b1/4.0d0/pi/(1.0d0-nu)*((1.0d0-2.0d0*nu)*((2.0d0*(1.0d0-nu)*cotB**2-nu)*log(rb+y3b)-(2.0d0*&
       (1.0d0-nu)*cotB**2+1.0d0-2.0d0*nu)*cosB*log(rb+z3b))-(1.0d0-2.0d0*nu)/(rb+y3b)*(y1*&
       cotB*(1.0d0-2.0d0*nu-a/rb)+nu*y3b-a+y2**2/(rb+y3b)*(nu+a/rb))-(1.0d0-2.0d0*&
       nu)*z1b*cotB/(rb+z3b)*(cosB+a/rb)-a*y1*(y3b-a)*cotB/rb**3+&
       (y3b-a)/(rb+y3b)*(-2*nu+1.0d0/rb*((1.0d0-2.0d0*nu)*y1*cotB-a)+y2**2/(rb*&
       (rb+y3b))*(2.0d0*nu+a/rb)+a*y2**2/rb**3)+(y3b-a)/(rb+z3b)*(cosB**2-&
       1.0d0/rb*((1.0d0-2.0d0*nu)*z1b*cotB+a*cosB)+a*y3b*z1b*cotB/rb**3-1/(rb*&
       (rb+z3b))*(y2**2*cosB**2-a*z1b*cotB/rb*(rb*cosB+y3b))));
  
  v3cb1 = b1/4.0d0/pi/(1.0d0-nu)*(2.0d0*(1.0d0-nu)*(((1.0d0-2.0d0*nu)*Fib*cotB)+(y2/(rb+y3b)*(2.0d0*&
       nu+a/rb))-(y2*cosB/(rb+z3b)*(cosB+a/rb)))+y2*(y3b-a)/rb*(2.0d0*&
       nu/(rb+y3b)+a/rb**2)+y2*(y3b-a)*cosB/(rb*(rb+z3b))*(1.0d0-2.0d0*nu-&
       (rb*cosB+y3b)/(rb+z3b)*(cosB+a/rb)-a*y3b/rb**2));
  
  v1cb2 = b2/4.0d0/pi/(1.0d0-nu)*((1.0d0-2.0d0*nu)*((2.0d0*(1.0d0-nu)*cotB**2+nu)*log(rb+y3b)-(2.0d0*&
       (1.0d0-nu)*cotB**2+1)*cosB*log(rb+z3b))+(1.0d0-2.0d0*nu)/(rb+y3b)*(-(1.0d0-2.0d0*nu)*&
       y1*cotB+nu*y3b-a+a*y1*cotB/rb+y1**2/(rb+y3b)*(nu+a/rb))-(1-2*&
       nu)*cotB/(rb+z3b)*(z1b*cosB-a*(rb*sinB-y1)/(rb*cosB))-a*y1*&
       (y3b-a)*cotB/rb**3+(y3b-a)/(rb+y3b)*(2.0d0*nu+1.0d0/rb*((1.0d0-2.0d0*nu)*y1*&
       cotB+a)-y1**2/(rb*(rb+y3b))*(2.0d0*nu+a/rb)-a*y1**2/rb**3)+(y3b-a)*&
       cotB/(rb+z3b)*(-cosB*sinB+a*y1*y3b/(rb**3*cosB)+(rb*sinB-y1)/&
       rb*(2.0d0*(1.0d0-nu)*cosB-(rb*cosB+y3b)/(rb+z3b)*(1.d0+a/(rb*cosB)))));
  
  v2cb2 = b2/4.0d0/pi/(1.0d0-nu)*(2.0d0*(1.0d0-nu)*(1.0d0-2.0d0*nu)*Fib*cotB**2+(1.0d0-2.0d0*nu)*y2/&
       (rb+y3b)*(-(1.0d0-2.0d0*nu-a/rb)*cotB+y1/(rb+y3b)*(nu+a/rb))-(1.0d0-2.0d0*nu)*&
       y2*cotB/(rb+z3b)*(1.d0+a/(rb*cosB))-a*y2*(y3b-a)*cotB/rb**3+y2*&
       (y3b-a)/(rb*(rb+y3b))*((1.0d0-2.0d0*nu)*cotB-2*nu*y1/(rb+y3b)-a*y1/rb*&
       (1.0d0/rb+1.0d0/(rb+y3b)))+y2*(y3b-a)*cotB/(rb*(rb+z3b))*(-2.0d0*(1.0d0-nu)*&
       cosB+(rb*cosB+y3b)/(rb+z3b)*(1.d0+a/(rb*cosB))+a*y3b/(rb**2*cosB)));
  
  v3cb2 = b2/4.0d0/pi/(1.0d0-nu)*(-2*(1.0d0-nu)*(1.0d0-2.0d0*nu)*cotB*(log(rb+y3b)-cosB*&
       log(rb+z3b))-2.0d0*(1.0d0-nu)*y1/(rb+y3b)*(2.0d0*nu+a/rb)+2.0d0*(1.0d0-nu)*z1b/(rb+&
       z3b)*(cosB+a/rb)+(y3b-a)/rb*((1.0d0-2.0d0*nu)*cotB-2*nu*y1/(rb+y3b)-a*&
       y1/rb**2)-(y3b-a)/(rb+z3b)*(cosB*sinB+(rb*cosB+y3b)*cotB/rb*&
       (2.0d0*(1.0d0-nu)*cosB-(rb*cosB+y3b)/(rb+z3b))+a/rb*(sinB-y3b*z1b/&
       rb**2-z1b*(rb*cosB+y3b)/(rb*(rb+z3b)))));
  
  v1cb3 = b3/4.0d0/pi/(1.0d0-nu)*((1.0d0-2.0d0*nu)*(y2/(rb+y3b)*(1.d0+a/rb)-y2*cosB/(rb+&
       z3b)*(cosB+a/rb))-y2*(y3b-a)/rb*(a/rb**2+1.0d0/(rb+y3b))+y2*&
       (y3b-a)*cosB/(rb*(rb+z3b))*((rb*cosB+y3b)/(rb+z3b)*(cosB+a/&
       rb)+a*y3b/rb**2));
  
  v2cb3 = b3/4.0d0/pi/(1.0d0-nu)*((1.0d0-2.0d0*nu)*(-sinB*log(rb+z3b)-y1/(rb+y3b)*(1.d0+a/&
       rb)+z1b/(rb+z3b)*(cosB+a/rb))+y1*(y3b-a)/rb*(a/rb**2+1.0d0/(rb+&
       y3b))-(y3b-a)/(rb+z3b)*(sinB*(cosB-a/rb)+z1b/rb*(1.d0+a*y3b/&
       rb**2)-1/(rb*(rb+z3b))*(y2**2*cosB*sinB-a*z1b/rb*(rb*cosB+y3b))));
  
  v3cb3 = b3/4.0d0/pi/(1.0d0-nu)*(2.0d0*(1.0d0-nu)*Fib+2.0d0*(1.0d0-nu)*(y2*sinB/(rb+z3b)*(cosB+&
       a/rb))+y2*(y3b-a)*sinB/(rb*(rb+z3b))*(1.d0+(rb*cosB+y3b)/(rb+&
       z3b)*(cosB+a/rb)+a*y3b/rb**2));
  
  v1 = v1cb1+v1cb2+v1cb3;
  v2 = v2cb1+v2cb2+v2cb3;
  v3 = v3cb1+v3cb2+v3cb3;
end  subroutine AngDisDispFSC
subroutine matrix_from_3vec(x,y,z,A)
  implicit none
  real*8,intent(in),dimension(3):: x,y,z
  real*8,intent(out),dimension(3,3) :: A
  a(:,1) = x(:)
  a(:,2) = y(:)
  a(:,3) = z(:)
end  subroutine matrix_from_3vec
subroutine dcross(x,y,out)
  implicit none
  real*8, dimension(3), intent(in) :: x, y
  real*8, dimension(3), intent(out) :: out
  
  out(1) = x(2)*y(3) - x(3)*y(2)
  out(2) = x(3)*y(1) - x(1)*y(3)
  out(3) = x(1)*y(2) - x(2)*y(1)
end subroutine dcross

subroutine normalize_vec(x,y)
  implicit none
  real*8,intent(in),dimension(3) :: x
  real*8,intent(out),dimension(3) :: y
  y = x/ norm2(x)
end subroutine normalize_vec

subroutine setup_geometry(x,y,z,Ts,Ss,Ds,p1,p2,p3,bx,by,bz,&
     x_,y_,z_,p1_,p2_,p3_,e12,e13,e23,aA,aB,aC,Ar,trimode)
  implicit none
  real*8,intent(in) :: x,y,z,Ts,Ss,Ds
  real*8,intent(in),dimension(3) :: p1,p2,p3
  integer,intent(out) :: trimode
  real*8,intent(out) :: bx,by,bz,x_,y_,z_,aA,aB,aC
  real*8,intent(out),dimension(3,3) :: Ar
  real*8,intent(out),dimension(3) :: p1_,p2_,p3_,e12,e13,e23
  real*8,dimension(3,3) :: At
  real*8,dimension(3) :: vnorm,vstrike,vdip
  ! burgers vector
  bx = Ts; ! Tensile-slip
  by = Ss; ! Strike-slip
  bz = Ds; ! Dip-slip
  !print *,bx,by,bz
  call get_base_vectors(p1,p2,p3,vstrike,vdip,vnorm)
  
  ! Transform coordinates and slip vector components from EFCS into TDCS

  p1_ = 0.d0
  p2_ = 0.d0
  p3_ = 0.d0

  call matrix_from_3vec(vnorm,vstrike,vdip,Ar)
  At = transpose(Ar)
  
  call CoordTrans(X-P2(1),Y-P2(2),Z-P2(3),At,x_,y_,z_);
  call CoordTrans(P1(1)-P2(1),P1(2)-P2(2),P1(3)-P2(3),At,p1_(1),p1_(2),p1_(3));
  call CoordTrans(P3(1)-P2(1),P3(2)-P2(2),P3(3)-P2(3),At,p3_(1),p3_(2),p3_(3));
  ! Calculate the unit vectors along TD sides in TDCS

  call normalize_vec(p2_ - p1_,e12)
  call normalize_vec(p3_ - p1_,e13)
  call normalize_vec(p3_ - p2_,e23)

  ! Calculate the TD angles
  aA = acos(dot_product( e12,e13));
  aB = acos(dot_product(-e12,e23));
  aC = acos(dot_product( e23,e13));

  ! Determine the best arteact-free configuration for each calculation point
  call trimodefinder(y_,z_,x_,p1_(2:3),p2_(2:3),p3_(2:3),trimode);
  
end subroutine setup_geometry

subroutine get_base_vectors(p1,p2,p3,vstrike,vdip,vnorm)
  implicit none 
  ! Calculate unit strike, dip and normal to TD vectors: For a horizontal TD 
  ! as an exception, if the normal vector points upward, the strike and dip 
  ! vectors point Northward and Westward, whereas if the normal vector points
  ! downward, the strike and dip vectors point Southward and Westward, 
  ! respectively.
  real*8,intent(in),dimension(3) :: p1,p2,p3
  real*8,intent(out),dimension(3) :: vstrike,vdip,vnorm
  real*8,dimension(3) :: ey,ez
  
  call dcross(P2-P1,P3-P1,vnorm);
  call normalize_vec(vnorm,vnorm)
  
  eZ = (/0.d0, 0.0d0, 1.0d0/);
  call dcross(eZ,Vnorm,vstrike);
  if(norm2(Vstrike)==0.d0)then
     eY = (/0.d0, 1.0d0, 0.0d0/);
     vstrike = ey*vnorm(3);
     ! For horizontal elements in case of half-space calculation!!!
     ! Correct the strike vector of image dislocation only
     if(p1(3) .gt. 0.0d0)then 
        vstrike = -vstrike;
     end if
  end if
  call normalize_vec(vstrike,vstrike)
  call dcross(vnorm,vstrike,vdip);
end subroutine get_base_vectors

subroutine get_tdcd_adcs(sidevec,trivertex,y,z,by,bz,A,y1,z1,by1,bz1)
  implicit none
  real*8,intent(in),dimension(3) :: sidevec,trivertex
  real*8,intent(in) :: y,z,by,bz
  real*8,intent(out) :: by1,bz1
  real*8,intent(out),dimension(2,2) :: A
  real*8,intent(out) :: y1,z1

  ! Transformation matrix
  !A = [[SideVec(3);-SideVec(2)] SideVec(2:3)];
  !print *,x,y,z,alpha,bx,by,bz,nu,TriVertex,SideVec
  A(1, 1) =  SideVec(3)
  A(1, 2) = -SideVec(2)
  A(2, 1) =  SideVec(2)
  A(2, 2) =  SideVec(3)

  ! Transform coordinates of the calculation points from TDCS into ADCS
  !r1 = A*[y'-TriVertex(2);z'-TriVertex(3)];
  !y1 = r1(1,:)';
  !z1 = r1(2,:)';
  y1 = A(1, 1) * (y - TriVertex(2)) + A(1, 2) * (z -TriVertex(3))
  z1 = A(2, 1) * (y - TriVertex(2)) + A(2, 2) * (z -TriVertex(3))

  ! Transform the in-plane slip vector components from TDCS into ADCS
  !r2 = A*[by;bz];
  !by1 = r2(1,:)';
  !bz1 = r2(2,:)';
  by1 = A(1, 1) * by + A(1, 2) * bz
  bz1 = A(2, 1) * by + A(2, 2) * bz
end subroutine get_tdcd_adcs
