#define ANG_STRAIN_ROUTINE AngDisStrainFSC
#define ANG_DISS_ROUTINE AngDisStrain

subroutine TDstressHS(loc,P1,P2,P3,Ss,Ds,Ts,mu,lambda,stress,strain)
  ! TDstressHS 
  ! calculates stresses and strains associated with a triangular dislocation 
  ! in an elastic half-space.
  !
  ! TD: Triangular Dislocation
  ! EFCS: Earth-Fixed Coordinate System
  ! TDCS: Triangular Dislocation Coordinate System
  ! ADCS: Angular Dislocation Coordinate System
  ! 
  ! INPUTS 
  ! X, Y and Z:  ---> changed to location vector
  ! Coordinates of calculation points in EFCS (East, North, Up). X, Y and Z 
  ! must have the same size. - THIS NEEDS TO BE 1D here
  !
  ! P1,P2 and P3:
  ! Coordinates of TD vertices in EFCS.
  ! 
  ! Ss, Ds and Ts:
  ! TD slip vector components (Strike-slip, Dip-slip, Tensile-slip).
  !
  ! mu and lambda:
  ! Lame constants.
  !
  ! OUTPUTS
  ! Stress:
  ! Calculated stress tensor components in EFCS. The six columns of Stress 
  ! are Sxx, Syy, Szz, Sxy, Sxz and Syz, respectively. The stress components 
  ! have the same unit as Lame constants.
  !
  ! Strain:
  ! Calculated strain tensor components in EFCS. The six columns of Strain 
  ! are Exx, Eyy, Ezz, Exy, Exz and Eyz, respectively. The strain components 
  ! are dimensionless.
  ! 
  ! 

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
  ! 1) Bug description: The equation for Poisson's ratio (Line 130 and 399 in
  ! the previous version) was valid only for a Poisson solid with 
  ! mu = lambda. 
  ! Bug fixed: The issue has been fixed and the correct equation has been 
  ! replaced (Line 147 and 416 in this version).
  !   
  ! 2) The "Configuration I" (Line 446-457 in the first version) was correct,
  ! but had not been simplified like "Configuration II" in the same
  ! sub-function. The simplified form (Lines 463-468 in this version) has 
  ! been replaced.
  ! 
  ! 3) The three vertices of a TD must not be located on the free surface 
  ! simultaneously. The strains and stresses corresponding to such a case are
  ! set to zero (see Lines 112-115 here). Accordingly, the "if" expression at
  ! Lines 114-119 in the previous version was removed.

  ! Mehdi Nikkhoo
  ! created: 2013.1.28
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
  implicit none 
  real*8, intent(in) :: mu,lambda,Ss,Ds,Ts
  real*8, intent(out),dimension(6) :: stress,strain
  real*8, dimension(3), intent(in) :: p1, p2, p3,loc
  real*8, dimension(3) :: p1n, p2n, p3n
  real*8, dimension(6) :: StsMS,StrMS,StsFSC,StrFSC,StsIS,StrIS

  if ((loc(3) .gt. 0.d0).or. (P1(3) > 0.d0).or.( P2(3) > 0d0) .or.( P3(3)>0d0))then 
     write(*,*)'Half-space solution: Z coordinates must be negative!'
     stop 
  else if ((P1(3)==0.d0).and.(P2(3)==0.d0).and.(P3(3)==0.d0))then
     Stress = 0.d0
     Strain = 0.d0
     return 
  endif

  !print *,Ts,Ss,Ds
  ! Calculate main dislocation contribution to strains and stresses
  call TDstressFS(loc(1),loc(2),loc(3),P1,P2,P3,Ss,Ds,Ts,mu,lambda,StsMS,StrMS)
  !print *,stsms
  !print *,strms
  ! Calculate harmonic function contribution to strains and stresses
  call TDstress_HarFunc(loc(1),loc(2),loc(3),P1,P2,P3,Ss,Ds,Ts,mu,lambda,StsFSC,StrFSC)

  ! Calculate image dislocation contribution to strains and stresses
  p1n(1:2) = p1(1:2);p2n(1:2)=p2(1:2);p3n(1:2) =  p3(1:2);
  p1n(3)  = -p1(3);  p2n(3)  = -p2(3);p3n(3)   = -p3(3);

  call TDstressFS(loc(1),loc(2),loc(3),P1n,P2n,P3n,Ss,Ds,Ts,mu,lambda,StsIS,StrIS)

  ! Calculate the complete stress and strain tensor components in EFCS
  Stress = StsMS+StsIS+StsFSC
  Strain = StrMS+StrIS+StrFSC

end subroutine TDstressHS

subroutine TDstressFS(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda,Stress,Strain)
  USE, INTRINSIC :: IEEE_ARITHMETIC
  ! TDstressFS 
  ! Calculates stresses and strains associated with a triangular dislocation 
  ! in an elastic full-space.
  implicit none
  real*8,intent(in) :: x,y,z,Ss,Ds,Ts,mu,lambda
  real*8,dimension(3) :: p1,p2,p3
  real*8,intent(out),dimension(6) :: stress,strain

  real*8 nu,bx,by,bz,x_,y_,z_,aA,aB,aC,&
       Exx,Eyy,Ezz,Exy,Exz,Eyz,&
       Exxr,Eyyr,Ezzr,Exyr,Exzr,Eyzr,&
       Sxx,Syy,Szz,Sxy,Sxz,Syz,two_mu
  real*8 Exx1,Eyy1,Ezz1,Exy1,Exz1,Eyz1,&
       Exx2,Eyy2,Ezz2,Exy2,Exz2,Eyz2,&
       Exx3,Eyy3,Ezz3,Exy3,Exz3,Eyz3,&
       etracel
  real*8,dimension(3) :: e12,e13,e23,p1_,p2_,p3_
  real*8,dimension(3,3) :: Ar
  INTEGER :: trimode
  logical :: caseplog,casenlog,casezlog

  nu = lambda/((lambda+mu)*2.d0) 
  two_mu = 2.0d0*mu

  ! this is the part shared between deformation and stress compution
  !print *,Ts,Ss,Ds
  call setup_geometry(x,y,z,Ts,Ss,Ds,p1,p2,p3,bx,by,bz,&
       x_,y_,z_,p1_,p2_,p3_,e12,e13,e23,aA,aB,aC,Ar,trimode)
  !print *,bx,by,bz
  !stop
  !print *,aA,aB,aC

  casepLog = (Trimode ==  1)
  casenLog = (Trimode == -1)
  casezLog = (Trimode ==  0)
  !print *,caseplog,casenlog,casezlog
  ! Configuration I
  if (casepLog) then
     !print *,x_,y_,z_,aA,bx,by,bz,nu,p1_,-e13
     ! Calculate first angular dislocation contribution
     call TDSetupS(x_,y_,z_,aA,bx,by,bz,nu,p1_,-e13,Exx1,Eyy1,Ezz1,Exy1,Exz1,Eyz1);
     ! Calculate second angular dislocation contribution
     call TDSetupS(x_,y_,z_,aB,bx,by,bz,nu,p2_, e12,Exx2,Eyy2,Ezz2,Exy2,Exz2,Eyz2);
     ! Calculate third angular dislocation contribution
     call TDSetupS(x_,y_,z_,aC,bx,by,bz,nu,p3_, e23,Exx3,Eyy3,Ezz3,Exy3,Exz3,Eyz3);
     !print *,'exx1'
     !print *,Exx1,Eyy1,Ezz1,Exy1,Exz1,Eyz1
     !print *,'exx2'
     !print *,Exx2,Eyy2,Ezz2,Exy2,Exz2,Eyz2
     !print *,'exx3'
     !print *,Exx3,Eyy3,Ezz3,Exy3,Exz3,Eyz3
     !stop
  else if (casenLog) then 
     ! Configuration II
     ! Calculate first angular dislocation contribution
     call TDSetupS(x_,y_,z_,aA,bx,by,bz,nu,p1_, e13,Exx1,Eyy1,Ezz1,Exy1,Exz1,Eyz1);
     ! Calculate second angular dislocation contribution
     call TDSetupS(x_,y_,z_,aB,bx,by,bz,nu,p2_,-e12,Exx2,Eyy2,Ezz2,Exy2,Exz2,Eyz2);
     ! Calculate third angular dislocation contribution
     call TDSetupS(x_,y_,z_,aC,bx,by,bz,nu,p3_,-e23,Exx3,Eyy3,Ezz3,Exy3,Exz3,Eyz3);    
  end if

  ! Calculate the strain tensor components in TDCS
  if (casepLog.or.casenLog)then 
     exx = Exx1+Exx2+Exx3
     eyy = Eyy1+Eyy2+Eyy3
     ezz = Ezz1+Ezz2+Ezz3
     exy = Exy1+Exy2+Exy3
     exz = Exz1+Exz2+Exz3
     eyz = Eyz1+Eyz2+Eyz3
  else if (casezLog)then 
     exx = IEEE_VALUE(exx, IEEE_QUIET_NAN)
     eyy = IEEE_VALUE(eyy, IEEE_QUIET_NAN)
     ezz = IEEE_VALUE(ezz, IEEE_QUIET_NAN)
     exy = IEEE_VALUE(exy, IEEE_QUIET_NAN)
     exz = IEEE_VALUE(exz, IEEE_QUIET_NAN)
     eyz = IEEE_VALUE(eyz, IEEE_QUIET_NAN)
  else
     print *,'logic error'
     stop
  endif
  !print *,'strain_pre',exx,eyy,ezz,exy,exz,eyz
  ! Transform the strain tensor components from TDCS into EFCS
  call TensTrans(exx,eyy,ezz,exy,exz,eyz,Ar,exxr,eyyr,ezzr,exyr,exzr,eyzr)

  ! Calculate the stress tensor components in EFCS
  etracel = (Exxr+Eyyr+Ezzr)*lambda
  Sxx = two_mu*Exxr+etracel
  Syy = two_mu*Eyyr+etracel
  Szz = two_mu*Ezzr+etracel
  Sxy = two_mu*Exyr;
  Sxz = two_mu*Exzr;
  Syz = two_mu*Eyzr;

  Strain = (/Exxr,Eyyr,Ezzr,Exyr,Exzr,Eyzr/)
  !print *,'strainFS',strain
  !stop
  Stress = (/Sxx,Syy,Szz,Sxy,Sxz,Syz/)
end subroutine TDstressFS

subroutine TDstress_HarFunc(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda,stress,strain)
  implicit none
  real*8,intent(in) :: x,y,z,ss,ds,ts,mu,lambda
  real*8,intent(in),dimension(3) :: p1,p2,p3
  real*8,intent(out),dimension(6) :: stress,strain
  real*8,dimension(6) :: stress1,strain1,stress2,strain2,stress3,strain3
  real*8,dimension(3) :: vstrike,vnorm,vdip
  real*8,dimension(3,3) :: A
  real*8 bx,by,bz,bx_,by_,bz_
  ! TDstress_HarFunc calculates the harmonic function contribution to the
  ! strains and stresses associated with a triangular dislocation in a 
  ! half-space. The function cancels the surface normal tractions induced by 
  ! the main and image dislocations.

  bx = Ts; ! Tensile-slip
  by = Ss; ! Strike-slip
  bz = Ds; ! Dip-slip

  call get_base_vectors(p1,p2,p3,vstrike,vdip,vnorm)

  ! Transform slip vector components from TDCS into EFCS
  call matrix_from_3vec(Vnorm, Vstrike, Vdip,A)
  call CoordTrans(bx,by,bz,A,bx_,by_,bz_) 
  !print *,'b_',bx_,by_,bz_
  ! Calculate contribution of angular dislocation pair on each TD side
  !print *,X,Y,Z,bX_,bY_,bZ_,P1,P2,mu,lambda
  !print *,P2,P3
  call AngSetupFSC_S(x,y,z,bx_,by_,bz_,p1,p2,mu,lambda,stress1,strain1); ! P1P2
  call AngSetupFSC_S(x,y,z,bx_,by_,bz_,p2,P3,mu,lambda,stress2,strain2); ! P2P3
  call AngSetupFSC_S(x,y,z,bx_,by_,bz_,p3,P1,mu,lambda,stress3,strain3); ! P3P1
  !print *,'strainhar1',strain1/1e-2
  !print *,'strainhar1',strain2/1e-2
  !print *,'strainhar1',strain3/1e-2
  !stop
  ! Calculate total harmonic function contribution to strains and stresses
  stress = stress1+stress2+stress3;
  strain = strain1+strain2+strain3;
  !print *,'strainhar',strain
end subroutine TDstress_HarFunc

subroutine TensTrans(Txx1,Tyy1,Tzz1,Txy1,Txz1,Tyz1,Am,Txx2,Tyy2,Tzz2,Txy2,Txz2,Tyz2)
  implicit none
  real*8,intent(in) :: Txx1,Tyy1,Tzz1,Txy1,Txz1,Tyz1
  real*8,intent(in),dimension(3,3) :: Am
  real*8,intent(out) :: Txx2,Tyy2,Tzz2,Txy2,Txz2,Tyz2
  real*8, dimension(9) :: a
  A(1:3) = Am(:,1);A(4:6) = Am(:,2);A(7:9) = Am(:,3)

  ! TensTrans Transforms the coordinates of tensors,from x1y1z1 coordinate
  ! system to x2y2z2 coordinate system. "A" is the transformation matrix, 
  ! whose columns e1,e2 and e3 are the unit base vectors of the x1y1z1. The 
  ! coordinates of e1,e2 and e3 in A must be given in x2y2z2. The transpose 
  ! of A (i.e., A') does the transformation from x2y2z2 into x1y1z1.
  Txx2 = A(1)**2*Txx1+2.0d0*A(1)*A(4)*Txy1+2.0d0*A(1)*A(7)*Txz1+2.0d0*A(4)*A(7)*Tyz1+&
       A(4)**2*Tyy1+A(7)**2*Tzz1;
  Tyy2 = A(2)**2*Txx1+2.0d0*A(2)*A(5)*Txy1+2.0d0*A(2)*A(8)*Txz1+2.0d0*A(5)*A(8)*Tyz1+&
       A(5)**2*Tyy1+A(8)**2*Tzz1;
  Tzz2 = A(3)**2*Txx1+2.0d0*A(3)*A(6)*Txy1+2.0d0*A(3)*A(9)*Txz1+2.0d0*A(6)*A(9)*Tyz1+&
       A(6)**2*Tyy1+A(9)**2*Tzz1;
  Txy2 = A(1)*A(2)*Txx1+(A(1)*A(5)+A(2)*A(4))*Txy1+(A(1)*A(8)+&
       A(2)*A(7))*Txz1+(A(8)*A(4)+A(7)*A(5))*Tyz1+A(5)*A(4)*Tyy1+&
       A(7)*A(8)*Tzz1;
  Txz2 = A(1)*A(3)*Txx1+(A(1)*A(6)+A(3)*A(4))*Txy1+(A(1)*A(9)+&
       A(3)*A(7))*Txz1+(A(9)*A(4)+A(7)*A(6))*Tyz1+A(6)*A(4)*Tyy1+&
       A(7)*A(9)*Tzz1;
  Tyz2 = A(2)*A(3)*Txx1+(A(3)*A(5)+A(2)*A(6))*Txy1+(A(3)*A(8)+&
       A(2)*A(9))*Txz1+(A(8)*A(6)+A(9)*A(5))*Tyz1+A(5)*A(6)*Tyy1+&
       A(8)*A(9)*Tzz1;
end subroutine TensTrans



subroutine TDSetupS(x,y,z,alpha,bx,by,bz,nu,TriVertex,SideVec,exxt,eyyt,ezzt,exyt,exzt,eyzt)
  ! TDSetupS transforms coordinates of the calculation points as well as 
  ! slip vector components from ADCS into TDCS. It then calculates the 
  ! strains in ADCS and transforms them into TDCS.
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: alpha, bx, by, bz, nu, x, y, z
  REAL*8, DIMENSION(3), INTENT(IN) :: SideVec, TriVertex
  REAL*8, INTENT(OUT) :: exxt,eyyt,ezzt,exyt,exzt,eyzt
  REAL*8 :: exx,eyy,ezz,exy,exz,eyz
  REAL*8, PARAMETER :: pi = 3.14159265358979D0
  REAL*8, DIMENSION(2, 2) :: A
  REAL*8, DIMENSION(3, 3) :: B
  real*8 by1,bz1,y1,z1


  call get_tdcd_adcs(sidevec,trivertex,y,z,by,bz,A,y1,z1,by1,bz1)


  ! Calculate strains associated with an angular dislocation in ADCS
  call ANG_DISS_ROUTINE(x,y1,z1,-pi+alpha,bx,by1,bz1,nu,exx,eyy,ezz,exy,exz,eyz);

  ! Transform strains from ADCS into TDCS
  B(1,1) = 1.0d0;B(1,2:3)=0.0d0;
  B(2:3,1) = 0.d0;
  !B(2:3,2:3) = transpose(A)
  B(2,2)=A(1,1);B(2,3)=A(2,1)
  B(3,2)=A(1,2);B(3,3)=A(2,2)
  ! B matrix should look like
  ! 1   0   0 
  ! 0 a11 a21
  ! 0 a12 a22
  !
  !B = [[1 0 0];[zeros(2,1),A']]; ! 3x3 Transformation matrix
  call TensTrans(exx,eyy,ezz,exy,exz,eyz,B,exxt,eyyt,ezzt,exyt,exzt,eyzt)
end subroutine TDSetupS

subroutine AngSetupFSC_S(X,Y,Z,bX,bY,bZ,PA,PB,mu,lambda,Stress,Strain)
  IMPLICIT NONE
  real*8, intent(in) :: X,Y,Z,bX,bY,bZ,mu,lambda
  real*8, dimension(6), intent(out) :: stress,strain
  REAL*8, DIMENSION(3), INTENT(IN) :: PA, PB
  REAL*8, PARAMETER :: pi = 3.14159265358979D0
  REAL*8 :: beta,two_mu,ltrace_E,nu,Exx,Eyy,Ezz,Exy,Exz,Eyz,&
       Sxx,Syy,Szz,Sxy,Sxz,Syz,&
       v11,v22,v33,v12,v23,v13,v11A,v22A,v33A,v12A,v23A,v13A,&
       v11B,v22B,v33B,v12B,v23B,v13B
  real*8, dimension(3) :: SideVec,eZ,ey1,ey2,ey3
  real*8 y1A,y2A,y3A,y1B,y2B,y3B,y1AB,y2AB,y3AB,b1,b2,b3
  real*8, dimension(3,3) :: A,At

  REAL*8, PARAMETER :: eps = 1.0D-15 ! a crude approximation of the MatLab "eps" constant.
  LOGICAL :: I
  ! AngSetupFSC_S calculates the Free Surface Correction to strains and 
  ! stresses associated with angular dislocation pair on each TD side.

  nu = lambda/(lambda+mu)/2.0d0; ! Poisson's ratio
  two_mu = 2.0d0 * mu

  ! Calculate TD side vector and the angle of the angular dislocation pair
  SideVec = PB - PA;
  eZ = (/ 0.d0, 0.d0, 1.d0 /)
  beta = acos(-dot_product(SideVec,eZ)/norm2(SideVec));
  !print *,beta
  if ((abs(beta).lt.eps).or.(abs(pi-beta).lt.eps))then
     Stress = 0.d0
     Strain = 0.d0
  else
     ey1 = (/ SideVec(1), SideVec(2), 0.d0 /)
     call normalize_vec(ey1,ey1)
     ey3 = -eZ;
     call Dcross(ey3,ey1,ey2);
     call matrix_from_3vec(ey1,ey2,ey3,A)
     !print *,A
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
     I = ((beta*y1A) >= 0.d0)
     !print *,I,b1,b2,b3
     ! For singularities at surface
     v11A = 0.d0
     v22A = 0.d0
     v33A = 0.d0
     v12A = 0.d0
     v13A = 0.d0
     v23A = 0.d0

     v11B = 0.d0
     v22B = 0.d0
     v33B = 0.d0
     v12B = 0.d0
     v13B = 0.d0
     v23B = 0.d0
     IF (I) THEN
        ! Configuration I
        call ANG_STRAIN_ROUTINE(y1A,y2A,y3A,-pi+beta,b1,b2,b3,nu,-PA(3),v11A,v22A,v33A,v12A,v13A,v23A);
        call ANG_STRAIN_ROUTINE(y1B,y2B,y3B,-pi+beta,b1,b2,b3,nu,-PB(3),v11B,v22B,v33B,v12B,v13B,v23B);
     else
        ! Configuration II
        call ANG_STRAIN_ROUTINE(y1A,y2A,y3A,beta,b1,b2,b3,nu,    -PA(3),v11A,v22A,v33A,v12A,v13A,v23A);
        call ANG_STRAIN_ROUTINE(y1B,y2B,y3B,beta,b1,b2,b3,nu,    -PB(3),v11B,v22B,v33B,v12B,v13B,v23B);
     endif

     !print *,v11A,v22A,v33A,v12A,v13A,v23A
     !print *,v11B,v22B,v33B,v12B,v13B,v23B
     ! Calculate total Free Surface Correction to strains in ADCS
     v11 = v11B-v11A;
     v22 = v22B-v22A;
     v33 = v33B-v33A;
     v12 = v12B-v12A;
     v13 = v13B-v13A;
     v23 = v23B-v23A;

     At = transpose(A)
     ! Transform total Free Surface Correction to strains from ADCS to EFCS
     call TensTrans(v11,v22,v33,v12,v13,v23,At,Exx,Eyy,Ezz,Exy,Exz,Eyz);

     ! Calculate total Free Surface Correction to stresses in EFCS
     ltrace_E = lambda*(Exx+Eyy+Ezz)

     Sxx = two_mu*Exx+ltrace_E
     Syy = two_mu*Eyy+ltrace_E
     Szz = two_mu*Ezz+ltrace_E
     Sxy = two_mu*Exy;
     Sxz = two_mu*Exz;
     Syz = two_mu*Eyz;

     Strain = (/ Exx,Eyy,Ezz,Exy,Exz,Eyz /)
     Stress = (/ Sxx,Syy,Szz,Sxy,Sxz,Syz /)
  end if

end subroutine AngSetupFSC_S

subroutine AngDisStrain(x,y,z,alpha,bx,by,bz,nu,Exx,Eyy,Ezz,Exy,Exz,Eyz)
  IMPLICIT NONE
  ! AngDisStrain calculates the strains associated with an angular 
  ! dislocation in an elastic full-space.
  real*8,intent(in) :: x,y,z,alpha,bx,by,bz,nu
  real*8,intent(out) :: exx,eyy,ezz,exy,exz,eyz

  real*8 sinA,cosA,eta,zeta,x2,y2,z2,r2,r,r3,rz,r3z,W,W2,Wr,W2r,Wr3,W2r2, C,S
  real*8 rFi_rx,rFi_ry,rFi_rz,r2z2
  REAL*8, PARAMETER :: pi = 3.14159265358979D0

  sinA = sin(alpha);
  cosA = cos(alpha);
  eta = y*cosA-z*sinA;
  zeta = y*sinA+z*cosA;

  x2 = x**2;
  y2 = y**2;
  z2 = z**2;
  r2 = x2+y2+z2;
  r = sqrt(r2);
  r3 = r*r2;
  rz = r*(r-z);
  r2z2 = r2*(r-z)**2;
  r3z = r3*(r-z);

  W = zeta-r;
  W2 = W**2;
  Wr = W*r;
  W2r = W2*r;
  Wr3 = W*r3;
  W2r2 = W2*r2;

  C = (r*cosA-z)/Wr;
  S = (r*sinA-y)/Wr;

  ! Partial derivatives of the Burgers' function
  rFi_rx = (eta/r/(r-zeta)-y/r/(r-z))/4.0d0/pi;
  rFi_ry = (x/r/(r-z)-cosA*x/r/(r-zeta))/4.0d0/pi;
  rFi_rz = (sinA*x/r/(r-zeta))/4.0d0/pi;

  Exx = bx*(rFi_rx)+&
       bx/8.0d0/pi/(1.0d0-nu)*(eta/Wr+eta*x2/W2r2-eta*x2/Wr3+y/rz-&
       x2*y/r2z2-x2*y/r3z)-&
       by*x/8.0d0/pi/(1.0d0-nu)*(((2.0d0*nu+1.0d0)/Wr+x2/W2r2-x2/Wr3)*cosA+&
       (2.0d0*nu+1.0d0)/rz-x2/r2z2-x2/r3z)+&
       bz*x*sinA/8.0d0/pi/(1.0d0-nu)*((2.0d0*nu+1.0d0)/Wr+x2/W2r2-x2/Wr3);

  Eyy = by*(rFi_ry)+&
       bx/8.0d0/pi/(1.0d0-nu)*((1.0d0/Wr+S**2-y2/Wr3)*eta+(2.0d0*nu+1.0d0)*y/rz-y**3/r2z2-&
       y**3/r3z-2.0d0*nu*cosA*S)-&
       by*x/8.0d0/pi/(1.0d0-nu)*(1.0d0/rz-y2/r2z2-y2/r3z+&
       (1.0d0/Wr+S**2-y2/Wr3)*cosA)+&
       bz*x*sinA/8.0d0/pi/(1.0d0-nu)*(1.0d0/Wr+S**2-y2/Wr3);

  Ezz = bz*(rFi_rz)+&
       bx/8.0d0/pi/(1.0d0-nu)*(eta/W/r+eta*C**2-eta*z2/Wr3+y*z/r3+&
       2.0d0*nu*sinA*C)-&
       by*x/8.0d0/pi/(1.0d0-nu)*((1.0d0/Wr+C**2-z2/Wr3)*cosA+z/r3)+&
       bz*x*sinA/8.0d0/pi/(1.0d0-nu)*(1.0d0/Wr+C**2-z2/Wr3);

  Exy = bx*(rFi_ry)/2+by*(rFi_rx)/2-&
       bx/8.0d0/pi/(1.0d0-nu)*(x*y2/r2z2-nu*x/rz+x*y2/r3z-nu*x*cosA/Wr+&
       eta*x*S/Wr+eta*x*y/Wr3)+&
       by/8.0d0/pi/(1.0d0-nu)*(x2*y/r2z2-nu*y/rz+x2*y/r3z+nu*cosA*S+&
       x2*y*cosA/Wr3+x2*cosA*S/Wr)-&
       bz*sinA/8.0d0/pi/(1.0d0-nu)*(nu*S+x2*S/Wr+x2*y/Wr3);

  Exz = bx*(rFi_rz)/2+bz*(rFi_rx)/2-&
       bx/8.0d0/pi/(1.0d0-nu)*(-x*y/r3+nu*x*sinA/Wr+eta*x*C/Wr+&
       eta*x*z/Wr3)+&
       by/8.0d0/pi/(1.0d0-nu)*(-x2/r3+nu/r+nu*cosA*C+x2*z*cosA/Wr3+&
       x2*cosA*C/Wr)-&
       bz*sinA/8.0d0/pi/(1.0d0-nu)*(nu*C+x2*C/Wr+x2*z/Wr3);

  Eyz = by*(rFi_rz)/2+bz*(rFi_ry)/2+&
       bx/8.0d0/pi/(1.0d0-nu)*(y2/r3-nu/r-nu*cosA*C+nu*sinA*S+eta*sinA*cosA/W2-&
       eta*(y*cosA+z*sinA)/W2r+eta*y*z/W2r2-eta*y*z/Wr3)-&
       by*x/8.0d0/pi/(1.0d0-nu)*(y/r3+sinA*cosA**2/W2-cosA*(y*cosA+z*sinA)/&
       W2r+y*z*cosA/W2r2-y*z*cosA/Wr3)-&
       bz*x*sinA/8.0d0/pi/(1.0d0-nu)*(y*z/Wr3-sinA*cosA/W2+(y*cosA+z*sinA)/&
       W2r-y*z/W2r2);
end subroutine AngDisStrain



subroutine AngDisStrainFSC(y1,y2,y3,beta,b1,b2,b3,nu,a,v11,v22,v33,v12,v13,v23)
  ! AngDisStrainFSC calculates the harmonic function contribution to the 
  ! strains associated with an angular dislocation in an elastic half-space.
  IMPLICIT NONE
  real*8, intent(in) :: y1,y2,y3,beta,b1,b2,b3,nu,a
  real*8, intent (out) :: v11,v22,v33,v12,v13,v23
  real*8 y3b,z1b,z3b,rb2,rb,W1,W2,W3,W4,W5,W6,W7,W8,W9,y2s,rb2s,rbc,W6s,W7s
  real*8 rFib_ry2,rFib_ry1,rFib_ry3,cosB,sinB,cotB,N1,N2,N3,N4,two_nu,&
       one_over_rb,one_over_rb2,one_over_rb_p3,one_over_W6,one_over_W6_p2,w1_over_w7,w1_over_w7s,one_over_w7
  REAL*8, PARAMETER :: pi = 3.14159265358979D0, one_quarter = 1.0d0/4.0d0, one_half = 1.0d0/2.0d0

  two_nu = 2.0d0*nu
  
  sinB = sin(beta);
  cosB = cos(beta);
  !cotB = cot(beta);
  cotB = 1.0d0/tan(beta)
  y3b = y3 + 2.0d0*a;
  z1b =  y1*cosB+y3b*sinB;
  z3b = -y1*sinB+y3b*cosB;
  y2s = y2**2

  rb2 =  y1**2+y2s+y3b**2;
  rb2s = rb2**2
  rb = sqrt(rb2);
  one_over_rb = 1.0d0/rb
  one_over_rb2 = 1.0d0/rb2
  rbc = rb**3
  
  one_over_rb_p3 = 1.0d0/rbc
  
  W1 = rb*cosB+y3b;
  W2 = cosB+a/rb;
  W3 = cosB+y3b/rb;
  W4 = nu+a/rb;
  W5 = two_nu+a/rb;
  W6 = rb+y3b;
  W7 = rb+z3b;
  W8 = y3+a;
  W9 = 1.0d0+a/rb/cosB;

  W6s = W6**2
  W7s = W7**2
  
  w1_over_w7 = W1/W7
  w1_over_w7s = W1/W7s
  
  one_over_W6 = 1.0d0/W6
  one_over_W6_p2 = 1.0d0/W6s

  one_over_w7 = 1.0d0/w7
  
  N1 =  1.0d0 - two_nu;
  N2 = -2.0d0 + two_nu;
  N3 =  2.0d0 - two_nu;
  N4 =  1.0d0 - nu;
  
  ! Partial derivatives of the Burgers' function
  rFib_ry2 = z1b/rb/(rb+z3b)-y1/rb/(rb+y3b); ! y2 = x in ADCS
  rFib_ry1 = y2/rb/(rb+y3b)-cosB*y2/rb/(rb+z3b); ! y1 = y in ADCS
  rFib_ry3 = -sinB*y2/rb/(rb+z3b); ! y3 = z in ADCS

  !print *,W1,W2,W3,W4,W5,W6,W7,W8,W9,N1,b1,b2,b3
  v11 = b1*(one_quarter*((N2)*N1*rFib_ry1*cotB**2-N1*y2/W6s*((1.0d0-W5)*cotB-&
       y1/W6*W4)/rb*y1+N1*y2/W6*(a/rbc*y1*cotB-one_over_W6*W4+y1**2/&
       W6s*W4/rb+y1**2/W6*a/rbc)-N1*y2*cosB*cotB/W7s*W2*(y1/&
       rb-sinB)-N1*y2*cosB*cotB/W7*a/rbc*y1-3*a*y2*W8*cotB/rb**5*&
       y1-y2*W8/rbc/W6*(-N1*cotB+y1/W6*W5+a*y1/rb2)*y1-y2*W8/&
       rb2/W6s*(-N1*cotB+y1/W6*W5+a*y1/rb2)*y1+y2*W8/rb/W6*&
       (one_over_W6*W5-y1**2/W6s*W5/rb-y1**2/W6*a/rbc+a/rb2-2.0d0*a*y1**&
       2/rb2s)-y2*W8/rbc/W7*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+&
       (N3)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/rb2)*y1-y2*W8/rb/&
       W7s*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+(N3)*(rb*sinB-y1)*&
       cosB)-a*y3b*cosB*cotB/rb2)*(y1/rb-sinB)+y2*W8/rb/W7*(-cosB/&
       W7s*(W1*(N1*cosB-a/rb)*cotB+(N3)*(rb*sinB-y1)*cosB)*(y1/&
       rb-sinB)+cosB/W7*(one_over_rb*cosB*y1*(N1*cosB-a/rb)*cotB+W1*a/rb**&
       3*y1*cotB+(N3)*(one_over_rb*sinB*y1-1.0d0)*cosB)+2.0d0*a*y3b*cosB*cotB/&
       rb2s*y1))/pi/(N4))+&
       b2*(one_quarter*(N1*(((N3)*cotB**2+nu)/rb*y1/W6-((N3)*cotB**2+1.0d0)*&
       cosB*(y1/rb-sinB)/W7)-N1/W6s*(-N1*y1*cotB+nu*y3b-a+a*y1*&
       cotB/rb+y1**2/W6*W4)/rb*y1+N1/W6*(-N1*cotB+a*cotB/rb-a*&
       y1**2*cotB/rbc+2.0d0*y1/W6*W4-y1**3/W6s*W4/rb-y1**3/W6*a/&
       rbc)+N1*cotB/W7s*(z1b*cosB-a*(rb*sinB-y1)/rb/cosB)*(y1/&
       rb-sinB)-N1*cotB/W7*(cosB**2-a*(one_over_rb*sinB*y1-1.0d0)/rb/cosB+a*&
       (rb*sinB-y1)/rbc/cosB*y1)-a*W8*cotB/rbc+3*a*y1**2*W8*&
       cotB/rb**5-W8/W6s*(two_nu+one_over_rb*(N1*y1*cotB+a)-y1**2/rb/W6*&
       W5-a*y1**2/rbc)/rb*y1+W8/W6*(-one_over_rb_p3*(N1*y1*cotB+a)*y1+&
       one_over_rb*N1*cotB-2.0d0*y1/rb/W6*W5+y1**3/rbc/W6*W5+y1**3/rb2/&
       W6s*W5+y1**3/rb2s/W6*a-2.0d0*a/rbc*y1+3*a*y1**3/rb**5)-W8*&
       cotB/W7s*(-cosB*sinB+a*y1*y3b/rbc/cosB+(rb*sinB-y1)/rb*&
       ((N3)*cosB-w1_over_w7*W9))*(y1/rb-sinB)+W8*cotB/W7*(a*y3b/&
       rbc/cosB-3*a*y1**2*y3b/rb**5/cosB+(one_over_rb*sinB*y1-1.0d0)/rb*&
       ((N3)*cosB-w1_over_w7*W9)-(rb*sinB-y1)/rbc*((N3)*cosB-W1/&
       W7*W9)*y1+(rb*sinB-y1)/rb*(-one_over_rb*cosB*y1/W7*W9+w1_over_w7s*&
       W9*(y1/rb-sinB)+w1_over_w7*a/rbc/cosB*y1)))/pi/(N4))+&
       b3*(one_quarter*(N1*(-y2/W6s*(1.0d0+a/rb)/rb*y1-y2/W6*a/rbc*y1+y2*&
       cosB/W7s*W2*(y1/rb-sinB)+y2*cosB/W7*a/rbc*y1)+y2*W8/&
       rbc*(a/rb2+one_over_W6)*y1-y2*W8/rb*(-2.0d0*a/rb2s*y1-one_over_W6_p2/&
       rb*y1)-y2*W8*cosB/rbc/W7*(w1_over_w7*W2+a*y3b/rb2)*y1-y2*W8*&
       cosB/rb/W7s*(w1_over_w7*W2+a*y3b/rb2)*(y1/rb-sinB)+y2*W8*&
       cosB/rb/W7*(one_over_rb*cosB*y1/W7*W2-w1_over_w7s*W2*(y1/rb-sinB)-&
       w1_over_w7*a/rbc*y1-2.0d0*a*y3b/rb2s*y1))/pi/(N4));
  !print *,v11
  v22 = b1*(one_quarter*(N1*(((N3)*cotB**2-nu)/rb*y2/W6-((N3)*cotB**2+1.0d0-&
       two_nu)*cosB/rb*y2/W7)+N1/W6s*(y1*cotB*(1.0d0-W5)+nu*y3b-a+y2**&
       2/W6*W4)/rb*y2-N1/W6*(a*y1*cotB/rbc*y2+2.0d0*y2/W6*W4-y2**3&
       /W6s*W4/rb-y2**3/W6*a/rbc)+N1*z1b*cotB/W7s*W2/rb*&
       y2+N1*z1b*cotB/W7*a/rbc*y2+3*a*y2*W8*cotB/rb**5*y1-W8/&
       W6s*(-two_nu+one_over_rb*(N1*y1*cotB-a)+y2s/rb/W6*W5+a*y2s/&
       rbc)/rb*y2+W8/W6*(-one_over_rb_p3*(N1*y1*cotB-a)*y2+2.0d0*y2/rb/&
       W6*W5-y2**3/rbc/W6*W5-y2**3/rb2/W6s*W5-y2**3/rb2s/W6*&
       a+2.0d0*a/rbc*y2-3*a*y2**3/rb**5)-W8/W7s*(cosB**2-one_over_rb*(N1*&
       z1b*cotB+a*cosB)+a*y3b*z1b*cotB/rbc-one_over_rb/W7*(y2s*cosB**2-&
       a*z1b*cotB/rb*W1))/rb*y2+W8/W7*(one_over_rb_p3*(N1*z1b*cotB+a*&
       cosB)*y2-3*a*y3b*z1b*cotB/rb**5*y2+one_over_rb_p3/W7*(y2s*cosB**2-&
       a*z1b*cotB/rb*W1)*y2+one_over_rb2/W7s*(y2s*cosB**2-a*z1b*cotB/&
       rb*W1)*y2-one_over_rb/W7*(2.0d0*y2*cosB**2+a*z1b*cotB/rbc*W1*y2-a*&
       z1b*cotB/rb2*cosB*y2)))/pi/(N4))+&
       b2*(one_quarter*((N3)*N1*rFib_ry2*cotB**2+N1/W6*((W5-1.0d0)*cotB+y1/W6*&
       W4)-N1*y2s/W6s*((W5-1.0d0)*cotB+y1/W6*W4)/rb+N1*y2/W6*(-a/&
       rbc*y2*cotB-y1/W6s*W4/rb*y2-y2/W6*a/rbc*y1)-N1*cotB/&
       W7*W9+N1*y2s*cotB/W7s*W9/rb+N1*y2s*cotB/W7*a/rbc/&
       cosB-a*W8*cotB/rbc+3*a*y2s*W8*cotB/rb**5+W8/rb/W6*(N1*&
       cotB-two_nu*y1/W6-a*y1/rb*(one_over_rb+one_over_W6))-y2s*W8/rbc/W6*&
       (N1*cotB-two_nu*y1/W6-a*y1/rb*(one_over_rb+one_over_W6))-y2s*W8/rb2/W6s&
       *(N1*cotB-two_nu*y1/W6-a*y1/rb*(one_over_rb+one_over_W6))+y2*W8/rb/W6*&
       (two_nu*y1/W6s/rb*y2+a*y1/rbc*(one_over_rb+one_over_W6)*y2-a*y1/rb*&
       (-one_over_rb_p3*y2-one_over_W6_p2/rb*y2))+W8*cotB/rb/W7*((N2)*cosB+&
       w1_over_w7*W9+a*y3b/rb2/cosB)-y2s*W8*cotB/rbc/W7*((N2)*&
       cosB+w1_over_w7*W9+a*y3b/rb2/cosB)-y2s*W8*cotB/rb2/W7s*((-2+&
       two_nu)*cosB+w1_over_w7*W9+a*y3b/rb2/cosB)+y2*W8*cotB/rb/W7*(1/&
       rb*cosB*y2/W7*W9-w1_over_w7s*W9/rb*y2-w1_over_w7*a/rbc/cosB*y2-&
       2.0d0*a*y3b/rb2s/cosB*y2))/pi/(N4))+&
       b3*(one_quarter*(N1*(-sinB/rb*y2/W7+y2/W6s*(1.0d0+a/rb)/rb*y1+y2/W6*&
       a/rbc*y1-z1b/W7s*W2/rb*y2-z1b/W7*a/rbc*y2)-y2*W8/&
       rbc*(a/rb2+one_over_W6)*y1+y1*W8/rb*(-2.0d0*a/rb2s*y2-one_over_W6_p2/&
       rb*y2)+W8/W7s*(sinB*(cosB-a/rb)+z1b/rb*(1.0d0+a*y3b/rb2)-1.0d0/&
       rb/W7*(y2s*cosB*sinB-a*z1b/rb*W1))/rb*y2-W8/W7*(sinB*a/&
       rbc*y2-z1b/rbc*(1.0d0+a*y3b/rb2)*y2-2.0d0*z1b/rb**5*a*y3b*y2+&
       one_over_rb**3/W7*(y2s*cosB*sinB-a*z1b/rb*W1)*y2+one_over_rb2/W7s*&
       (y2s*cosB*sinB-a*z1b/rb*W1)*y2-one_over_rb/W7*(2.0d0*y2*cosB*sinB+a*&
       z1b/rbc*W1*y2-a*z1b/rb2*cosB*y2)))/pi/(N4));

  v33 = b1*(one_quarter*((N3)*(N1*rFib_ry3*cotB-y2/W6s*W5*(y3b/rb+1.0d0)-&
       one_half*y2/W6*a/rbc*2.0d0*y3b+y2*cosB/W7s*W2*W3+one_half*y2*cosB/W7*&
       a/rbc*2.0d0*y3b)+y2/rb*(two_nu/W6+a/rb2)-one_half*y2*W8/rbc*(2.0d0*&
       nu/W6+a/rb2)*2.0d0*y3b+y2*W8/rb*(-two_nu/W6s*(y3b/rb+1.0d0)-a/&
       rb2s*2*y3b)+y2*cosB/rb/W7*(1.0d0-two_nu-w1_over_w7*W2-a*y3b/rb2)-&
       one_half*y2*W8*cosB/rbc/W7*(1.0d0-two_nu-w1_over_w7*W2-a*y3b/rb2)*2.0d0*&
       y3b-y2*W8*cosB/rb/W7s*(1.0d0-two_nu-w1_over_w7*W2-a*y3b/rb2)*W3+y2*&
       W8*cosB/rb/W7*(-(cosB*y3b/rb+1.0d0)/W7*W2+w1_over_w7s*W2*W3+one_half*&
       w1_over_w7*a/rbc*2.0d0*y3b-a/rb2+a*y3b/rb2s*2.0d0*y3b))/pi/(N4))+&
       b2*(one_quarter*((N2)*N1*cotB*((y3b/rb+1.0d0)/W6-cosB*W3/W7)+(N3)*&
       y1/W6s*W5*(y3b/rb+1.0d0)+one_half*(N3)*y1/W6*a/rbc*2.0d0*y3b+(2.0d0-&
       two_nu)*sinB/W7*W2-(N3)*z1b/W7s*W2*W3-one_half*(N3)*z1b/&
       W7*a/rbc*2.0d0*y3b+one_over_rb*(N1*cotB-two_nu*y1/W6-a*y1/rb2)-one_half*&
       W8/rbc*(N1*cotB-two_nu*y1/W6-a*y1/rb2)*2.0d0*y3b+W8/rb*(two_nu*&
       y1/W6s*(y3b/rb+1.0d0)+a*y1/rb2s*2*y3b)-one_over_w7*(cosB*sinB+W1*&
       cotB/rb*((N3)*cosB-w1_over_w7)+a/rb*(sinB-y3b*z1b/rb2-z1b*&
       W1/rb/W7))+W8/W7s*(cosB*sinB+W1*cotB/rb*((N3)*cosB-W1/&
       W7)+a/rb*(sinB-y3b*z1b/rb2-z1b*W1/rb/W7))*W3-W8/W7*((cosB*&
       y3b/rb+1.0d0)*cotB/rb*((N3)*cosB-w1_over_w7)-one_half*W1*cotB/rbc*&
       ((N3)*cosB-w1_over_w7)*2.0d0*y3b+W1*cotB/rb*(-(cosB*y3b/rb+1.0d0)/W7+&
       w1_over_w7s*W3)-one_half*a/rbc*(sinB-y3b*z1b/rb2-z1b*W1/rb/W7)*&
       2*y3b+a/rb*(-z1b/rb2-y3b*sinB/rb2+y3b*z1b/rb2s*2.0d0*y3b-&
       sinB*W1/rb/W7-z1b*(cosB*y3b/rb+1.0d0)/rb/W7+one_half*z1b*W1/rbc/&
       W7*2.0d0*y3b+z1b*W1/rb/W7s*W3)))/pi/(N4))+&
       b3*(one_quarter*((N3)*rFib_ry3-(N3)*y2*sinB/W7s*W2*W3-one_half*&
       (N3)*y2*sinB/W7*a/rbc*2.0d0*y3b+y2*sinB/rb/W7*(1.0d0+w1_over_w7*&
       W2+a*y3b/rb2)-one_half*y2*W8*sinB/rbc/W7*(1.0d0+w1_over_w7*W2+a*y3b/&
       rb2)*2.0d0*y3b-y2*W8*sinB/rb/W7s*(1.0d0+w1_over_w7*W2+a*y3b/rb2)*W3+&
       y2*W8*sinB/rb/W7*((cosB*y3b/rb+1.0d0)/W7*W2-w1_over_w7s*W2*W3-&
       one_half*w1_over_w7*a/rbc*2.0d0*y3b+a/rb2-a*y3b/rb2s*2*y3b))/pi/(N4));

  v12 = b1/2.0d0*(one_quarter*((N2)*N1*rFib_ry2*cotB**2+N1/W6*((1.0d0-W5)*cotB-y1/&
       W6*W4)-N1*y2s/W6s*((1.0d0-W5)*cotB-y1/W6*W4)/rb+N1*y2/W6*&
       (a/rbc*y2*cotB+y1/W6s*W4/rb*y2+y2/W6*a/rbc*y1)+N1*&
       cosB*cotB/W7*W2-N1*y2s*cosB*cotB/W7s*W2/rb-N1*y2s*cosB*&
       cotB/W7*a/rbc+a*W8*cotB/rbc-3*a*y2s*W8*cotB/rb**5+W8/&
       rb/W6*(-N1*cotB+y1/W6*W5+a*y1/rb2)-y2s*W8/rbc/W6*(-N1*&
       cotB+y1/W6*W5+a*y1/rb2)-y2s*W8/rb2/W6s*(-N1*cotB+y1/&
       W6*W5+a*y1/rb2)+y2*W8/rb/W6*(-y1/W6s*W5/rb*y2-y2/W6*&
       a/rbc*y1-2.0d0*a*y1/rb2s*y2)+W8/rb/W7*(cosB/W7*(W1*(N1*&
       cosB-a/rb)*cotB+(N3)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/&
       rb2)-y2s*W8/rbc/W7*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+(2.0d0-&
       two_nu)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/rb2)-y2s*W8/rb2/&
       W7s*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+(N3)*(rb*sinB-y1)*&
       cosB)-a*y3b*cosB*cotB/rb2)+y2*W8/rb/W7*(-cosB/W7s*(W1*&
       (N1*cosB-a/rb)*cotB+(N3)*(rb*sinB-y1)*cosB)/rb*y2+cosB/&
       W7*(one_over_rb*cosB*y2*(N1*cosB-a/rb)*cotB+W1*a/rbc*y2*cotB+(2.0d0-2.0d0*&
       nu)/rb*sinB*y2*cosB)+2.0d0*a*y3b*cosB*cotB/rb2s*y2))/pi/(N4))+&
       b2/2.0d0*(one_quarter*(N1*(((N3)*cotB**2+nu)/rb*y2/W6-((N3)*cotB**2+1.0d0)*&
       cosB/rb*y2/W7)-N1/W6s*(-N1*y1*cotB+nu*y3b-a+a*y1*cotB/rb+&
       y1**2/W6*W4)/rb*y2+N1/W6*(-a*y1*cotB/rbc*y2-y1**2/W6s&
       *W4/rb*y2-y1**2/W6*a/rbc*y2)+N1*cotB/W7s*(z1b*cosB-a*&
       (rb*sinB-y1)/rb/cosB)/rb*y2-N1*cotB/W7*(-a/rb2*sinB*y2/&
       cosB+a*(rb*sinB-y1)/rbc/cosB*y2)+3*a*y2*W8*cotB/rb**5*y1-&
       W8/W6s*(two_nu+one_over_rb*(N1*y1*cotB+a)-y1**2/rb/W6*W5-a*y1**2/&
       rbc)/rb*y2+W8/W6*(-one_over_rb_p3*(N1*y1*cotB+a)*y2+y1**2/rb**&
       3/W6*W5*y2+y1**2/rb2/W6s*W5*y2+y1**2/rb2s/W6*a*y2+3*&
       a*y1**2/rb**5*y2)-W8*cotB/W7s*(-cosB*sinB+a*y1*y3b/rbc/&
       cosB+(rb*sinB-y1)/rb*((N3)*cosB-w1_over_w7*W9))/rb*y2+W8*cotB/&
       W7*(-3*a*y1*y3b/rb**5/cosB*y2+one_over_rb2*sinB*y2*((N3)*cosB-&
       w1_over_w7*W9)-(rb*sinB-y1)/rbc*((N3)*cosB-w1_over_w7*W9)*y2+(rb*&
       sinB-y1)/rb*(-one_over_rb*cosB*y2/W7*W9+w1_over_w7s*W9/rb*y2+w1_over_w7*&
       a/rbc/cosB*y2)))/pi/(N4))+&
       b3/2.0d0*(one_quarter*(N1*(one_over_W6*(1.0d0+a/rb)-y2s/W6s*(1.0d0+a/rb)/rb-y2s/&
       W6*a/rbc-cosB/W7*W2+y2s*cosB/W7s*W2/rb+y2s*cosB/W7*&
       a/rbc)-W8/rb*(a/rb2+one_over_W6)+y2s*W8/rbc*(a/rb2+one_over_W6)-&
       y2*W8/rb*(-2.0d0*a/rb2s*y2-one_over_W6_p2/rb*y2)+W8*cosB/rb/W7*&
       (w1_over_w7*W2+a*y3b/rb2)-y2s*W8*cosB/rbc/W7*(w1_over_w7*W2+a*&
       y3b/rb2)-y2s*W8*cosB/rb2/W7s*(w1_over_w7*W2+a*y3b/rb2)+y2*&
       W8*cosB/rb/W7*(one_over_rb*cosB*y2/W7*W2-w1_over_w7s*W2/rb*y2-W1/&
       W7*a/rbc*y2-2.0d0*a*y3b/rb2s*y2))/pi/(N4))+&
       b1/2.0d0*(one_quarter*(N1*(((N3)*cotB**2-nu)/rb*y1/W6-((N3)*cotB**2+1.0d0-&
       two_nu)*cosB*(y1/rb-sinB)/W7)+N1/W6s*(y1*cotB*(1.0d0-W5)+nu*y3b-&
       a+y2s/W6*W4)/rb*y1-N1/W6*((1.0d0-W5)*cotB+a*y1**2*cotB/rbc-&
       y2s/W6s*W4/rb*y1-y2s/W6*a/rbc*y1)-N1*cosB*cotB/W7*&
       W2+N1*z1b*cotB/W7s*W2*(y1/rb-sinB)+N1*z1b*cotB/W7*a/rb**&
       3*y1-a*W8*cotB/rbc+3*a*y1**2*W8*cotB/rb**5-W8/W6s*(-2.0d0*&
       nu+one_over_rb*(N1*y1*cotB-a)+y2s/rb/W6*W5+a*y2s/rbc)/rb*&
       y1+W8/W6*(-one_over_rb_p3*(N1*y1*cotB-a)*y1+one_over_rb*N1*cotB-y2s/&
       rbc/W6*W5*y1-y2s/rb2/W6s*W5*y1-y2s/rb2s/W6*a*y1-&
       3*a*y2s/rb**5*y1)-W8/W7s*(cosB**2-one_over_rb*(N1*z1b*cotB+a*&
       cosB)+a*y3b*z1b*cotB/rbc-one_over_rb/W7*(y2s*cosB**2-a*z1b*cotB/&
       rb*W1))*(y1/rb-sinB)+W8/W7*(one_over_rb_p3*(N1*z1b*cotB+a*cosB)*&
       y1-one_over_rb*N1*cosB*cotB+a*y3b*cosB*cotB/rbc-3*a*y3b*z1b*cotB/&
       rb**5*y1+one_over_rb_p3/W7*(y2s*cosB**2-a*z1b*cotB/rb*W1)*y1+1.0d0/&
       rb/W7s*(y2s*cosB**2-a*z1b*cotB/rb*W1)*(y1/rb-sinB)-one_over_rb/&
       W7*(-a*cosB*cotB/rb*W1+a*z1b*cotB/rbc*W1*y1-a*z1b*cotB/&
       rb2*cosB*y1)))/pi/(N4))+&
       b2/2.0d0*(one_quarter*((N3)*N1*rFib_ry1*cotB**2-N1*y2/W6s*((W5-1.0d0)*cotB+&
       y1/W6*W4)/rb*y1+N1*y2/W6*(-a/rbc*y1*cotB+one_over_W6*W4-y1**&
       2/W6s*W4/rb-y1**2/W6*a/rbc)+N1*y2*cotB/W7s*W9*(y1/&
       rb-sinB)+N1*y2*cotB/W7*a/rbc/cosB*y1+3*a*y2*W8*cotB/rb**&
       5*y1-y2*W8/rbc/W6*(N1*cotB-two_nu*y1/W6-a*y1/rb*(one_over_rb+1.0d0/&
       W6))*y1-y2*W8/rb2/W6s*(N1*cotB-two_nu*y1/W6-a*y1/rb*(1/&
       rb+one_over_W6))*y1+y2*W8/rb/W6*(-two_nu/W6+two_nu*y1**2/W6s/rb-a/&
       rb*(one_over_rb+one_over_W6)+a*y1**2/rbc*(one_over_rb+one_over_W6)-a*y1/rb*(-1.0d0/&
       rbc*y1-one_over_W6_p2/rb*y1))-y2*W8*cotB/rbc/W7*((N2)*&
       cosB+w1_over_w7*W9+a*y3b/rb2/cosB)*y1-y2*W8*cotB/rb/W7s*((-2+&
       two_nu)*cosB+w1_over_w7*W9+a*y3b/rb2/cosB)*(y1/rb-sinB)+y2*W8*&
       cotB/rb/W7*(one_over_rb*cosB*y1/W7*W9-w1_over_w7s*W9*(y1/rb-sinB)-&
       w1_over_w7*a/rbc/cosB*y1-2.0d0*a*y3b/rb2s/cosB*y1))/pi/(N4))+&
       b3/2.0d0*(one_quarter*(N1*(-sinB*(y1/rb-sinB)/W7-one_over_W6*(1.0d0+a/rb)+y1**2/W6s&
       *(1.0d0+a/rb)/rb+y1**2/W6*a/rbc+cosB/W7*W2-z1b/W7s*W2*&
       (y1/rb-sinB)-z1b/W7*a/rbc*y1)+W8/rb*(a/rb2+one_over_W6)-y1**2*&
       W8/rbc*(a/rb2+one_over_W6)+y1*W8/rb*(-2.0d0*a/rb2s*y1-one_over_W6_p2/&
       rb*y1)+W8/W7s*(sinB*(cosB-a/rb)+z1b/rb*(1.0d0+a*y3b/rb2)-1.0d0/&
       rb/W7*(y2s*cosB*sinB-a*z1b/rb*W1))*(y1/rb-sinB)-W8/W7*&
       (sinB*a/rbc*y1+cosB/rb*(1.0d0+a*y3b/rb2)-z1b/rbc*(1.0d0+a*y3b/&
       rb2)*y1-2.0d0*z1b/rb**5*a*y3b*y1+one_over_rb_p3/W7*(y2s*cosB*sinB-a*&
       z1b/rb*W1)*y1+one_over_rb/W7s*(y2s*cosB*sinB-a*z1b/rb*W1)*&
       (y1/rb-sinB)-one_over_rb/W7*(-a*cosB/rb*W1+a*z1b/rbc*W1*y1-a*&
       z1b/rb2*cosB*y1)))/pi/(N4));

  v13 = b1/2.0d0*(one_quarter*((N2)*N1*rFib_ry3*cotB**2-N1*y2/W6s*((1.0d0-W5)*&
       cotB-y1/W6*W4)*(y3b/rb+1.0d0)+N1*y2/W6*(one_half*a/rbc*2.0d0*y3b*cotB+&
       y1/W6s*W4*(y3b/rb+1.0d0)+one_half*y1/W6*a/rbc*2.0d0*y3b)-N1*y2*cosB*&
       cotB/W7s*W2*W3-one_half*N1*y2*cosB*cotB/W7*a/rbc*2.0d0*y3b+a/&
       rbc*y2*cotB-3/2.0d0*a*y2*W8*cotB/rb**5*2.0d0*y3b+y2/rb/W6*(-N1*&
       cotB+y1/W6*W5+a*y1/rb2)-one_half*y2*W8/rbc/W6*(-N1*cotB+y1/&
       W6*W5+a*y1/rb2)*2.0d0*y3b-y2*W8/rb/W6s*(-N1*cotB+y1/W6*W5+&
       a*y1/rb2)*(y3b/rb+1.0d0)+y2*W8/rb/W6*(-y1/W6s*W5*(y3b/rb+&
       1)-one_half*y1/W6*a/rbc*2.0d0*y3b-a*y1/rb2s*2*y3b)+y2/rb/W7*&
       (cosB/W7*(W1*(N1*cosB-a/rb)*cotB+(N3)*(rb*sinB-y1)*cosB)-&
       a*y3b*cosB*cotB/rb2)-one_half*y2*W8/rbc/W7*(cosB/W7*(W1*(N1*&
       cosB-a/rb)*cotB+(N3)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/&
       rb2)*2.0d0*y3b-y2*W8/rb/W7s*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+&
       (N3)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/rb2)*W3+y2*W8/rb/&
       W7*(-cosB/W7s*(W1*(N1*cosB-a/rb)*cotB+(N3)*(rb*sinB-y1)*&
       cosB)*W3+cosB/W7*((cosB*y3b/rb+1.0d0)*(N1*cosB-a/rb)*cotB+one_half*W1*&
       a/rbc*2.0d0*y3b*cotB+one_half*(N3)/rb*sinB*2.0d0*y3b*cosB)-a*cosB*&
       cotB/rb2+a*y3b*cosB*cotB/rb2s*2*y3b))/pi/(N4))+&
       b2/2.0d0*(one_quarter*(N1*(((N3)*cotB**2+nu)*(y3b/rb+1.0d0)/W6-((N3)*cotB**&
       2+1.0d0)*cosB*W3/W7)-N1/W6s*(-N1*y1*cotB+nu*y3b-a+a*y1*cotB/&
       rb+y1**2/W6*W4)*(y3b/rb+1.0d0)+N1/W6*(nu-one_half*a*y1*cotB/rbc*2.0d0*&
       y3b-y1**2/W6s*W4*(y3b/rb+1.0d0)-one_half*y1**2/W6*a/rbc*2.0d0*y3b)+&
       N1*cotB/W7s*(z1b*cosB-a*(rb*sinB-y1)/rb/cosB)*W3-N1*cotB/&
       W7*(cosB*sinB-one_half*a/rb2*sinB*2.0d0*y3b/cosB+one_half*a*(rb*sinB-y1)/&
       rbc/cosB*2.0d0*y3b)-a/rbc*y1*cotB+3/2.0d0*a*y1*W8*cotB/rb**5*2.0d0*&
       y3b+one_over_W6*(two_nu+one_over_rb*(N1*y1*cotB+a)-y1**2/rb/W6*W5-a*y1**2/&
       rbc)-W8/W6s*(two_nu+one_over_rb*(N1*y1*cotB+a)-y1**2/rb/W6*W5-a*&
       y1**2/rbc)*(y3b/rb+1.0d0)+W8/W6*(-one_half/rbc*(N1*y1*cotB+a)*2.0d0*&
       y3b+one_half*y1**2/rbc/W6*W5*2.0d0*y3b+y1**2/rb/W6s*W5*(y3b/rb+&
       1)+one_half*y1**2/rb2s/W6*a*2.0d0*y3b+3/2.0d0*a*y1**2/rb**5*2.0d0*y3b)+&
       cotB/W7*(-cosB*sinB+a*y1*y3b/rbc/cosB+(rb*sinB-y1)/rb*((2.0d0-&
       two_nu)*cosB-w1_over_w7*W9))-W8*cotB/W7s*(-cosB*sinB+a*y1*y3b/rb**&
       3/cosB+(rb*sinB-y1)/rb*((N3)*cosB-w1_over_w7*W9))*W3+W8*cotB/&
       W7*(a/rbc/cosB*y1-3/2.0d0*a*y1*y3b/rb**5/cosB*2.0d0*y3b+one_half/&
       rb2*sinB*2.0d0*y3b*((N3)*cosB-w1_over_w7*W9)-one_half*(rb*sinB-y1)/rb**&
       3*((N3)*cosB-w1_over_w7*W9)*2.0d0*y3b+(rb*sinB-y1)/rb*(-(cosB*y3b/&
       rb+1.0d0)/W7*W9+w1_over_w7s*W9*W3+one_half*w1_over_w7*a/rbc/cosB*2.0d0*&
       y3b)))/pi/(N4))+&
       b3/2.0d0*(one_quarter*(N1*(-y2/W6s*(1.0d0+a/rb)*(y3b/rb+1.0d0)-one_half*y2/W6*a/&
       rbc*2.0d0*y3b+y2*cosB/W7s*W2*W3+one_half*y2*cosB/W7*a/rbc*2.0d0*&
       y3b)-y2/rb*(a/rb2+one_over_W6)+one_half*y2*W8/rbc*(a/rb2+one_over_W6)*2.0d0*&
       y3b-y2*W8/rb*(-a/rb2s*2*y3b-one_over_W6_p2*(y3b/rb+1.0d0))+y2*cosB/&
       rb/W7*(w1_over_w7*W2+a*y3b/rb2)-one_half*y2*W8*cosB/rbc/W7*(W1/&
       W7*W2+a*y3b/rb2)*2.0d0*y3b-y2*W8*cosB/rb/W7s*(w1_over_w7*W2+a*&
       y3b/rb2)*W3+y2*W8*cosB/rb/W7*((cosB*y3b/rb+1.0d0)/W7*W2-W1/&
       W7s*W2*W3-1.0d0/2.d0*w1_over_w7*a/rbc*2.0d0*y3b+a/rb2-a*y3b/rb2s*2*&
       y3b))/pi/(N4))+&
       b1/2*(one_quarter*((N3)*(N1*rFib_ry1*cotB-y1/W6s*W5/rb*y2-y2/W6*&
       a/rbc*y1+y2*cosB/W7s*W2*(y1/rb-sinB)+y2*cosB/W7*a/rb**&
       3*y1)-y2*W8/rbc*(two_nu/W6+a/rb2)*y1+y2*W8/rb*(-two_nu/W6s&
       /rb*y1-2.0d0*a/rb2s*y1)-y2*W8*cosB/rbc/W7*(1.0d0-two_nu-w1_over_w7*&
       W2-a*y3b/rb2)*y1-y2*W8*cosB/rb/W7s*(1.0d0-two_nu-w1_over_w7*W2-a*&
       y3b/rb2)*(y1/rb-sinB)+y2*W8*cosB/rb/W7*(-one_over_rb*cosB*y1/W7*&
       W2+w1_over_w7s*W2*(y1/rb-sinB)+w1_over_w7*a/rbc*y1+2.0d0*a*y3b/rb2**&
       2*y1))/pi/(N4))+&
       b2/2*(one_quarter*((N2)*N1*cotB*(one_over_rb*y1/W6-cosB*(y1/rb-sinB)/W7)-&
       (N3)/W6*W5+(N3)*y1**2/W6s*W5/rb+(N3)*y1**2/W6*&
       a/rbc+(N3)*cosB/W7*W2-(N3)*z1b/W7s*W2*(y1/rb-&
       sinB)-(N3)*z1b/W7*a/rbc*y1-W8/rbc*(N1*cotB-two_nu*y1/&
       W6-a*y1/rb2)*y1+W8/rb*(-two_nu/W6+two_nu*y1**2/W6s/rb-a/rb2+&
       2*a*y1**2/rb2s)+W8/W7s*(cosB*sinB+W1*cotB/rb*((N3)*&
       cosB-w1_over_w7)+a/rb*(sinB-y3b*z1b/rb2-z1b*W1/rb/W7))*(y1/rb-&
       sinB)-W8/W7*(one_over_rb2*cosB*y1*cotB*((N3)*cosB-w1_over_w7)-W1*&
       cotB/rbc*((N3)*cosB-w1_over_w7)*y1+W1*cotB/rb*(-one_over_rb*cosB*&
       y1/W7+w1_over_w7s*(y1/rb-sinB))-a/rbc*(sinB-y3b*z1b/rb2-&
       z1b*W1/rb/W7)*y1+a/rb*(-y3b*cosB/rb2+2.0d0*y3b*z1b/rb2s*y1-&
       cosB*W1/rb/W7-z1b/rb2*cosB*y1/W7+z1b*W1/rbc/W7*y1+z1b*&
       W1/rb/W7s*(y1/rb-sinB))))/pi/(N4))+&
       b3/2*(one_quarter*((N3)*rFib_ry1-(N3)*y2*sinB/W7s*W2*(y1/rb-&
       sinB)-(N3)*y2*sinB/W7*a/rbc*y1-y2*W8*sinB/rbc/W7*(1.0d0+&
       w1_over_w7*W2+a*y3b/rb2)*y1-y2*W8*sinB/rb/W7s*(1.0d0+w1_over_w7*W2+&
       a*y3b/rb2)*(y1/rb-sinB)+y2*W8*sinB/rb/W7*(one_over_rb*cosB*y1/&
       W7*W2-w1_over_w7s*W2*(y1/rb-sinB)-w1_over_w7*a/rbc*y1-2.0d0*a*y3b/&
       rb2s*y1))/pi/(N4));

  v23 = b1/2.0d0*(one_quarter*(N1*(((N3)*cotB**2-nu)*(y3b/rb+1.0d0)/W6-((N3)*&
       cotB**2+1.0d0-two_nu)*cosB*W3/W7)+N1/W6s*(y1*cotB*(1.0d0-W5)+nu*y3b-a+&
       y2s/W6*W4)*(y3b/rb+1.0d0)-N1/W6*(one_half*a*y1*cotB/rbc*2.0d0*y3b+&
       nu-y2s/W6s*W4*(y3b/rb+1.0d0)-one_half*y2s/W6*a/rbc*2.0d0*y3b)-N1*&
       sinB*cotB/W7*W2+N1*z1b*cotB/W7s*W2*W3+one_half*N1*z1b*cotB/W7*&
       a/rbc*2.0d0*y3b-a/rbc*y1*cotB+3/2.0d0*a*y1*W8*cotB/rb**5*2.0d0*y3b+&
       1/W6*(-two_nu+one_over_rb*(N1*y1*cotB-a)+y2s/rb/W6*W5+a*y2s/&
       rbc)-W8/W6s*(-two_nu+one_over_rb*(N1*y1*cotB-a)+y2s/rb/W6*W5+&
       a*y2s/rbc)*(y3b/rb+1.0d0)+W8/W6*(-one_half/rbc*(N1*y1*cotB-a)*&
       2*y3b-one_half*y2s/rbc/W6*W5*2.0d0*y3b-y2s/rb/W6s*W5*(y3b/&
       rb+1.0d0)-one_half*y2s/rb2s/W6*a*2.0d0*y3b-3/2.0d0*a*y2s/rb**5*2.0d0*y3b)+&
       1/W7*(cosB**2-one_over_rb*(N1*z1b*cotB+a*cosB)+a*y3b*z1b*cotB/rb**&
       3-one_over_rb/W7*(y2s*cosB**2-a*z1b*cotB/rb*W1))-W8/W7s*(cosB**2-&
       1/rb*(N1*z1b*cotB+a*cosB)+a*y3b*z1b*cotB/rbc-one_over_rb/W7*&
       (y2s*cosB**2-a*z1b*cotB/rb*W1))*W3+W8/W7*(one_half/rbc*(N1*&
       z1b*cotB+a*cosB)*2.0d0*y3b-one_over_rb*N1*sinB*cotB+a*z1b*cotB/rbc+a*&
       y3b*sinB*cotB/rbc-3/2.0d0*a*y3b*z1b*cotB/rb**5*2.0d0*y3b+one_half/rbc&
       /W7*(y2s*cosB**2-a*z1b*cotB/rb*W1)*2.0d0*y3b+one_over_rb/W7s*(y2**&
       2*cosB**2-a*z1b*cotB/rb*W1)*W3-one_over_rb/W7*(-a*sinB*cotB/rb*W1+&
       one_half*a*z1b*cotB/rbc*W1*2.0d0*y3b-a*z1b*cotB/rb*(cosB*y3b/rb+&
       1.0d0))))/pi/(N4))+&
       b2/2.0d0*(one_quarter*((N3)*N1*rFib_ry3*cotB**2-N1*y2/W6s*((W5-1.0d0)*cotB+&
       y1/W6*W4)*(y3b/rb+1.0d0)+N1*y2/W6*(-one_half*a/rbc*2.0d0*y3b*cotB-y1/&
       W6s*W4*(y3b/rb+1.0d0)-one_half*y1/W6*a/rbc*2.0d0*y3b)+N1*y2*cotB/&
       W7s*W9*W3+one_half*N1*y2*cotB/W7*a/rbc/cosB*2.0d0*y3b-a/rbc*&
       y2*cotB+3.0d0/2.0d0*a*y2*W8*cotB/rb**5*2.0d0*y3b+y2/rb/W6*(N1*cotB-2.0d0*&
       nu*y1/W6-a*y1/rb*(one_over_rb+one_over_W6))-one_half*y2*W8/rbc/W6*(N1*&
       cotB-two_nu*y1/W6-a*y1/rb*(one_over_rb+one_over_W6))*2.0d0*y3b-y2*W8/rb/W6s&
       *(N1*cotB-two_nu*y1/W6-a*y1/rb*(one_over_rb+one_over_W6))*(y3b/rb+1.0d0)+y2*&
       W8/rb/W6*(two_nu*y1/W6s*(y3b/rb+1.0d0)+one_half*a*y1/rbc*(one_over_rb+&
       1.0d0/W6)*2.0d0*y3b-a*y1/rb*(-one_half/rbc*2.0d0*y3b-one_over_W6_p2*(y3b/rb+&
       1.0d0)))+y2*cotB/rb/W7*((N2)*cosB+w1_over_w7*W9+a*y3b/rb2/cosB)-&
       one_half*y2*W8*cotB/rbc/W7*((N2)*cosB+w1_over_w7*W9+a*y3b/&
       rb2/cosB)*2.0d0*y3b-y2*W8*cotB/rb/W7s*((N2)*cosB+w1_over_w7*&
       W9+a*y3b/rb2/cosB)*W3+y2*W8*cotB/rb/W7*((cosB*y3b/rb+1.0d0)/&
       W7*W9-w1_over_w7s*W9*W3-one_half*w1_over_w7*a/rbc/cosB*2.0d0*y3b+a/rb2/&
       cosB-a*y3b/rb2s/cosB*2.0d0*y3b))/pi/(N4))+&
       b3/2.0d0*(one_quarter*(N1*(-sinB*W3/W7+y1/W6s*(1.0d0+a/rb)*(y3b/rb+1.0d0)+&
       one_half*y1/W6*a/rbc*2.0d0*y3b+sinB/W7*W2-z1b/W7s*W2*W3-one_half*&
       z1b/W7*a/rbc*2.0d0*y3b)+y1/rb*(a/rb2+one_over_W6)-one_half*y1*W8/rb**&
       3*(a/rb2+one_over_W6)*2.0d0*y3b+y1*W8/rb*(-a/rb2s*2.0d0*y3b-one_over_W6_p2*&
       (y3b/rb+1.0d0))-one_over_w7*(sinB*(cosB-a/rb)+z1b/rb*(1.0d0+a*y3b/rb2)-1.0d0/&
       rb/W7*(y2s*cosB*sinB-a*z1b/rb*W1))+W8/W7s*(sinB*(cosB-&
       a/rb)+z1b/rb*(1.0d0+a*y3b/rb2)-one_over_rb/W7*(y2s*cosB*sinB-a*z1b/&
       rb*W1))*W3-W8/W7*(one_half*sinB*a/rbc*2.0d0*y3b+sinB/rb*(1.0d0+a*y3b/&
       rb2)-one_half*z1b/rbc*(1.0d0+a*y3b/rb2)*2.0d0*y3b+z1b/rb*(a/rb2-a*&
       y3b/rb2s*2.0d0*y3b)+one_half/rbc/W7*(y2s*cosB*sinB-a*z1b/rb*&
       W1)*2.0d0*y3b+one_over_rb/W7s*(y2s*cosB*sinB-a*z1b/rb*W1)*W3-1.0d0/&
       rb/W7*(-a*sinB/rb*W1+one_half*a*z1b/rbc*W1*2.0d0*y3b-a*z1b/rb*&
       (cosB*y3b/rb+1.0d0))))/pi/(N4))+&
       b1/2.0d0*(one_quarter*((N3)*(N1*rFib_ry2*cotB+one_over_W6*W5-y2s/W6s*W5/&
       rb-y2s/W6*a/rbc-cosB/W7*W2+y2s*cosB/W7s*W2/rb+y2s*&
       cosB/W7*a/rbc)+W8/rb*(two_nu/W6+a/rb2)-y2s*W8/rbc*(2.0d0*&
       nu/W6+a/rb2)+y2*W8/rb*(-two_nu/W6s/rb*y2-2.0d0*a/rb2s*y2)+&
       W8*cosB/rb/W7*(1.0d0-two_nu-w1_over_w7*W2-a*y3b/rb2)-y2s*W8*cosB/&
       rbc/W7*(1.0d0-two_nu-w1_over_w7*W2-a*y3b/rb2)-y2s*W8*cosB/rb2/W7s&
       *(1.0d0-two_nu-w1_over_w7*W2-a*y3b/rb2)+y2*W8*cosB/rb/W7*(-one_over_rb*&
       cosB*y2/W7*W2+w1_over_w7s*W2/rb*y2+w1_over_w7*a/rbc*y2+2.0d0*a*&
       y3b/rb2s*y2))/pi/(N4))+&
       b2/2.0d0*(one_quarter*((N2)*N1*cotB*(one_over_rb*y2/W6-cosB/rb*y2/W7)+(2.0d0-&
       two_nu)*y1/W6s*W5/rb*y2+(N3)*y1/W6*a/rbc*y2-(2.0d0-2.0d0*&
       nu)*z1b/W7s*W2/rb*y2-(N3)*z1b/W7*a/rbc*y2-W8/rb**&
       3*(N1*cotB-two_nu*y1/W6-a*y1/rb2)*y2+W8/rb*(two_nu*y1/W6s/&
       rb*y2+2.0d0*a*y1/rb2s*y2)+W8/W7s*(cosB*sinB+W1*cotB/rb*((2.0d0-&
       two_nu)*cosB-w1_over_w7)+a/rb*(sinB-y3b*z1b/rb2-z1b*W1/rb/W7))/&
       rb*y2-W8/W7*(one_over_rb2*cosB*y2*cotB*((N3)*cosB-w1_over_w7)-W1*&
       cotB/rbc*((N3)*cosB-w1_over_w7)*y2+W1*cotB/rb*(-cosB/rb*&
       y2/W7+w1_over_w7s/rb*y2)-a/rbc*(sinB-y3b*z1b/rb2-z1b*W1/&
       rb/W7)*y2+a/rb*(2.0d0*y3b*z1b/rb2s*y2-z1b/rb2*cosB*y2/W7+&
       z1b*W1/rbc/W7*y2+z1b*W1/rb2/W7s*y2)))/pi/(N4))+&
       b3/2.0d0*(one_quarter*((N3)*rFib_ry2+(N3)*sinB/W7*W2-(N3)*y2s*&
       sinB/W7s*W2/rb-(N3)*y2s*sinB/W7*a/rbc+W8*sinB/rb/&
       W7*(1.0d0+w1_over_w7*W2+a*y3b/rb2)-y2s*W8*sinB/rbc/W7*(1.0d0+W1/&
       W7*W2+a*y3b/rb2)-y2s*W8*sinB/rb2/W7s*(1.0d0+w1_over_w7*W2+a*&
       y3b/rb2)+y2*W8*sinB/rb/W7*(one_over_rb*cosB*y2/W7*W2-w1_over_w7s*&
       W2/rb*y2-w1_over_w7*a/rbc*y2-2.0d0*a*y3b/rb2s*y2))/pi/(N4));
end subroutine AngDisStrainFSC
