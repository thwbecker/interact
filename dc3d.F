c
c this routine always double precision
c
#define C_PREC real*8
#define EPS_DC3D_PREC 7.0D-15
c
c
c     This routine is based on the original source code as 
c     published in Okada, BSSA, 1992, distributed by Y. Okada 
c     on disk in summer 2000
c
c     
c     small modifications, mostly to improve speed by Thorsten Becker, 
c     becker@eps.harvard.edu
c
c     - floating point precision is set in include files, C-style
c       all references to real*4 are replaced with C_PREC
c
c     - introduced tmpd[1-5] variables to speed up computations
c
c     - changed value of EPS to reflect double precision (was 1.0d-6)
c
c     - added call to sincos functions
c
c     see also rcs logs!
c
c $Id: dc3d.F,v 2.25 2011/01/09 02:02:43 becker Exp $ 
c
c
c
c
      SUBROUTINE  DC3D0(ALPHA,X,Y,Z,DEPTH,DIP,POT1,POT2,POT3,POT4,      
     *               UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET) 
      IMPLICIT C_PREC (A-H,O-Z)                                         
      C_PREC   ALPHA,X,Y,Z,DEPTH,DIP,POT1,POT2,POT3,POT4,               
     *         UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ             
C                                                                       
C********************************************************************   
C*****                                                          *****   
C*****    DISPLACEMENT AND STRAIN AT DEPTH                      *****   
C*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****   
C*****                         CODED BY  Y.OKADA ... SEP.1991   *****   
C*****                         REVISED   Y.OKADA ... NOV.1991   *****   
C*****                                                          *****   
C********************************************************************   
C                                                                       
C***** INPUT                                                            
C*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)           
C*****   X,Y,Z : COORDINATE OF OBSERVING POINT                          
C*****   DEPTH : SOURCE DEPTH                                           
C*****   DIP   : DIP-ANGLE (DEGREE)                                     
C*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY        
C*****       POTENCY=(  MOMENT OF DOUBLE-COUPLE  )/MYU     FOR POT1,2   
C*****       POTENCY=(INTENSITY OF ISOTROPIC PART)/LAMBDA  FOR POT3     
C*****       POTENCY=(INTENSITY OF LINEAR DIPOLE )/MYU     FOR POT4     
C                                                                       
C***** OUTPUT                                                           
C*****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF POTENCY) /          
C*****               :                     (UNIT OF X,Y,Z,DEPTH)**2  )  
C*****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT= UNIT OF POTENCY) /          
C*****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH)**3  )  
C*****   UXZ,UYZ,UZZ : Z-DERIVATIVE                                     
C*****   IRET        : RETURN CODE  ( =0....NORMAL,   =1....SINGULAR )  
C                                                                       
      COMMON /C1/DUMMY(8),R,dummy2(15)
      DIMENSION  U(12),DUA(12),DUB(12),DUC(12)                          
      DATA  F0/0.D0/                                                    
C-----                                                                  
      IF(Z.GT.0.0) then
         print *,'POSITIVE Z WAS GIVEN IN SUB-DC3D0'
         print *,z
         stop
      endif
      DO 111 I=1,12                                                     
        U(I)=F0                                                         
        DUA(I)=F0                                                       
        DUB(I)=F0                                                       
        DUC(I)=F0                                                       
  111 CONTINUE                                                          
      AALPHA=ALPHA                                                      
      DDIP=DIP                                                          
      CALL DCCON0(AALPHA,DDIP)                                          
C======================================                                 
C=====  REAL-SOURCE CONTRIBUTION  =====                                 
C======================================                                 
      XX=X                                                              
      YY=Y                                                              
      ZZ=Z                                                              
      DD=DEPTH+Z                                                        
      CALL DCCON1(XX,YY,DD)                                             
      IF(R.EQ.F0) GO TO 99                                              
      PP1=POT1                                                          
      PP2=POT2                                                          
      PP3=POT3                                                          
      PP4=POT4                                                          
      CALL UA0(XX,YY,DD,PP1,PP2,PP3,PP4,DUA)                            
C-----                                                                  
      DO 222 I=1,12                                                     
        IF(I.LT.10) U(I)=U(I)-DUA(I)                                    
        IF(I.GE.10) U(I)=U(I)+DUA(I)                                    
  222 CONTINUE                                                          
C=======================================                                
C=====  IMAGE-SOURCE CONTRIBUTION  =====                                
C=======================================                                
      DD=DEPTH-Z                                                        
      CALL DCCON1(XX,YY,DD)                                             
      CALL UA0(XX,YY,DD,PP1,PP2,PP3,PP4,DUA)                            
      CALL UB0(XX,YY,DD,ZZ,PP1,PP2,PP3,PP4,DUB)                         
      CALL UC0(XX,YY,DD,ZZ,PP1,PP2,PP3,PP4,DUC)                         
C-----                                                                  
      DO 333 I=1,12                                                     
        DU=DUA(I)+DUB(I)+ZZ*DUC(I)                                      
        IF(I.GE.10) DU=DU+DUC(I-9)                                      
        U(I)=U(I)+DU                                                    
  333 CONTINUE                                                          
C=====                                                                  
      UX=U(1)                                                           
      UY=U(2)                                                           
      UZ=U(3)                                                           
      UXX=U(4)                                                          
      UYX=U(5)                                                          
      UZX=U(6)                                                          
      UXY=U(7)                                                          
      UYY=U(8)                                                          
      UZY=U(9)                                                          
      UXZ=U(10)                                                         
      UYZ=U(11)                                                         
      UZZ=U(12)                                                         
      IRET=0                                                            
      RETURN                                                            
C=======================================                                
C=====  IN CASE OF SINGULAR (R=0)  =====                                
C=======================================                                
   99 UX=F0                                                             
      UY=F0                                                             
      UZ=F0                                                             
      UXX=F0                                                            
      UYX=F0                                                            
      UZX=F0                                                            
      UXY=F0                                                            
      UYY=F0                                                            
      UZY=F0                                                            
      UXZ=F0                                                            
      UYZ=F0                                                            
      UZZ=F0                                                            
      IRET=1                                                            
      RETURN                                                            
      END                                                               
      SUBROUTINE  UA0(X,Y,D,POT1,POT2,POT3,POT4,U)                      
      IMPLICIT C_PREC (A-H,O-Z)      
c     added this, th
      C_PREC tmpd1,tmpd2,tmpd3,tmpd4,tmpd5

      DIMENSION U(12),DU(12)                                            
C                                                                       
C********************************************************************   
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-A)             *****   
C*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****   
C********************************************************************   
C                                                                       
C***** INPUT                                                            
C*****   X,Y,D : STATION COORDINATES IN FAULT SYSTEM                    
C*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY        
C***** OUTPUT                                                           
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     
C                                                                       
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  
      COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,     
     *           UY,VY,WY,UZ,VZ,WZ                                      
      DATA F0,F1,F3/0.D0,1.D0,3.D0/                                     
      DATA PI2/6.283185307179586D0/                                     
C-----                                                                  
      DO 111  I=1,12                                                    
  111 U(I)=F0                                                           
C======================================                                 
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 
C======================================                                 
      IF(POT1.NE.F0) THEN                                               
         tmpd4=alp1/r3
         tmpd1=tmpd4*x
         tmpd2=tmpd4*a3
         tmpd3=F3*X/R5
         tmpd4=alp2*qr
         DU( 1)= ALP1*Q/R3    +tmpd4*X2                                
         DU( 2)= tmpd1*SD +tmpd4*XY                                
         DU( 3)=-tmpd1*CD +tmpd4*X*D                               
         DU( 4)= X*QR*(-ALP1 +ALP2*(F1+A5) )    
         DU( 5)= tmpd2*SD +tmpd4*Y*A5                             
         DU( 6)=-tmpd2*CD +tmpd4*D*A5                             
         DU( 7)= ALP1*(SD/R3-Y*QR) +ALP2*F3*X2/R5*UY   
         DU( 8)= tmpd3*(-ALP1*Y*SD +ALP2*(Y*UY+Q) )                    
         DU( 9)= tmpd3*( ALP1*Y*CD +ALP2*D*UY )                        
         DU(10)= ALP1*(CD/R3+D*QR) +ALP2*F3*X2/R5*UZ  
         DU(11)= tmpd3*( ALP1*D*SD +ALP2*Y*UZ )                        
         DU(12)= tmpd3*(-ALP1*D*CD +ALP2*(D*UZ-Q) )                    
         tmpd5 = POT1/PI2
         DO 222 I=1,12     
 222        U(I) = U(I) + tmpd5 * DU(I)      
      ENDIF                                                             
C=================================== 
C=====  DIP-SLIP CONTRIBUTION  =====                                    
C===================================                                    
      IF(POT2.NE.F0) THEN 
         tmpd1=alp2*qr*p
         tmpd2=alp2*f3/r5
         tmpd3=alp1*f3*x

         DU( 1)=            tmpd1*X                                  
         DU( 2)= ALP1*S/R3 +tmpd1*Y                                  
         DU( 3)=-ALP1*T/R3 +tmpd1*D                                  
         DU( 4)=                 tmpd1*A5                            


        DU( 5)=-tmpd3*S/R5 -ALP2*Y*P*QRX                            
        DU( 6)= tmpd3*T/R5 -ALP2*D*P*QRX                            
        DU( 7)=                          tmpd2*X*VY                
        DU( 8)= ALP1*(S2D/R3-F3*Y*S/R5) +ALP2*(F3*Y/R5*VY+P*QR)         
        DU( 9)=-ALP1*(C2D/R3-F3*Y*T/R5) +tmpd2*D*VY                
        DU(10)=                          tmpd2*X*VZ                
        DU(11)= ALP1*(C2D/R3+F3*D*S/R5) +tmpd2*Y*VZ                
        DU(12)= ALP1*(S2D/R3-F3*D*T/R5) +ALP2*(F3*D/R5*VZ-P*QR)         
        tmpd5 = POT2/PI2
        DO 333 I=1,12                                                   
  333   U(I)=U(I) + tmpd5 * DU(I)                                        
      ENDIF                                                             
C========================================                               
C=====  TENSILE-FAULT CONTRIBUTION  =====                               
C========================================                               
      IF(POT3.NE.F0) THEN                                               
        tmpd1=alp1*f3*x
        tmpd2=alp2*q
        tmpd3=alp2*qr
        tmpd4=tmpd2*qr

        DU( 1)= ALP1*X/R3 -tmpd4*X                                  
        DU( 2)= ALP1*T/R3 -tmpd4*Y                                 
        DU( 3)= ALP1*S/R3 -tmpd4*D
        DU( 4)= ALP1*A3/R3     -tmpd4*A5                            
        DU( 5)=-tmpd1*T/R5 +tmpd2*Y*QRX                            
        DU( 6)=-tmpd1*S/R5 +tmpd2*D*QRX                            
        DU( 7)=-ALP1*F3*XY/R5           -tmpd3*X*WY                   
        DU( 8)= ALP1*(C2D/R3-F3*Y*T/R5) -tmpd3*(Y*WY+Q)
        DU( 9)= ALP1*(S2D/R3-F3*Y*S/R5) -tmpd3*D*WY                   
        DU(10)= ALP1*F3*X*D/R5          -tmpd3*X*WZ                   
        DU(11)=-ALP1*(S2D/R3-F3*D*T/R5) -tmpd3*Y*WZ                   
        DU(12)= ALP1*(C2D/R3+F3*D*S/R5) -tmpd3*(D*WZ-Q)
        tmpd5 = POT3/PI2
        DO 444 I=1,12                                                   
 444       U(I)=U(I) + tmpd5 * DU(I)                                        
      ENDIF                                                             
C=========================================                              
C=====  INFLATE SOURCE CONTRIBUTION  =====                              
C=========================================                              
      IF(POT4.NE.F0) THEN                                               
         tmpd1=alp1/r3
        DU( 1)=-tmpd1*X                                               
        DU( 2)=-tmpd1*Y                                               
        DU( 3)=-tmpd1*D                                               
        DU( 4)=-tmpd1*A3                                              
        DU( 5)= ALP1*F3*XY/R5                                           
        DU( 6)= ALP1*F3*X*D/R5                                          
        DU( 7)= DU(5)                                                   
        DU( 8)=-tmpd1*B3
        DU( 9)= ALP1*F3*Y*D/R5                                          
        DU(10)=-DU(6)                                                   
        DU(11)=-DU(9)                                                   
        DU(12)= tmpd1*C3
        tmpd5 = POT4/PI2
        DO 555 I=1,12                                                   
 555       U(I)=U(I)+tmpd5 * DU(I)   
        ENDIF                  
      RETURN                                                            
      END                    


                                           
      SUBROUTINE  UB0(X,Y,D,Z,POT1,POT2,POT3,POT4,U)                    
      IMPLICIT C_PREC (A-H,O-Z)                                         
      C_PREC tmpd1,tmpd2,tmpd3,tmpd4,tmpd5
      DIMENSION U(12),DU(12)                                            
C                                                                       
C********************************************************************   
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****   
C*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****   
C********************************************************************   
C                                                                       
C***** INPUT                                                            
C*****   X,Y,D,Z : STATION COORDINATES IN FAULT SYSTEM                  
C*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY        
C***** OUTPUT                                                           
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     
C                                                                       
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  
      COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,     
     *           UY,VY,WY,UZ,VZ,WZ                                      
      DATA F0,F1,F2,F3,F4,F5,F8,F9                                      
     *        /0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,8.D0,9.D0/                 
      DATA PI2/6.283185307179586D0/                                     
C-----                                                                  
      C=D+Z                                                             
      RD=R+D                                                            
      D12=F1/(R*RD*RD)                                                  
      D32=D12*(F2*R+D)/R2                                               
      D33=D12*(F3*R+D)/(R2*RD)                                          
      D53=D12*(F8*R2+F9*R*D+F3*D2)/(R2*R2*RD)                           
      D54=D12*(F5*R2+F4*R*D+D2)/R3*D12                                  
C-----                                                                  
      FI1= Y*(D12-X2*D33)                                               
      FI2= X*(D12-Y2*D33)                                               
      FI3= X/R3-FI2                                                     
      FI4=-XY*D32                                                       
      FI5= F1/(R*RD)-X2*D32                                             
      FJ1=-F3*XY*(D33-X2*D54)                                           
      FJ2= F1/R3-F3*D12+F3*X2*Y2*D54                                    
      FJ3= A3/R3-FJ2                                                    
      FJ4=-F3*XY/R5-FJ1                                                 
      FK1=-Y*(D32-X2*D53)                                               
      FK2=-X*(D32-Y2*D53)                                               
      FK3=-F3*X*D/R5-FK2                                                
C-----                                                                  
      DO 111  I=1,12                                                    
 111     U(I)=F0                                                           
C======================================                                 
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 
C======================================                                 
      IF(POT1.NE.F0) THEN           
         tmpd1=alp3*sd
         tmpd2=f3/r5
         DU( 1)=-X2*QR  -tmpd1*FI1                                     
         DU( 2)=-XY*QR  -tmpd1*FI2                                     
         DU( 3)=-C*X*QR -tmpd1*FI4                                     
         DU( 4)=-X*QR*(F1+A5) -tmpd1*FJ1
         DU( 5)=-Y*QR*A5      -tmpd1*FJ2
         DU( 6)=-C*QR*A5      -tmpd1*FK1
         DU( 7)=-tmpd2*X2*UY      -tmpd1*FJ2
         DU( 8)=-tmpd2*XY*UY-X*QR -tmpd1*FJ4
         DU( 9)=-tmpd2*C*X*UY     -tmpd1*FK2
         DU(10)=-tmpd2*X2*UZ  +tmpd1*FK1
         DU(11)=-tmpd2*XY*UZ  +tmpd1*FK2
         DU(12)= tmpd2*X*(-C*UZ +tmpd1*Y)                              
         tmpd5 = POT1/PI2
         DO 222 I=1,12                                                   
 222        U(I)=U(I)+tmpd5 * DU(I)                                        
      ENDIF                                                             
C===================================                                    
C=====  DIP-SLIP CONTRIBUTION  =====                                    
C===================================                                    
      IF(POT2.NE.F0) THEN                                               
         tmpd1=alp3*sdcd
         tmpd2=f3/r5
         tmpd3=tmpd2*vy
         tmpd4=tmpd2*vz
         
         DU( 1)=-X*P*QR +tmpd1*FI3
         DU( 2)=-Y*P*QR +tmpd1*FI1                                   
         DU( 3)=-C*P*QR +tmpd1*FI5                                   
         DU( 4)=-P*QR*A5 +tmpd1*FJ3                                  
         DU( 5)= Y*P*QRX +tmpd1*FJ1                                  
         DU( 6)= C*P*QRX +tmpd1*FK3                                  
         DU( 7)=-tmpd3*X      +tmpd1*FJ1                          
         DU( 8)=-tmpd3*Y-P*QR +tmpd1*FJ2                          
         DU( 9)=-tmpd3*C      +tmpd1*FK1                          
         DU(10)=-tmpd4*X -tmpd1*FK3                               
         DU(11)=-tmpd4*Y -tmpd1*FK1
         DU(12)=-tmpd4*C +tmpd1*A3/R3                             
         tmpd5 = POT2/PI2
         DO 333 I=1,12                                                   
 333        U(I)=U(I)+tmpd5 *DU(I)                                        
      ENDIF                                                             
C========================================                               
C=====  TENSILE-FAULT CONTRIBUTION  =====                               
C========================================                               
      IF(POT3.NE.F0) THEN                                               
         tmpd1=alp3*sdsd
         tmpd2=q*qr
         DU( 1)= X*tmpd2 -tmpd1*FI3                                   
         DU( 2)= Y*tmpd2 -tmpd1*FI1                                   
         DU( 3)= C*tmpd2 -tmpd1*FI5                                   
         DU( 4)= tmpd2*A5 -tmpd1*FJ3                                  
         DU( 5)=-Y*Q*QRX -tmpd1*FJ1                                  
         DU( 6)=-C*Q*QRX -tmpd1*FK3                                  
         DU( 7)= X*QR*WY     -tmpd1*FJ1                              
         DU( 8)= QR*(Y*WY+Q) -tmpd1*FJ2                              
         DU( 9)= C*QR*WY     -tmpd1*FK1                              
         DU(10)= X*QR*WZ +tmpd1*FK3                                  
         DU(11)= Y*QR*WZ +tmpd1*FK1                                  
         DU(12)= C*QR*WZ -tmpd1*A3/R3                                
         tmpd5 = POT3/PI2
         DO 444 I=1,12                                                   
 444        U(I)=U(I)+tmpd5*DU(I)                                        
      ENDIF                                                             
C=========================================                              
C=====  INFLATE SOURCE CONTRIBUTION  =====                              
C=========================================                              
      IF(POT4.NE.F0) THEN                                               
         tmpd1=alp3/r3
         tmpd2=alp3*f3/r5 
         DU( 1)= tmpd1*X
         DU( 2)= tmpd1*Y                                               
         DU( 3)= tmpd1*D                                               
         DU( 4)= tmpd1*A3                                              
         DU( 5)=-tmpd2*XY                                           
         DU( 6)=-tmpd2*X*D                                          
         DU( 7)= DU(5)                                                   
         DU( 8)= tmpd1*B3                                              
         DU( 9)=-tmpd2*Y*D                                          
         DU(10)=-DU(6)                                                   
         DU(11)=-DU(9)                                                   
         DU(12)=-tmpd1*C3               
         tmpd5 = POT4/PI2
         DO 555 I=1,12                                                   
 555        U(I)=U(I)+tmpd5 * DU(I)                                        
      ENDIF                                                             
      RETURN                                                            
      END                     

                                          
      SUBROUTINE  UC0(X,Y,D,Z,POT1,POT2,POT3,POT4,U)                    
      IMPLICIT C_PREC (A-H,O-Z)                                         
c
      C_PREC tmpd1,tmpd2,tmpd3,tmpd4
      DIMENSION U(12),DU(12)                                            
C                                                                       
C********************************************************************   
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****   
C*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****   
C********************************************************************   
C                                                                       
C***** INPUT                                                            
C*****   X,Y,D,Z : STATION COORDINATES IN FAULT SYSTEM                  
C*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY        
C***** OUTPUT                                                           
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     
C                                                                       
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  
      COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,
     &     dummy(6)
      DATA F0,F1,F2,F3,F5,F7,F10,F15                                    
     *        /0.D0,1.D0,2.D0,3.D0,5.D0,7.D0,10.D0,15.D0/               
      DATA PI2/6.283185307179586D0/                                     
C-----                                                                  
      C=D+Z                                                             
      Q2=Q*Q                                                            
      R7=R5*R2                                                          
      A7=F1-F7*X2/R2                                                    
      B5=F1-F5*Y2/R2                                                    
      B7=F1-F7*Y2/R2                                                    
      C5=F1-F5*D2/R2                                                    
      C7=F1-F7*D2/R2                                                    
      D7=F2-F7*Q2/R2                                                    
      QR5=F5*Q/R2                                                       
      QR7=F7*Q/R2                                                       
      DR5=F5*D/R2                                                       
C-----                                                                  
      DO 111  I=1,12                                                    
 111     U(I)=F0                                                           
C======================================                                 
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 
C======================================                                 
      IF(POT1.NE.F0) THEN     
         tmpd1=alp4*y
         tmpd2=alp5*c
         DU( 1)=-ALP4*A3/R3*CD  +tmpd2*QR*A5                            
         DU( 2)= F3*X/R5*( tmpd1*CD +tmpd2*(SD-Y*QR5) )                
         DU( 3)= F3*X/R5*(-tmpd1*SD +tmpd2*(CD+D*QR5) )                
         DU( 4)= ALP4*F3*X/R5*(F2+A5)*CD   -tmpd2*QRX*(F2+A7)           
         DU( 5)= F3/R5*( tmpd1*A5*CD +tmpd2*(A5*SD-Y*QR5*A7) )         
         DU( 6)= F3/R5*(-tmpd1*A5*SD +tmpd2*(A5*CD+D*QR5*A7) )         
         DU( 7)= DU(5)                                                   
         DU( 8)= F3*X/R5*( ALP4*B5*CD -tmpd2*F5/R2*(F2*Y*SD+Q*B7) )     
         DU( 9)= F3*X/R5*(-ALP4*B5*SD +tmpd2*F5/R2*(D*B7*SD-Y*C7*CD) )  
         DU(10)= F3/R5*   (-ALP4*D*A5*CD +tmpd2*(A5*CD+D*QR5*A7) )      
         DU(11)= F15*X/R7*( tmpd1*D*CD  +tmpd2*(D*B7*SD-Y*C7*CD) )     
         DU(12)= F15*X/R7*(-tmpd1*D*SD  +tmpd2*(F2*D*CD-Q*C7) )        
         tmpd4 = POT1/PI2
         DO 222 I=1,12                                                   
 222        U(I)=U(I)+tmpd4 *DU(I)                                        
      ENDIF                                                             
C===================================                                    
C=====  DIP-SLIP CONTRIBUTION  =====                                    
C===================================                                    
      IF(POT2.NE.F0) THEN         
         tmpd3=alp5*c
         tmpd1=tmpd3*f3/r5
         tmpd2=tmpd3*f5

         DU( 1)= ALP4*F3*X*T/R5          -tmpd3*P*QRX                   
         DU( 2)=-ALP4/R3*(C2D-F3*Y*T/R2) + tmpd1*(S-Y*P*QR5)       
         DU( 3)=-ALP4*A3/R3*SDCD         + tmpd1*(T+D*P*QR5)       
         DU( 4)= ALP4*F3*T/R5*A5              -tmpd2*P*QR/R2*A7      
         DU( 5)= F3*X/R5*(ALP4*(C2D-F5*Y*T/R2)-tmpd2/R2*(S-Y*P*QR7)) 
         DU( 6)= F3*X/R5*(ALP4*(F2+A5)*SDCD   -tmpd2/R2*(T+D*P*QR7)) 
         DU( 7)= DU(5)                                                   
         DU( 8)= F3/R5*(ALP4*(F2*Y*C2D+T*B5)                             
     *        +tmpd3*(S2D-F10*Y*S/R2-P*QR5*B7)) 
         DU( 9)= F3/R5*(ALP4*Y*A5*SDCD-tmpd3*((F3+A5)*C2D+Y*P*DR5*QR7)) 
         DU(10)= F3*X/R5*(-ALP4*(S2D-T*DR5) -tmpd2/R2*(T+D*P*QR7))   
         DU(11)= F3/R5*(-ALP4*(D*B5*C2D+Y*C5*S2D)                        
     *        -tmpd3*((F3+A5)*C2D+Y*P*DR5*QR7))
         DU(12)= F3/R5*(-ALP4*D*A5*SDCD-tmpd3*(S2D-F10*D*T/R2+P*QR5*C7))
         tmpd4 = POT2/PI2
         DO 333 I=1,12                                                   
 333        U(I)=U(I)+tmpd4*DU(I)                                        
      ENDIF                                                             
C========================================                               
C=====  TENSILE-FAULT CONTRIBUTION  =====                               
C========================================                               
      IF(POT3.NE.F0) THEN                   
         tmpd1=alp5*f3
         tmpd2=tmpd1/r5

         DU( 1)= F3*X/R5*(-ALP4*S +ALP5*(C*Q*QR5-Z))                     
         DU( 2)= ALP4/R3*(S2D-F3*Y*S/R2)+tmpd2*(C*(T-Y+Y*Q*QR5)-Y*Z)
         DU( 3)=-ALP4/R3*(F1-A3*SDSD)   -tmpd2*(C*(S-D+D*Q*QR5)-D*Z)
         DU( 4)=-ALP4*F3*S/R5*A5 +ALP5*(C*QR*QR5*A7-F3*Z/R5*A5)          
         DU( 5)= F3*X/R5*(-ALP4*(S2D-F5*Y*S/R2)                          
     *        -ALP5*F5/R2*(C*(T-Y+Y*Q*QR7)-Y*Z)) 
         DU( 6)= F3*X/R5*( ALP4*(F1-(F2+A5)*SDSD)                        
     *        +ALP5*F5/R2*(C*(S-D+D*Q*QR7)-D*Z)) 
         DU( 7)= DU(5)                                                   
         DU( 8)= F3/R5*(-ALP4*(F2*Y*S2D+S*B5)                            
     *        -ALP5*(C*(F2*SDSD+F10*Y*(T-Y)/R2-Q*QR5*B7)+Z*B5)) 
         DU( 9)= F3/R5*( ALP4*Y*(F1-A5*SDSD)                             
     *        +ALP5*(C*(F3+A5)*S2D-Y*DR5*(C*D7+Z)))             
         DU(10)= F3*X/R5*(-ALP4*(C2D+S*DR5)                              
     *        +ALP5*(F5*C/R2*(S-D+D*Q*QR7)-F1-Z*DR5))            
         DU(11)= F3/R5*( ALP4*(D*B5*S2D-Y*C5*C2D)                        
     *        +ALP5*(C*((F3+A5)*S2D-Y*DR5*D7)-Y*(F1+Z*DR5)))     
         DU(12)= F3/R5*(-ALP4*D*(F1-A5*SDSD)                             
     *        -ALP5*(C*(C2D+F10*D*(S-D)/R2-Q*QR5*C7)+Z*(F1+C5))) 
         tmpd4 = POT3/PI2
         DO 444 I=1,12                                                   
 444        U(I)=U(I)+ tmpd4*DU(I)                                        
      ENDIF                                                             
C=========================================                              
C=====  INFLATE SOURCE CONTRIBUTION  =====                              
C=========================================                              
      IF(POT4.NE.F0) THEN                                               
         tmpd1=alp4*f3
         tmpd2=tmpd1*d/r5

         DU( 1)= tmpd2*x
         DU( 2)= tmpd2*y
         DU( 3)= ALP4*C3/R3                                              
         DU( 4)= tmpd2*A5                                         
         DU( 5)=-ALP4*F15*XY*D/R7                                        
         DU( 6)=-tmpd1*X/R5*C5                                         
         DU( 7)= DU(5)                                                   
         DU( 8)= tmpd2*B5                                         
         DU( 9)=-tmpd1*Y/R5*C5                                         
         DU(10)= DU(6)                                                   
         DU(11)= DU(9)                                                   
         DU(12)= tmpd2*(F2+C5)    
         tmpd4 = POT4/PI2
         DO 555 I=1,12                                                   
 555        U(I)=U(I)+tmpd4  *DU(I)                                        
      ENDIF                                                             
      RETURN                                                            
      END                    

      
      SUBROUTINE  DC3D(ALPHA,X,Y,Z,DEPTH,DIP,                           
     *     AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3,                  
     *     UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)  
      IMPLICIT C_PREC (A-H,O-Z)                                         
      C_PREC   ALPHA,X,Y,Z,DEPTH,DIP,AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3, 
     *         UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ        
c     my addtions, th
      C_PREC tmpd1,tmpd2,tmpd3,tmpd4
c     
C                                                                       
C********************************************************************   
C*****                                                          *****   
C*****    DISPLACEMENT AND STRAIN AT DEPTH                      *****   
C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   
C*****              CODED BY  Y.OKADA ... SEP.1991              *****   
C*****              REVISED ... NOV.1991, APR.1992, MAY.1993,   *****   
C*****                          JUL.1993                        *****   
C********************************************************************   
C                                                                       
C***** INPUT                                                            
C*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)           
C*****   X,Y,Z : COORDINATE OF OBSERVING POINT                          
C*****   DEPTH : DEPTH OF REFERENCE POINT                               
C*****   DIP   : DIP-ANGLE (DEGREE)                                     
C*****   AL1,AL2   : FAULT LENGTH RANGE                                 
C*****   AW1,AW2   : FAULT WIDTH RANGE                                  
C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              
C                                                                       
C***** OUTPUT                                                           
C*****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF DISL)               
C*****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT=(UNIT OF DISL) /             
C*****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH,AL,AW) )
C*****   UXZ,UYZ,UZZ : Z-DERIVATIVE                                     
C*****   IRET        : RETURN CODE  ( =0....NORMAL,   =1....SINGULAR )  
C                                                                       
      COMMON /C0/DUMMY(5),SD,CD,dummy2(5)
      DIMENSION  XI(2),ET(2),KXI(2),KET(2)                              
      DIMENSION  U(12),DU(12),DUA(12),DUB(12),DUC(12)                   
      DATA  F0,EPS/ 0.D0, EPS_DC3D_PREC /                                       
C-----                                                                  
      IF(Z.GT.0.) then
         print *,' ** POSITIVE Z WAS GIVEN IN SUB-DC3D'
         print *,z
         stop
      endif
c      print *,ALPHA,X,Y,Z,DEPTH,DIP,                           
c     *     AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3
c      stop
      DO 111 I=1,12                                                     
         U  (I)=F0                                                       
c     added this line
         du(i)=f0
         DUA(I)=F0                                                       
         DUB(I)=F0                                                       
         DUC(I)=F0                                                       
 111  CONTINUE                                                          
      AALPHA=ALPHA                                                      
      DDIP=DIP                                                          
      CALL DCCON0(AALPHA,DDIP)                                          
C-----                                                                  
      ZZ=Z                                                              
      DD1=DISL1                                                         
      DD2=DISL2                                                         
      DD3=DISL3                                                         
      XI(1)=X-AL1                                                       
      XI(2)=X-AL2                                                       
      IF(DABS(XI(1)).LT.EPS) XI(1)=F0                                   
      IF(DABS(XI(2)).LT.EPS) XI(2)=F0                                   


C======================================                                 
C=====  REAL-SOURCE CONTRIBUTION  =====                                 
C======================================                                 
      D=DEPTH+Z                                                         
      P=Y*CD+D*SD                                                       
      Q=Y*SD-D*CD                                                       
      ET(1)=P-AW1                                                       
      ET(2)=P-AW2                                                       
      IF(DABS(Q).LT.EPS)  Q=F0                                          
      IF(DABS(ET(1)).LT.EPS) ET(1)=F0                                   
      IF(DABS(ET(2)).LT.EPS) ET(2)=F0                                   


C--------------------------------
C-----REJECT SINGULAR CASE -----                                       
C--------------------------------
C----- ON FAULT EDGE                                                    
      IF(Q.EQ.F0                                                        
     *     .AND.(    (XI(1)*XI(2).LE.F0 .AND. ET(1)*ET(2).EQ.F0)          
     *     .OR.(ET(1)*ET(2).LE.F0 .AND. XI(1)*XI(2).EQ.F0) ))       
     *     GO TO 99                                                       
C-----ON NEGATIVE EXTENSION OF FAULT EDGE                              
      KXI(1)=0                                                          
      KXI(2)=0                                                          
      KET(1)=0                                                          
      KET(2)=0                                                          
c     changed this, th
      tmpd1=XI(2)*XI(2)
      tmpd2=Q*Q
      tmpd3=ET(2)*ET(2)
      R12=DSQRT(XI(1)*XI(1)+tmpd3+tmpd2)                            
      R21=DSQRT(tmpd1+ET(1)*ET(1)+tmpd2)                            
      R22=DSQRT(tmpd1+tmpd3+tmpd2)                            
c     
      IF(XI(1).LT.F0 .AND. R21+XI(2).LT.EPS) KXI(1)=1                   
      IF(XI(1).LT.F0 .AND. R22+XI(2).LT.EPS) KXI(2)=1                   
      IF(ET(1).LT.F0 .AND. R12+ET(2).LT.EPS) KET(1)=1                   
      IF(ET(1).LT.F0 .AND. R22+ET(2).LT.EPS) KET(2)=1                   
C=====                                                                  
      DO 223 K=1,2                                                      
      DO 222 J=1,2                                                      
         CALL DCCON2(XI(J),ET(K),Q,SD,CD,KXI(K),KET(J))                  
         CALL UA(XI(J),ET(K),Q,DD1,DD2,DD3,DUA)                          
C-----
         DO 220 I=1,10,3                                                 
            DU(I)  =-DUA(I)                                               
            DU(I+1)=-DUA(I+1)*CD+DUA(I+2)*SD                              
            DU(I+2)=-DUA(I+1)*SD-DUA(I+2)*CD                              
            IF(I.LT.10) GO TO 220                                         
            DU(I)  =-DU(I)                                                
            DU(I+1)=-DU(I+1)                                              
            DU(I+2)=-DU(I+2)                                              
 220     CONTINUE                                                        
         DO 221 I=1,12                                                   
            IF(J+K.NE.3) U(I)=U(I)+DU(I)                                  
            IF(J+K.EQ.3) U(I)=U(I)-DU(I)                                  
 221     CONTINUE                                                        
C-----                                                                  
 222  CONTINUE                                                          
 223  CONTINUE                                                          
C=======================================                                
C=====  IMAGE-SOURCE CONTRIBUTION  =====                                
C=======================================
      D=DEPTH-Z                                                         
      P=Y*CD+D*SD                                                       
      Q=Y*SD-D*CD                                                       
      ET(1)=P-AW1                                                       
      ET(2)=P-AW2                                                       
      IF(DABS(Q).LT.EPS)  Q=F0                                          
      IF(DABS(ET(1)).LT.EPS) ET(1)=F0                                   
      IF(DABS(ET(2)).LT.EPS) ET(2)=F0                                   
C--------------------------------                                       
C----- REJECT SINGULAR CASE -----                                       
C--------------------------------                                       
C----- ON FAULT EDGE                                                    
      IF(Q.EQ.F0                                                        
     *     .AND.(    (XI(1)*XI(2).LE.F0 .AND. ET(1)*ET(2).EQ.F0)          
     *     .OR.(ET(1)*ET(2).LE.F0 .AND. XI(1)*XI(2).EQ.F0) ))       
     *     GO TO 99                                                       
C----- ON NEGATIVE EXTENSION OF FAULT EDGE                              
      KXI(1)=0                                                          
      KXI(2)=0                                                          
      KET(1)=0                                                          
      KET(2)=0       
c     changed this, th
      tmpd1=XI(2)*XI(2)
      tmpd2=Q*Q
      tmpd3=ET(2)*ET(2)
c                                                   
      R12=DSQRT(XI(1)*XI(1)+tmpd3+tmpd2)                            
      R21=DSQRT(tmpd1+ET(1)*ET(1)+tmpd2)                            
      R22=DSQRT(tmpd1+tmpd3+tmpd2)                            
      IF(XI(1).LT.F0 .AND. R21+XI(2).LT.EPS) KXI(1)=1                   
      IF(XI(1).LT.F0 .AND. R22+XI(2).LT.EPS) KXI(2)=1                   
      IF(ET(1).LT.F0 .AND. R12+ET(2).LT.EPS) KET(1)=1                   
      IF(ET(1).LT.F0 .AND. R22+ET(2).LT.EPS) KET(2)=1                   
C=====                                                                  
      DO 334 K=1,2                                                      
      DO 333 J=1,2                                                      
         CALL DCCON2(XI(J),ET(K),Q,SD,CD,KXI(K),KET(J))                  
         CALL UA(XI(J),ET(K),Q,DD1,DD2,DD3,DUA)                          
         CALL UB(XI(J),ET(K),Q,DD1,DD2,DD3,DUB)                          
         CALL UC(XI(J),ET(K),Q,ZZ,DD1,DD2,DD3,DUC)                       
C-----
         DO 330 I=1,10,3                                                 
            DU(I)=DUA(I)+DUB(I)+Z*DUC(I)
            tmpd1=z*duc(i+1)
            tmpd2=z*duc(i+2)
            tmpd3=dua(i+1)+dub(i+1)
            tmpd4=dua(i+2)+dub(i+2)
            DU(I+1)=(tmpd3+tmpd1)*CD                     
     *           -(tmpd4+tmpd2)*SD                     
            DU(I+2)=(tmpd3-tmpd1)*SD                     
     *           +(tmpd4-tmpd2)*CD                     
            IF(I.LT.10) GO TO 330                                         
            DU(10)=DU(10)+DUC(1)                                          
            DU(11)=DU(11)+DUC(2)*CD-DUC(3)*SD                             
            DU(12)=DU(12)-DUC(2)*SD-DUC(3)*CD                             
 330     CONTINUE                                                        
         DO 331 I=1,12                                                   
            IF(J+K.NE.3) U(I)=U(I)+DU(I)                                  
            IF(J+K.EQ.3) U(I)=U(I)-DU(I)                                  
 331     CONTINUE                                                        
C-----                                                                  
 333  CONTINUE                                                          
 334  CONTINUE                                                          
C=====
      UX=U(1)                                                           
      UY=U(2)                                                           
      UZ=U(3)                                                           
      UXX=U(4)                                                          
      UYX=U(5)                                                          
      UZX=U(6)                                                          
      UXY=U(7)                                                          
      UYY=U(8)                                                          
      UZY=U(9)                                                          
      UXZ=U(10)                                                         
      UYZ=U(11)                                                         
      UZZ=U(12)                                                         
c      print *,u(1),u(2),u(3)
c      stop
      IRET=0                                                            
      RETURN                                                            
C===========================================                            
C=====  IN CASE OF SINGULAR (ON EDGE)  =====                            
C===========================================                            
   99 UX=F0                                                             
      UY=F0                                                             
      UZ=F0                                                             
      UXX=F0                                                            
      UYX=F0                                                            
      UZX=F0                                                            
      UXY=F0                                                            
      UYY=F0                                                            
      UZY=F0                                                            
      UXZ=F0                                                            
      UYZ=F0                                                            
      UZZ=F0                                                            
      IRET=1                                                            
      RETURN                                                            
      END                                     
                          
      SUBROUTINE  UA(XI,ET,Q,DISL1,DISL2,DISL3,U)                       
      IMPLICIT C_PREC (A-H,O-Z)                                         
      C_PREC tmpd1,tmpd2,tmpd5
      DIMENSION U(12),DU(12)                                            
C                                                                       
C********************************************************************   
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-A)             *****   
C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   
C********************************************************************   
C                                                                       
C***** INPUT                                                            
C*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  
C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              
C***** OUTPUT                                                           
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     
C                                                                       
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  
      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,  
     *           EY,EZ,FY,FZ,GY,GZ,HY,HZ                                
      DATA F0,F2,PI2/0.D0,2.D0,6.283185307179586477D0/   
C-----                                                                  
      DO 111  I=1,12                                                    
 111     U(I)=F0                                                           
      XY=XI*Y11                                                         
      QX=Q *X11                                                         
      QY=Q *Y11                                                         
C======================================
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 
C======================================                                 
      IF(DISL1.NE.F0) THEN          
         tmpd1=alp2*q
         tmpd2=alp1*xy
         DU( 1)=    TT/F2 +ALP2*XI*QY                                    
         DU( 2)=           tmpd1/R                                      
         DU( 3)= ALP1*ALE -tmpd1*QY                                     
         DU( 4)=-ALP1*QY  -tmpd1*XI2*Y32                                
         DU( 5)=          -tmpd1*XI/R3                                  
         DU( 6)= tmpd2  +ALP2*XI*Q2*Y32                                
         DU( 7)= tmpd2*SD        +ALP2*XI*FY+D/F2*X11                  
         DU( 8)=                    ALP2*EY                              
         DU( 9)= ALP1*(CD/R+QY*SD) -tmpd1*FY                            
         DU(10)= tmpd2*CD        +ALP2*XI*FZ+Y/F2*X11                  
         DU(11)=                    ALP2*EZ                              
         DU(12)=-ALP1*(SD/R-QY*CD) -tmpd1*FZ                            
         tmpd5 = DISL1/PI2
         DO 222 I=1,12                                                   
 222        U(I)=U(I)+tmpd5*DU(I)                                       
      ENDIF                                                             
C======================================                                 
C=====    DIP-SLIP CONTRIBUTION   =====                                 
C======================================                                 
      IF(DISL2.NE.F0) THEN                                              
         tmpd1=alp2*q
         tmpd2=alp1*x11
         DU( 1)=           tmpd1/R                                      
         DU( 2)=    TT/F2 +ALP2*ET*QX                                    
         DU( 3)= ALP1*ALX -tmpd1*QX                                     
         DU( 4)=        -tmpd1*XI/R3                                    
         DU( 5)= -QY/F2 -tmpd1*ET/R3                                    
         DU( 6)= ALP1/R +ALP2*Q2/R3                                      
         DU( 7)=                      ALP2*EY                            
         DU( 8)= tmpd2*D+XY/F2*SD +ALP2*ET*GY                         
         DU( 9)= tmpd2*Y          -tmpd1*GY                          
         DU(10)=                      ALP2*EZ                            
         DU(11)= tmpd2*Y+XY/F2*CD +ALP2*ET*GZ                         
         DU(12)=-tmpd2*D          -tmpd1*GZ                          
         tmpd5 = DISL2/PI2
         DO 333 I=1,12                                                   
 333        U(I)=U(I)+tmpd5*DU(I)                                       
      ENDIF                                                             
C========================================                               
C=====  TENSILE-FAULT CONTRIBUTION  =====                               
C========================================                               
      IF(DISL3.NE.F0) THEN                                              
         tmpd1=alp2*q
         tmpd2=alp1*x11
         DU( 1)=-ALP1*ALE -tmpd1*QY                                     
         DU( 2)=-ALP1*ALX -tmpd1*QX                                     
         DU( 3)=    TT/F2 -ALP2*(ET*QX+XI*QY)                            
         DU( 4)=-ALP1*XY  +ALP2*XI*Q2*Y32                                
         DU( 5)=-ALP1/R   +ALP2*Q2/R3                                    
         DU( 6)=-ALP1*QY  -tmpd1*Q2*Y32                                 
         DU( 7)=-ALP1*(CD/R+QY*SD)  -tmpd1*FY                           
         DU( 8)=-tmpd2*Y         -tmpd1*GY                           
         DU( 9)= ALP1*(D*X11+XY*SD) +tmpd1*HY                           
         DU(10)= ALP1*(SD/R-QY*CD)  -tmpd1*FZ                           
         DU(11)= tmpd2*D         -tmpd1*GZ                           
         DU(12)= ALP1*(Y*X11+XY*CD) +tmpd1*HZ                           
         tmpd5 = DISL3/PI2
         DO 444 I=1,12                                                   
 444        U(I)=U(I)+tmpd5 *DU(I)                                       
      ENDIF                                                             
      RETURN                                                            
      END                   
                                            
      SUBROUTINE  UB(XI,ET,Q,DISL1,DISL2,DISL3,U)                       
      IMPLICIT C_PREC (A-H,O-Z)                                         
      DIMENSION U(12),DU(12)                                            
C                                                                       
C********************************************************************   
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****   
C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   
C********************************************************************   
C                                                                       
C***** INPUT                                                            
C*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  
C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              
C***** OUTPUT                                                           
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     
C                                                                       
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  
      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,  
     *           EY,EZ,FY,FZ,GY,GZ,HY,HZ                                
      DATA  F0,F1,F2,PI2/0.D0,1.D0,2.D0,6.283185307179586D0/            
      C_PREC tmpd1,tmpd2
C-----                                                                  
      RD=R+D                                                            
      D11=F1/(R*RD)                                                     
      AJ2=XI*Y/RD*D11                                                   
      AJ5=-(D+Y*Y/RD)*D11                                               
      IF(CD.NE.F0) THEN                                                 
         IF(XI.EQ.F0) THEN                                               
            AI4=F0                                                        
         ELSE                                                            
            X=DSQRT(XI2+Q2)                                               
            AI4=F1/CDCD*( XI/RD*SDCD                                      
     *           +F2*DATAN((ET*(X+Q*CD)+X*(R+X)*SD)/(XI*(R+X)*CD)) )        
         ENDIF                                                           
         AI3=(Y*CD/RD-ALE+SD*DLOG(RD))/CDCD                              
         AK1=XI*(D11-Y11*SD)/CD                                          
         AK3=(Q*Y11-Y*D11)/CD                                            
         AJ3=(AK1-AJ2*SD)/CD                                             
         AJ6=(AK3-AJ5*SD)/CD                                             
      ELSE                                                              
         RD2=RD*RD                                                       
         AI3=(ET/RD+Y*Q/RD2-ALE)/F2                                      
         AI4=XI*Y/RD2/F2                                                 
         AK1=XI*Q/RD*D11                                                 
         AK3=SD/RD*(XI2*D11-F1)                                          
         AJ3=-XI/RD2*(Q2*D11-F1/F2)                                      
         AJ6=-Y/RD2*(XI2*D11-F1/F2)                                      
      ENDIF                                                             
C-----                                                                  
      XY=XI*Y11                                                         
      AI1=-XI/RD*CD-AI4*SD                                              
      AI2= DLOG(RD)+AI3*SD                                              
      AK2= F1/R+AK3*SD                                                  
      AK4= XY*CD-AK1*SD                                                 
      AJ1= AJ5*CD-AJ6*SD                                                
      AJ4=-XY-AJ2*CD+AJ3*SD                                             
C=====                                                                  
      DO 111  I=1,12                                                    
 111     U(I)=F0                                                           
      QX=Q*X11                                                          
      QY=Q*Y11                                                          
C======================================                                 
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 
C======================================                                 
      IF(DISL1.NE.F0) THEN    
         tmpd2 = alp3 * sd
         DU( 1)=-XI*QY-TT -TMPD2*AI1                                   
         DU( 2)=-Q/R      +TMPD2*Y/RD                                  
         DU( 3)= Q*QY     -TMPD2*AI2                                   
         DU( 4)= XI2*Q*Y32 -TMPD2*AJ1                                  
         DU( 5)= XI*Q/R3   -TMPD2*AJ2                                  
         DU( 6)=-XI*Q2*Y32 -TMPD2*AJ3                                  
         DU( 7)=-XI*FY-D*X11 +TMPD2*(XY+AJ4)                           
         DU( 8)=-EY          +TMPD2*(F1/R+AJ5)                         
         DU( 9)= Q*FY        -TMPD2*(QY-AJ6)                           
         DU(10)=-XI*FZ-Y*X11 +TMPD2*AK1                                
         DU(11)=-EZ          +TMPD2*Y*D11                              
         DU(12)= Q*FZ        +TMPD2*AK2                 
         tmpd1 = DISL1/PI2
         DO 222 I=1,12                                                   
 222        U(I)=U(I)+ tmpd1 * DU(I)                                       
      ENDIF                                                             
C======================================                                 
C=====    DIP-SLIP CONTRIBUTION   =====                                 
C======================================                                 
      IF(DISL2.NE.F0) THEN                    
         tmpd2 = alp3 * sdcd
         DU( 1)=-Q/R      +TMPD2*AI3                                 
         DU( 2)=-ET*QX-TT -TMPD2*XI/RD                               
         DU( 3)= Q*QX     +TMPD2*AI4                                 
         DU( 4)= XI*Q/R3     +TMPD2*AJ4                              
         DU( 5)= ET*Q/R3+QY  +TMPD2*AJ5                              
         DU( 6)=-Q2/R3       +TMPD2*AJ6                              
         DU( 7)=-EY          +TMPD2*AJ1                              
         DU( 8)=-ET*GY-XY*SD +TMPD2*AJ2                              
         DU( 9)= Q*GY        +TMPD2*AJ3                              
         DU(10)=-EZ          -TMPD2*AK3                              
         DU(11)=-ET*GZ-XY*CD -TMPD2*XI*D11                           
         DU(12)= Q*GZ        -TMPD2*AK4                              
         tmpd1 = DISL2/PI2
         DO 333 I=1,12                                                   
 333        U(I)=U(I)+tmpd1*DU(I)                                       
      ENDIF                                                             
C========================================                               
C=====  TENSILE-FAULT CONTRIBUTION  =====                               
C========================================                               
      IF(DISL3.NE.F0) THEN                    
         tmpd2 = alp3 * sdsd
         DU( 1)= Q*QY           -TMPD2*AI3                           
         DU( 2)= Q*QX           +TMPD2*XI/RD                         
         DU( 3)= ET*QX+XI*QY-TT -TMPD2*AI4                           
         DU( 4)=-XI*Q2*Y32 -TMPD2*AJ4                                
         DU( 5)=-Q2/R3     -TMPD2*AJ5                                
         DU( 6)= Q*Q2*Y32  -TMPD2*AJ6                                
         DU( 7)= Q*FY -TMPD2*AJ1                                     
         DU( 8)= Q*GY -TMPD2*AJ2                                     
         DU( 9)=-Q*HY -TMPD2*AJ3                                     
         DU(10)= Q*FZ +TMPD2*AK3                                     
         DU(11)= Q*GZ +TMPD2*XI*D11
         DU(12)=-Q*HZ +TMPD2*AK4
         tmpd1 = DISL3/PI2
         DO 444 I=1,12                                                   
 444        U(I)=U(I)+tmpd1 * DU(I)                                       
      ENDIF                                                             
      RETURN                                                            
      END                     

                                          
      SUBROUTINE  UC(XI,ET,Q,Z,DISL1,DISL2,DISL3,U)                     
      IMPLICIT C_PREC (A-H,O-Z)                                         
      C_PREC tmpd1,tmpd2,tmpd5
      DIMENSION U(12),DU(12)                                            
C                                                                       
C********************************************************************   
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-C)             *****   
C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   
C********************************************************************   
C                                                                       
C***** INPUT                                                            
C*****   XI,ET,Q,Z   : STATION COORDINATES IN FAULT SYSTEM              
C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              
C***** OUTPUT                                                           
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     
C                                                                       
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  
      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,  
     *           EY,EZ,FY,FZ,GY,GZ,HY,HZ                                
      DATA F0,F1,F2,F3,PI2/0.D0,1.D0,2.D0,3.D0,6.283185307179586D0/     
C-----                                                                  
      C=D+Z                                                             
      X53=(8.D0*R2+9.D0*R*XI+F3*XI2)*X11*X11*X11/R2                     
      Y53=(8.D0*R2+9.D0*R*ET+F3*ET2)*Y11*Y11*Y11/R2                     
      H=Q*CD-Z                                                          
      Z32=SD/R3-H*Y32                                                   
      Z53=F3*SD/R5-H*Y53                                                
      Y0=Y11-XI2*Y32                                                    
      Z0=Z32-XI2*Z53                                                    
      PPY=CD/R3+Q*Y32*SD                                                
      PPZ=SD/R3-Q*Y32*CD                                                
      QQ=Z*Y32+Z32+Z0                                                   
      QQY=F3*C*D/R5-QQ*SD                                               
      QQZ=F3*C*Y/R5-QQ*CD+Q*Y32                                         
      XY=XI*Y11                                                         
      QX=Q*X11                                                          
      QY=Q*Y11                                                          
      QR=F3*Q/R5                                                        
      CQX=C*Q*X53                                                       
      CDR=(C+D)/R3                                                      
      YY0=Y/R3-Y0*CD                                                    
C=====                                                                  
      DO 111  I=1,12                                                    
  111 U(I)=F0                                                           
C======================================                                 
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 
C======================================                                 
      IF(DISL1.NE.F0) THEN          
         tmpd1=alp5*q
         tmpd2=alp4*xi
         DU( 1)= ALP4*XY*CD           -tmpd1*XI*Z32                     
         DU( 2)= ALP4*(CD/R+F2*QY*SD) -tmpd1*C/R3                       
         DU( 3)= ALP4*QY*CD           -ALP5*(C*ET/R3-Z*Y11+XI2*Z32)      
         DU( 4)= ALP4*Y0*CD                  -tmpd1*Z0                  
         DU( 5)=-tmpd2*(CD/R3+F2*Q*Y32*SD) +ALP5*C*XI*QR               
         DU( 6)=-tmpd2*Q*Y32*CD            +ALP5*XI*(F3*C*ET/R5-QQ)    
         DU( 7)=-tmpd2*PPY*CD    -ALP5*XI*QQY                          
         DU( 8)= ALP4*F2*(D/R3-Y0*SD)*SD-Y/R3*CD                         
     *        -ALP5*(CDR*SD-ET/R3-C*Y*QR)           
         DU( 9)=-ALP4*Q/R3+YY0*SD  +ALP5*(CDR*CD+C*D*QR-(Y0*CD+Q*Z0)*SD) 
         DU(10)= tmpd2*PPZ*CD    -ALP5*XI*QQZ                          
         DU(11)= ALP4*F2*(Y/R3-Y0*CD)*SD+D/R3*CD -ALP5*(CDR*CD+C*D*QR)   
         DU(12)=         YY0*CD    -ALP5*(CDR*SD-C*Y*QR-Y0*SDSD+Q*Z0*CD) 
         tmpd5 = DISL1/PI2
         DO 222 I=1,12                                                   
 222        U(I)=U(I)+tmpd5*DU(I)                                       
      ENDIF                                                             
C======================================                                 
C=====    DIP-SLIP CONTRIBUTION   =====                                 
C======================================                                 
      IF(DISL2.NE.F0) THEN                                              
         tmpd1=alp5*c
         DU( 1)= ALP4*CD/R -QY*SD -tmpd1*Q/R3                           
         DU( 2)= ALP4*Y*X11       -tmpd1*ET*Q*X32                       
         DU( 3)=     -D*X11-XY*SD -tmpd1*(X11-Q2*X32)                   
         DU( 4)=-ALP4*XI/R3*CD +tmpd1*XI*QR +XI*Q*Y32*SD                
         DU( 5)=-ALP4*Y/R3     +tmpd1*ET*QR                             
         DU( 6)=    D/R3-Y0*SD +tmpd1/R3*(F1-F3*Q2/R2)                  
         DU( 7)=-ALP4*ET/R3+Y0*SDSD -ALP5*(CDR*SD-C*Y*QR)                
         DU( 8)= ALP4*(X11-Y*Y*X32) -tmpd1*((D+F2*Q*CD)*X32-Y*ET*Q*X53) 
         DU( 9)=  XI*PPY*SD+Y*D*X32 +tmpd1*((Y+F2*Q*SD)*X32-Y*Q2*X53)   
         DU(10)=      -Q/R3+Y0*SDCD -ALP5*(CDR*CD+C*D*QR)                
         DU(11)= ALP4*Y*D*X32       -tmpd1*((Y-F2*Q*SD)*X32+D*ET*Q*X53) 
         DU(12)=-XI*PPZ*SD+X11-D*D*X32-tmpd1*((D-F2*Q*CD)*X32-D*Q2*X53) 
         tmpd5 = DISL2/PI2
         DO 333 I=1,12                                                   
 333        U(I)=U(I)+tmpd5*DU(I)                                       
      ENDIF                                                             
C========================================                               
C=====  TENSILE-FAULT CONTRIBUTION  =====                               
C========================================                               
      IF(DISL3.NE.F0) THEN                                              
         tmpd1=alp5*c
         tmpd2=alp4*xi
         DU( 1)=-ALP4*(SD/R+QY*CD)   -ALP5*(Z*Y11-Q2*Z32)                
         DU( 2)= ALP4*F2*XY*SD+D*X11 -tmpd1*(X11-Q2*X32)                
         DU( 3)= ALP4*(Y*X11+XY*CD)  +ALP5*Q*(C*ET*X32+XI*Z32)           
         DU( 4)= tmpd2/R3*SD+XI*Q*Y32*CD+ALP5*XI*(F3*C*ET/R5-F2*Z32-Z0)
         DU( 5)= ALP4*F2*Y0*SD-D/R3 +tmpd1/R3*(F1-F3*Q2/R2)             
         DU( 6)=-ALP4*YY0           -ALP5*(C*ET*QR-Q*Z0)                 
         DU( 7)= ALP4*(Q/R3+Y0*SDCD)   +ALP5*(Z/R3*CD+C*D*QR-Q*Z0*SD)    
         DU( 8)=-tmpd2*F2*PPY*SD-Y*D*X32                               
     *        +tmpd1*((Y+F2*Q*SD)*X32-Y*Q2*X53)            
         DU( 9)=-ALP4*(XI*PPY*CD-X11+Y*Y*X32)                            
     *        +ALP5*(C*((D+F2*Q*CD)*X32-Y*ET*Q*X53)+XI*QQY) 
         DU(10)=  -ET/R3+Y0*CDCD -ALP5*(Z/R3*SD-C*Y*QR-Y0*SDSD+Q*Z0*CD)  
         DU(11)= tmpd2*F2*PPZ*SD-X11+D*D*X32                           
     *        -tmpd1*((D-F2*Q*CD)*X32-D*Q2*X53)            
         DU(12)= ALP4*(XI*PPZ*CD+Y*D*X32)                                
     *        +ALP5*(C*((Y-F2*Q*SD)*X32+D*ET*Q*X53)+XI*QQZ) 
         tmpd5 = DISL3/PI2
         DO 444 I=1,12                                                   
 444        U(I)=U(I)+tmpd5*DU(I)                                       
      ENDIF                                                             
      RETURN                                                            
      END                  
                                             
      SUBROUTINE  DCCON0(ALPHA,DIP)                                     
      IMPLICIT C_PREC (A-H,O-Z)  
C                                                                       
C*******************************************************************    
C*****   CALCULATE MEDIUM CONSTANTS AND FAULT-DIP CONSTANTS    *****    
C*******************************************************************    
C                                                                       
C***** INPUT                                                            
C*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)           
C*****   DIP   : DIP-ANGLE (DEGREE)                                     
C### CAUTION ### IF COS(DIP) IS SUFFICIENTLY SMALL, IT IS SET TO ZERO   
C                                                                       
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  
      DATA F0,F1,F2/0.0D0,1.0D0,2.0D0/             
      DATA EPS/EPS_DC3D_PREC/                                                
C-----                                                                  
      ALP1=(F1-ALPHA)/F2                                                
      ALP2= ALPHA/F2                                                    
      ALP3=(F1-ALPHA)/ALPHA                                             
      ALP4= F1-ALPHA                                                    
      ALP5= ALPHA                                                       
C-----                          
c     obtain sin and cos of dip, given in degrees
      call my_sincos_deg_ftn(sd,cd,dip)
      IF(DABS(CD).LT.EPS) THEN                                          
         CD=F0                                                           
         IF(SD.GT.F0) SD= F1                                             
         IF(SD.LT.F0) SD=-F1                                             
      ENDIF                                                             
      SDSD=SD*SD                                                        
      CDCD=CD*CD                                                        
      SDCD=SD*CD                                                        
      S2D=F2*SDCD                                                       
      C2D=CDCD-SDSD                                                     
      RETURN                                                            
      END        
                                                       
      SUBROUTINE  DCCON1(X,Y,D)                                         
      IMPLICIT C_PREC (A-H,O-Z)                                         
      C_PREC tmpd1
C                                                                       
C********************************************************************** 
C*****   CALCULATE STATION GEOMETRY CONSTANTS FOR POINT SOURCE    ***** 
C********************************************************************** 
C                                                                       
C***** INPUT                                                            
C*****   X,Y,D : STATION COORDINATES IN FAULT SYSTEM                    
C### CAUTION ### IF X,Y,D ARE SUFFICIENTLY SMALL, THEY ARE SET TO ZERO  
C                                                                       
      COMMON /C0/DUMMY(5),SD,CD,dummy2(5)
      COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,     
     *           UY,VY,WY,UZ,VZ,WZ                                      
      DATA  F0,F1,F3,F5,EPS/0.D0,1.D0,3.D0,5.D0,EPS_DC3D_PREC/               
C-----                                                                  
      IF(DABS(X).LT.EPS) X=F0                                           
      IF(DABS(Y).LT.EPS) Y=F0                                           
      IF(DABS(D).LT.EPS) D=F0                                           
      P=Y*CD+D*SD                                                       
      Q=Y*SD-D*CD                                                       
      S=P*SD+Q*CD                                                       
      T=P*CD-Q*SD                                                       
      XY=X*Y                                                            
      X2=X*X                                                            
      Y2=Y*Y                                                            
      D2=D*D                                                            
      R2=X2+Y2+D2                                                       
      R =DSQRT(R2)                                                      
      IF(R.EQ.F0) RETURN                                                
      R3=R *R2                                                          
      R5=R3*R2                                                          
      R7=R5*R2                                                          
C-----                                                                  
      A3=F1-F3*X2/R2                                                    
      A5=F1-F5*X2/R2                                                    
      B3=F1-F3*Y2/R2                                                    
      C3=F1-F3*D2/R2                                                    
C-----                                                                  
      QR=F3*Q/R5                                                        
      QRX=F5*QR*X/R2                                                    
C-----                              
      tmpd1=q/r2*f5
      UY=SD-Y*tmpd1                                                   
      UZ=CD+D*tmpd1
      VY=S -Y*P*tmpd1                                                 
      VZ=T +D*P*tmpd1                                                 
      WY=UY+SD                                                          
      WZ=UZ+CD                                                          
      RETURN                                                            
      END                                                               

      SUBROUTINE  DCCON2(XI,ET,Q,SD,CD,KXI,KET)                         
      IMPLICIT C_PREC (A-H,O-Z)                                         
      C_PREC tmpd1
C                                                                       
C********************************************************************** 
C*****   CALCULATE STATION GEOMETRY CONSTANTS FOR FINITE SOURCE   ***** 
C********************************************************************** 
C                                                                       
C***** INPUT                                                            
C*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  
C*****   SD,CD   : SIN, COS OF DIP-ANGLE                                
C*****   KXI,KET : KXI=1, KET=1 MEANS R+XI<EPS, R+ET<EPS, RESPECTIVELY  
C                                                                       
C### CAUTION ### IF XI,ET,Q ARE SUFFICIENTLY SMALL, THEY ARE SET TO ZER0
C                                                                       
      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,  
     *           EY,EZ,FY,FZ,GY,GZ,HY,HZ                                
      DATA  F0,F1,F2,EPS/0.D0,1.D0,2.D0,EPS_DC3D_PREC/                       
C-----                                                                  
      IF(DABS(XI).LT.EPS) XI=F0                                         
      IF(DABS(ET).LT.EPS) ET=F0                                         
      IF(DABS( Q).LT.EPS)  Q=F0                                         
      XI2=XI*XI                                                         
      ET2=ET*ET                                                         
      Q2=Q*Q                                                            
      R2=XI2+ET2+Q2                                                     
      R =DSQRT(R2)                                                      
      IF(R.EQ.F0) RETURN                                                
      R3=R *R2                                                          
      R5=R3*R2                                                          
      Y =ET*CD+Q*SD                                                     
      D =ET*SD-Q*CD                                                     
C-----                                                                  
      IF(Q.EQ.F0) THEN                                                  
        TT=F0                                                           
      ELSE                                                              
        TT=DATAN(XI*ET/(Q*R))                                           
      ENDIF                                                             
C-----                                                                  
      IF(KXI.EQ.1) THEN                                                 
         ALX=-DLOG(R-XI)                                                 
         X11=F0                                                          
         X32=F0                                                          
      ELSE                                                              
         RXI=R+XI                                                        
         ALX=DLOG(RXI)                                                   
         X11=F1/(R*RXI)                                                  
         X32=(R+RXI)*X11*X11/R                                           
      ENDIF                                                             
C-----                                                                  
      IF(KET.EQ.1) THEN                                                 
         ALE=-DLOG(R-ET)                                                 
         Y11=F0                                                          
         Y32=F0                                                          
      ELSE                                                              
         RET=R+ET                                                        
         ALE=DLOG(RET)                                                   
         Y11=F1/(R*RET)                                                  
         Y32=(R+RET)*Y11*Y11/R                                           
      ENDIF                                                             
C-----                                                                  
      tmpd1=q/r3
      EY=SD/R-Y*tmpd1                     
      EZ=CD/R+D*tmpd1            
      FY=D/R3+XI2*Y32*SD                                                
      FZ=Y/R3+XI2*Y32*CD                                                
      GY=F2*X11*SD-Y*Q*X32                                              
      GZ=F2*X11*CD+D*Q*X32                                              
      HY=D*Q*X32+XI*Q*Y32*SD                                            
      HZ=Y*Q*X32+XI*Q*Y32*CD                                            
      RETURN                                                            
      END                                                               
