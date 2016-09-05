#ifdef USE_DOUBLE_PRECISION
#include "precision_double.h"
#else
#include "precision_single.h"
#endif
      SUBROUTINE  DCD3(ALPHA,X,Y,Z,DEPTH,DIP,
     & AL1,AL2,AW1,AW2,
     & UX1,UX2,UX3,UY1,UY2,UY3,UZ1,UZ2,UZ3,
     &    UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)   
      IMPLICIT REAL*8 (A-H,O-Z)                                          

      REAL*8 ALPHA,DEPTH,DIP,
     &    UX1,UX2,UX3,UY1,UY2,UY3,UZ1,UZ2,UZ3
	 
      DIMENSION UXX(3),UYX(3),UZX(3),UXY(3),UYY(3),UZY(3)
      DIMENSION UXZ(3),UYZ(3),UZZ(3)

C********************************************************************   
C*****                                                          *****   
C*****    DISPLACEMENT AND STRAIN AT DEPTH                      *****   
C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   
C*****                         CODED BY  Y.OKADA ... SEP 1991   *****   
C*****                         REVISED   Y.OKADA ... NOV 1991   *****   
C*****                         REVISED   Y.OKADA ... APR 1992   *****   
C*****                                                          *****   
C********************************************************************   
C                                                                       
C***** INPUT                                                            
C*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)           
C*****   X,Y,Z : COORDINATE OF OBSERVING POINT                          
C*****   DEPTH : DEPTH OF REFERENCE POINT                               
C*****   DIP   : DIP-ANGLE (DEGREE)                                     
C*****   AL1,AL2   : FAULT LENGTH RANGE                                 
C*****   AW1,AW2   : FAULT WIDTH RANGE                                  

C***** OUTPUT                                                           
C*****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF DISL)               
C*****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT=(UNIT OF DISL) /             
C*****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH,AL,AW) )
C*****   UXZ,UYZ,UZZ : Z-DERIVATIVE                                     
C*****   IRET        : RETURN CODE  ( =0....NORMAL,   =1....SINGULAR )  
C                                                                       
      COMMON /C0/DUMMY(5),SD,CD,DUMMY1(5)
      COMMON /C2/XI2,ET2,Q2,R,dum2(20)
      DIMENSION  U(12,3),DU(12,3),DUA(12,3),DUB(12,3),DUC(12,3)
      DATA  F0/0.D0/                                                    
      DATA EPS/1.D-6/                                                   
C-----                                                                  
      IF(Z.GT.0.) THEN
	  		WRITE(6,'('' ** POSITIVE Z WAS GIVEN IN SUB-DC3D'')') 
			IRET=1
	        RETURN
	  ENDIF
c     this was apparently added to use three components of slip
c     at the same time 
      DO 111 L=1,3
        DO 111 I=1,12                                                     
        U  (I,L)=F0                                                        
        DU (I,L)=F0                                                        
        DUA(I,L)=F0                                                       
        DUB(I,L)=F0                                                       
        DUC(I,L)=F0                                                       
  111 CONTINUE                                                         
      AALPHA=ALPHA                                                      
      DDIP=DIP                                                          
      CALL DCCON0(AALPHA,DDIP)                                          
C======================================                                 
C=====  REAL-SOURCE CONTRIBUTION  =====                                 
C======================================                                 
      D=DEPTH+Z                                                         
      P=Y*CD+D*SD                                                       
      Q=Y*SD-D*CD                                                       
      JXI=0                                                             
      JET=0                                                             
      IF((X-AL1)*(X-AL2).LE.0.) JXI=1                                   
      IF((P-AW1)*(P-AW2).LE.0.) JET=1

C-----                                                                  
      DO 223 K=1,2                                                      
      IF(K.EQ.1) ET=P-AW1                                               
      IF(K.EQ.2) ET=P-AW2                                               
      
      DO 222 J=1,2                                                      
        IF(J.EQ.1) XI=X-AL1                                           
        IF(J.EQ.2) XI=X-AL2                                             
C --   CALCULATE STATION GEOMETRY CONSTANTS FOR FINITE SOURCE  
        CALL DCCON2(XI,ET,Q,SD,CD)    
        IF(JXI.EQ.1 .AND. Q.EQ.F0 .AND. ET.EQ.F0) GO TO 99              
        IF(JET.EQ.1 .AND. Q.EQ.F0 .AND. XI.EQ.F0) GO TO 99              
        CALL UA(XI,ET,Q,DUA)                     
	
C-----   
	    DO 224 L=1,3
         DO 220 I=1,10,3                                                 
          DU(I,L)  =-DUA(I,L)                                               
          DU(I+1,L)=-DUA(I+1,L)*CD+DUA(I+2,L)*SD                              
          DU(I+2,L)=-DUA(I+1,L)*SD-DUA(I+2,L)*CD                              
          IF(I.LT.10) GO TO 220                                         
          DU(I,L)  =-DU(I,L)                                                
          DU(I+1,L)=-DU(I+1,L)                                              
          DU(I+2,L)=-DU(I+2,L)                                              
  220   CONTINUE                                                        
  
        DO 221 I=1,12                                                   
          IF(J+K.NE.3) U(I,L)=U(I,L)+DU(I,L)                                  
          IF(J+K.EQ.3) U(I,L)=U(I,L)-DU(I,L)                                  
  221   CONTINUE                                                        
  224  CONTINUE
C-----                                                                   
  222 CONTINUE                                                           
  223 CONTINUE                                                           
  
  
  
C=======================================                                
C=====  IMAGE-SOURCE CONTRIBUTION  =====                                
C=======================================                                
      ZZ=Z                                                              
      D=DEPTH-Z                                                         
      P=Y*CD+D*SD                                                       
      Q=Y*SD-D*CD                                                       
      JET=0                                                             
      IF((P-AW1)*(P-AW2).LE.0.) JET=1                                   
C-----                                                                  
      DO 334 K=1,2                                                      
      IF(K.EQ.1) ET=P-AW1                                               
      IF(K.EQ.2) ET=P-AW2                                               
      DO 333 J=1,2                                                      
        IF(J.EQ.1) XI=X-AL1                                             
        IF(J.EQ.2) XI=X-AL2                                             
        CALL DCCON2(XI,ET,Q,SD,CD)                                      
        CALL UA(XI,ET,Q,DUA)                                
        CALL UB(XI,ET,Q,DUB)                                
        CALL UC(XI,ET,Q,ZZ,DUC)                             
	
C-----                                                                  
	    DO 335 L=1,3
         DO 330 I=1,10,3                                                 
          DU(I,L)=DUA(I,L)+DUB(I,L)+Z*DUC(I,L)                                  
          DU(I+1,L)=(DUA(I+1,L)+DUB(I+1,L)+Z*DUC(I+1,L))*CD                     
     &           -(DUA(I+2,L)+DUB(I+2,L)+Z*DUC(I+2,L))*SD                     
          DU(I+2,L)=(DUA(I+1,L)+DUB(I+1,L)-Z*DUC(I+1,L))*SD                     
     &           +(DUA(I+2,L)+DUB(I+2,L)-Z*DUC(I+2,L))*CD                     
          IF(I.LT.10) GO TO 330                                         
          DU(10,L)=DU(10,L)+DUC(1,L)                                          
          DU(11,L)=DU(11,L)+DUC(2,L)*CD-DUC(3,L)*SD                             
          DU(12,L)=DU(12,L)-DUC(2,L)*SD-DUC(3,L)*CD                             
  330    CONTINUE                                                        
         DO 331 I=1,12                                                   
          IF(J+K.NE.3) U(I,L)=U(I,L)+DU(I,L)                                  
          IF(J+K.EQ.3) U(I,L)=U(I,L)-DU(I,L)                                  
  331    CONTINUE 
  335  CONTINUE
C-----                                                                  
  333 CONTINUE                                                          
  334 CONTINUE                                                          
  
  
C=====                                                                  
      DO 440 L=1,3
      UXX(L)	=U(4,L)                                                          
      UYX(L)	=U(5,L)                                                          
      UZX(L)	=U(6,L)                                                          
      UXY(L)	=U(7,L)                                                          
      UYY(L)	=U(8,L)                                                          
      UZY(L)	=U(9,L)                                                          
      UXZ(L)	=U(10,L)                                                         
      UYZ(L)	=U(11,L)                                                         
440   UZZ(L)	=U(12,L)                                                         

      UX1	=U(1,1)                                                           
      UY1	=U(2,1)                                                           
      UZ1	=U(3,1)                                                           
      UX2	=U(1,2)                                                           
      UY2	=U(2,2)                                                           
      UZ2	=U(3,2)                                                           
      UX3	=U(1,3)                                                           
      UY3	=U(2,3)                                                           
      UZ3	=U(3,3)                                                           
      IRET=0                                                            
      RETURN                                                      

C=======================================                                
C=====  IN CASE OF SINGULAR (R=0)  =====                                
C=======================================                                
99    UX1	=F0                                                           
      UY1	=F0                                                           
      UZ1	=F0                                                           
      UX2	=F0                                                          
      UY2	=F0                                                           
      UZ2	=F0                                                           
      UX3	=F0                                                           
      UY3	=F0                                                           
      UZ3	=F0                                                           
      DO 441 L=1,3
      UXX(L)	=F0                                                            
      UYX(L)	=F0                                                            
      UZX(L)	=F0                                                            
      UXY(L)	=F0                                                            
      UYY(L)	=F0                                                            
      UZY(L)	=F0                                                            
      UXZ(L)	=F0                                                            
      UYZ(L)	=F0                                                            
441   UZZ(L)	=F0                                                            
      write(*,'(a)') ' DEPTH,X,Y,Z,SD,CD,P,Q,ET,XI'
      write(*,*) DEPTH,X,Y,Z,SD,CD,P,Q,ET,XI
      IRET=1                                                            
      RETURN                                                            
      END                                                               
      
      SUBROUTINE  UA(XI,ET,Q,U)                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION U(12,3)                                            
C                                                                       
C********************************************************************   
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-A)             *****   
C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   
C********************************************************************   
C                                                                       
C***** INPUT                                                            
C*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  

C***** OUTPUT                                                           
C*****   U(12,3) : DISPLACEMENT AND THEIR DERIVATIVES                     
C                                                                       
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,
     &           SDSD,CDCD,SDCD,S2D,C2D  
      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,  
     &           X11,Y11,X32,Y32,EY,EZ,FY,FZ,GY,GZ,HY,HZ                                
      DATA F0,F2,PI2/0.D0,2.D0,6.283185307179586D0/                     
C-----                                                                  
      DO 111 L=1,3
       DO 111  I=1,12                                                    
  111 U(I,L)=F0                                                           
      XY=XI*Y11                                                     
      QX=Q *X11                                                         
      QY=Q *Y11                                                         
C======================================                                 
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 
C======================================                                 
      U(1,1)=    (TT/F2 +ALP2*XI*QY)/PI2                                    
      U(2,1)=           (ALP2*Q/R)/PI2                                      
      U(3,1)= (ALP1*ALE -ALP2*Q*QY)/PI2                                     
      U(4,1)= (-ALP1*QY  -ALP2*XI2*Q*Y32)/PI2                                
      U(5,1)=           (-ALP2*XI*Q/R3)/PI2                                  
      U(6,1)= (ALP1*XY  +ALP2*XI*Q2*Y32)/PI2                               
      U(7,1)= (ALP1*XY*SD        
     1          	+ALP2*XI*FY+D/F2*X11)/PI2                  
      U(8,1)=           (ALP2*EY)/PI2                              
      U(9,1)= (ALP1*(CD/R+QY*SD) 
     1		        -ALP2*Q*FY)/PI2                            
      U(10,1)=(ALP1*XY*CD        
     1			+ALP2*XI*FZ+Y/F2*X11)/PI2                  
      U(11,1)=          (ALP2*EZ)/PI2                              
      U(12,1)=(-ALP1*(SD/R-QY*CD) 
     1			-ALP2*Q*FZ)/PI2                            
C======================================                                 
C=====    DIP-SLIP CONTRIBUTION   =====                                 
C======================================                                 
      U(1,2)=           (ALP2*Q/R)/PI2                                      
      U(2,2)=    (TT/F2 +ALP2*ET*QX)/PI2                                    
      U(3,2)=  (ALP1*ALX -ALP2*Q*QX)/PI2                                     
      U(4,2)=           (-ALP2*XI*Q/R3)/PI2                                    
      U(5,2)=  (-QY/F2   -ALP2*ET*Q/R3)/PI2                                    
      U(6,2)=  (ALP1/R   +ALP2*Q2/R3)/PI2                                      
      U(7,2)=            (ALP2*EY)/PI2                            
      U(8,2)=  (ALP1*D*X11+XY/F2*SD 
     1		        +ALP2*ET*GY)/PI2                         
      U(9,2)=  (ALP1*Y*X11         
     1                  -ALP2*Q*GY)/PI2                          
      U(10,2)=          (ALP2*EZ)/PI2                            
      U(11,2)= (ALP1*Y*X11+XY/F2*CD 
     1		        +ALP2*ET*GZ)/PI2                         
      U(12,2)=(-ALP1*D*X11      
     1		        -ALP2*Q*GZ)/PI2                          
C========================================                               
C=====  TENSILE-FAULT CONTRIBUTION  =====                               
C========================================                               
      U(1,3)=  (-ALP1*ALE -ALP2*Q*QY)/PI2                                     
      U(2,3)=  (-ALP1*ALX -ALP2*Q*QX)/PI2                                     
      U(3,3)=    (TT/F2 -ALP2*(ET*QX+XI*QY))/PI2                            
      U(4,3)=  (-ALP1*XY  +ALP2*XI*Q2*Y32)/PI2                                
      U(5,3)=  (-ALP1/R   +ALP2*Q2/R3)/PI2                                    
      U(6,3)=  (-ALP1*QY  -ALP2*Q*Q2*Y32)/PI2                                 
      U(7,3)=  (-ALP1*(CD/R+QY*SD)  -ALP2*Q*FY)/PI2                           
      U(8,3)=  (-ALP1*Y*X11         -ALP2*Q*GY)/PI2                           
      U(9,3)=  (ALP1*(D*X11+XY*SD) +ALP2*Q*HY)/PI2                           
      U(10,3)=(ALP1*(SD/R-QY*CD)  -ALP2*Q*FZ)/PI2                           
      U(11,3)=(ALP1*D*X11         -ALP2*Q*GZ)/PI2                           
      U(12,3)=(ALP1*(Y*X11+XY*CD) +ALP2*Q*HZ)/PI2                           
      RETURN                                                            
      END                                                               

      SUBROUTINE  UB(XI,ET,Q,U)                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION U(12,3)                                           
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
C*****   U(12,3) : DISPLACEMENT AND THEIR DERIVATIVES                     
C                                                                       
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,
     &                   SDSD,CDCD,SDCD,S2D,C2D  
      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,
     &           X11,Y11,X32,Y32,EY,EZ,FY,FZ,GY,GZ,HY,HZ                                
      DATA  F0,F1,F2,PI2/0.D0,1.D0,2.D0,6.283185307179586D0/            
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
     *     +F2*DATAN((ET*(X+Q*CD)+X*(R+X)*SD)/(XI*(R+X)*CD)) )        
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
      DO 111 L=1,3 
      DO 111  I=1,12                                                    
  111 U(I,L)=F0                                                          
      QX=Q*X11                                                          
      QY=Q*Y11                                                          
      
C======================================                                 
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 
C======================================                                 
      U(1,1)= (-XI*QY-TT -ALP3*AI1*SD)/PI2                                 
      U(2,1)= (-Q/R      +ALP3*Y/RD*SD)/PI2                                  
      U(3,1)= ( Q*QY     -ALP3*AI2*SD)/PI2                                   
      U(4,1)= ( XI2*Q*Y32 -ALP3*AJ1*SD)/PI2                                  
      U(5,1)= ( XI*Q/R3   -ALP3*AJ2*SD)/PI2                                  
      U(6,1)= (-XI*Q2*Y32 -ALP3*AJ3*SD)/PI2                                  
      U(7,1)= (-XI*FY-D*X11 +ALP3*(XY+AJ4)*SD)/PI2                          
      U(8,1)= (-EY          +ALP3*(F1/R+AJ5)*SD)/PI2                         
      U(9,1)= ( Q*FY        -ALP3*(QY-AJ6)*SD)/PI2                           
      U(10,1)= (-XI*FZ-Y*X11 +ALP3*AK1*SD)/PI2                               
      U(11,1)= (-EZ          +ALP3*Y*D11*SD)/PI2                              
      U(12,1)= ( Q*FZ        +ALP3*AK2*SD)/PI2                                

C======================================                                 
C=====    DIP-SLIP CONTRIBUTION   =====                                 
C======================================                                 
      U(1,2)=(-Q/R      +ALP3*AI3*SDCD)/PI2                                 
      U(2,2)=(-ET*QX-TT -ALP3*XI/RD*SDCD)/PI2                               
      U(3,2)=( Q*QX     +ALP3*AI4*SDCD)/PI2                                 
      U(4,2)=( XI*Q/R3     +ALP3*AJ4*SDCD)/PI2                              
      U(5,2)=( ET*Q/R3+QY  +ALP3*AJ5*SDCD)/PI2                              
      U(6,2)=(-Q2/R3       +ALP3*AJ6*SDCD)/PI2                              
      U(7,2)=(-EY          +ALP3*AJ1*SDCD)/PI2                              
      U(8,2)=(-ET*GY-XY*SD +ALP3*AJ2*SDCD)/PI2                              
      U(9,2)=( Q*GY        +ALP3*AJ3*SDCD)/PI2                              
      U(10,2)=(-EZ          -ALP3*AK3*SDCD)/PI2                              
      U(11,2)=(-ET*GZ-XY*CD -ALP3*XI*D11*SDCD)/PI2                           
      U(12,2)=( Q*GZ        -ALP3*AK4*SDCD)/PI2                       
      
C========================================                               
C=====  TENSILE-FAULT CONTRIBUTION  =====                               
C========================================                               
      U(1,3)= (Q*QY           -ALP3*AI3*SDSD)/PI2                           
      U(2,3)= (Q*QX           +ALP3*XI/RD*SDSD)/PI2                         
      U(3,3)= (ET*QX+XI*QY-TT -ALP3*AI4*SDSD)/PI2                           
      U(4,3)=(-XI*Q2*Y32 -ALP3*AJ4*SDSD)/PI2                                
      U(5,3)=(-Q2/R3     -ALP3*AJ5*SDSD)/PI2                                
      U(6,3)= (Q*Q2*Y32  -ALP3*AJ6*SDSD)/PI2                                
      U(7,3)= (Q*FY -ALP3*AJ1*SDSD)/PI2                                     
      U(8,3)= (Q*GY -ALP3*AJ2*SDSD)/PI2                                     
      U(9,3)=(-Q*HY -ALP3*AJ3*SDSD)/PI2                                     
      U(10,3)= (Q*FZ +ALP3*AK3*SDSD)/PI2                                     
      U(11,3)= (Q*GZ +ALP3*XI*D11*SDSD)/PI2                                  
      U(12,3)=(-Q*HZ +ALP3*AK4*SDSD)/PI2                                     

      RETURN                                                            
      END                                                               
      
      SUBROUTINE  UC(XI,ET,Q,Z,U)                     
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION U(12,3)
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
C*****   U(12,3) : DISPLACEMENT AND THEIR DERIVATIVES                     
C                                                                       
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,
     &		 SDSD,CDCD,SDCD,S2D,C2D  
      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,  
     &           X11,Y11,X32,Y32,EY,EZ,FY,FZ,GY,GZ,HY,HZ                                
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
      DO 111 L=1,3 
      DO 111  I=1,12                                                    
  111 U(I,L)=F0                                                  
  
C======================================                                 
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 
C======================================                                 
      U(1,1)= (ALP4*XY*CD           
     &		  -ALP5*XI*Q*Z32)/PI2                     
      U(2,1)= (ALP4*(CD/R+F2*QY*SD) 
     &		  -ALP5*C*Q/R3)/PI2                       
      U(3,1)= (ALP4*QY*CD           
     &		  -ALP5*(C*ET/R3-Z*Y11+XI2*Z32))/PI2      
      U(4,1)= (ALP4*Y0*CD                  
     &		  -ALP5*Q*Z0)/PI2                  
      U(5,1)=(-ALP4*XI*(CD/R3+F2*Q*Y32*SD) 
     &		  +ALP5*C*XI*QR)/PI2               
      U(6,1)=(-ALP4*XI*Q*Y32*CD            
     &		  +ALP5*XI*(F3*C*ET/R5-QQ))/PI2    
      U(7,1)=(-ALP4*XI*PPY*CD    -ALP5*XI*QQY)/PI2                          
      U(8,1)= (ALP4*F2*(D/R3-Y0*SD)*SD-Y/R3*CD                         
     &            -ALP5*(CDR*SD-ET/R3-C*Y*QR))/PI2           
      U(9,1)=(-ALP4*Q/R3+YY0*SD  
     &	      +ALP5*(CDR*CD+C*D*QR-(Y0*CD+Q*Z0)*SD))/PI2 
      U(10,1)= (ALP4*XI*PPZ*CD    -ALP5*XI*QQZ)/PI2                          
      U(11,1)= (ALP4*F2*(Y/R3-Y0*CD)*SD+D/R3*CD 
     &		  -ALP5*(CDR*CD+C*D*QR))/PI2   
      U(12,1)=         (YY0*CD    
     &	     -ALP5*(CDR*SD-C*Y*QR-Y0*SDSD+Q*Z0*CD))/PI2  

C======================================                                 
C=====    DIP-SLIP CONTRIBUTION   =====                                 
C======================================                                 
      U(1,2)= (ALP4*CD/R -QY*SD -ALP5*C*Q/R3)/PI2                           
      U(2,2)= (ALP4*Y*X11       -ALP5*C*ET*Q*X32)/PI2                       
      U(3,2)=     (-D*X11-XY*SD -ALP5*C*(X11-Q2*X32))/PI2                   
      U(4,2)=(-ALP4*XI/R3*CD 
     &		+ALP5*C*XI*QR +XI*Q*Y32*SD)/PI2                
      U(5,2)=(-ALP4*Y/R3     +ALP5*C*ET*QR)/PI2                             
      U(6,2)= (D/R3-Y0*SD +ALP5*C/R3*(F1-F3*Q2/R2))/PI2                  
      U(7,2)=(-ALP4*ET/R3+Y0*SDSD 
     &		 -ALP5*(CDR*SD-C*Y*QR))/PI2                
      U(8,2)= (ALP4*(X11-Y*Y*X32) 
     &		 -ALP5*C*((D+F2*Q*CD)*X32-Y*ET*Q*X53))/PI2 
      U(9,2)=  (XI*PPY*SD+Y*D*X32 
     &		 +ALP5*C*((Y+F2*Q*SD)*X32-Y*Q2*X53))/PI2   
      U(10,2)=      (-Q/R3+Y0*SDCD 
     &		 -ALP5*(CDR*CD+C*D*QR))/PI2                
      U(11,2)= (ALP4*Y*D*X32       
     & 		-ALP5*C*((Y-F2*Q*SD)*X32+D*ET*Q*X53))/PI2 
      U(12,2)=(-XI*PPZ*SD+X11-D*D*X32
     &		-ALP5*C*((D-F2*Q*CD)*X32-D*Q2*X53))/PI2 

C========================================                               
C=====  TENSILE-FAULT CONTRIBUTION  =====                               
C========================================                               
      U(1,3)=(-ALP4*(SD/R+QY*CD)   -ALP5*(Z*Y11-Q2*Z32))/PI2                
      U(2,3)= (ALP4*F2*XY*SD+D*X11 -ALP5*C*(X11-Q2*X32))/PI2                
      U(3,3)= (ALP4*(Y*X11+XY*CD)  
     &		+ALP5*Q*(C*ET*X32+XI*Z32))/PI2           
      U(4,3)= (ALP4*XI/R3*SD+XI*Q*Y32*CD
     &		+ALP5*XI*(F3*C*ET/R5-F2*Z32-Z0))/PI2
      U(5,3)= (ALP4*F2*Y0*SD-D/R3 +ALP5*C/R3*(F1-F3*Q2/R2))/PI2             
      U(6,3)=(-ALP4*YY0           -ALP5*(C*ET*QR-Q*Z0))/PI2                 
      U(7,3)= (ALP4*(Q/R3+Y0*SDCD)   
     &		+ALP5*(Z/R3*CD+C*D*QR-Q*Z0*SD))/PI2    
      U(8,3)=(-ALP4*F2*XI*PPY*SD-Y*D*X32                               
     &       	+ALP5*C*((Y+F2*Q*SD)*X32-Y*Q2*X53))/PI2            
      U(9,3)=(-ALP4*(XI*PPY*CD-X11+Y*Y*X32)                            
     &          +ALP5*(C*((D+F2*Q*CD)*X32-Y*ET*Q*X53)+XI*QQY))/PI2 
      U(10,3)=  (-ET/R3+Y0*CDCD 
     &		-ALP5*(Z/R3*SD-C*Y*QR-Y0*SDSD+Q*Z0*CD))/PI2  
      U(11,3)= (ALP4*F2*XI*PPZ*SD-X11+D*D*X32                           
     &          -ALP5*C*((D-F2*Q*CD)*X32-D*Q2*X53))/PI2            
      U(12,3)= (ALP4*(XI*PPZ*CD+Y*D*X32)                                
     &          +ALP5*(C*((Y-F2*Q*SD)*X32+D*ET*Q*X53)+XI*QQZ))/PI2 
      RETURN                                                            
      END                                                      
      
      SUBROUTINE  DCCON0(ALPHA,DIP)                                    
      IMPLICIT REAL*8 (A-H,O-Z)                                         
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
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,
     &           SD,CD,SDSD,CDCD,SDCD,S2D,C2D  
      DATA F0,F1,F2,PI2/0.D0,1.D0,2.D0,6.283185307179586D0/             
      DATA EPS/1.D-6/                                                   
C-----                                                                  
      ALP1=(F1-ALPHA)/F2                                                
      ALP2= ALPHA/F2                                                    
      ALP3=(F1-ALPHA)/ALPHA                                             
      ALP4= F1-ALPHA                                                    
      ALP5= ALPHA                                                       
C-----                                                                  
      P18=PI2/360.D0                                                    
      SD=DSIN(DIP*P18)                                                  
      CD=DCOS(DIP*P18)                                                  
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
      
      SUBROUTINE  DCCON2(XI,ET,Q,SD,CD)                                 
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C                                                                       
C********************************************************************** 
C*****   CALCULATE STATION GEOMETRY CONSTANTS FOR FINITE SOURCE   ***** 
C********************************************************************** 
C                                                                       
C***** INPUT                                                            
C*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  
C*****   SD,CD   : SIN, COS OF DIP-ANGLE                                
C                                                                       
C### CAUTION ### IF XI,ET,Q ARE SUFFICIENTLY SMALL, THEY ARE SET TO ZERO
C                                                                       
      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,  
     *       X11,Y11,X32,Y32,EY,EZ,FY,FZ,GY,GZ,HY,HZ                                
      DATA  F0,F1,F2,EPS/0.D0,1.D0,2.D0,1.D-6/                          
      DATA PI2/6.283185307179586D0/                     
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
C        TT=F0                                                           
 		 TT=PI2*.25D0
 		 IF((XI*ET).LT.0.D0) TT=-TT
      ELSE                                                              
        TT=DATAN(XI*ET/(Q*R))                                           
      ENDIF                                                             
C-----                                                                  
      IF(XI.LT.F0 .AND. Q.EQ.F0 .AND. ET.EQ.F0) THEN                    
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
      IF(ET.LT.F0 .AND. Q.EQ.F0 .AND. XI.EQ.F0) THEN                    
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
      EY=SD/R-Y*Q/R3                                                    
      EZ=CD/R+D*Q/R3                                                    
      FY=D/R3+XI2*Y32*SD                                                
      FZ=Y/R3+XI2*Y32*CD                                                
      GY=F2*X11*SD-Y*Q*X32                                              
      GZ=F2*X11*CD+D*Q*X32                                              
      HY=D*Q*X32+XI*Q*Y32*SD                                            
      HZ=Y*Q*X32+XI*Q*Y32*CD                                            
      RETURN                                                            
      END                                                               
