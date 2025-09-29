!=====================================================================*
!                                                                     *
!   Software Name : HACApK                                            *
!         Version : 0.2.0                                             *
!                                                                     *
!   License                                                           *
!     This file is part of HACApK.                                    *
!     HACApK is a free software, you can use it under the terms       *
!     of The MIT License (MIT). See LICENSE file and User's guide     *
!     for more details.                                               *
!                                                                     *
!   ppOpen-HPC project:                                               *
!     Open Source Infrastructure for Development and Execution of     *
!     Large-Scale Scientific Applications on Post-Peta-Scale          *
!     Supercomputers with Automatic Tuning (AT).                      *
!                                                                     *
!   Sponsorship:                                                      *
!     Japan Science and Technology Agency (JST), Basic Research       *
!     Programs: CREST, Development of System Software Technologies    *
!     for post-Peta Scale High Performance Computing.                 *
!                                                                     *
!   Copyright (c) 2014 <Akihiro Ida and Takeshi Iwashita>             *
!                                                                     *
!=====================================================================*
module m_HACApK_calc_entry_ij

!*** type :: st_HACApK_calc_entry
  type :: st_HACApK_calc_entry
     integer :: nd,lp61
     real*8,pointer :: ao(:)
     ! user defined
     real*8 :: scale
     real*8,pointer :: xcol(:),ycol(:),zcol(:)
  
  end type st_HACApK_calc_entry

  public :: HACApK_entry_ij

contains
  !***HACApK_entry_ij
  real*8 function HACApK_entry_ij(i, j, st_bemv)
    implicit none
    type(st_HACApK_calc_entry), intent(in) :: st_bemv
    integer, intent(in) :: i,j
    !
    real*8 dist2

    dist2=      (st_bemv%xcol(i)-st_bemv%xcol(j))**2
    dist2=dist2+(st_bemv%ycol(i)-st_bemv%ycol(j))**2
    dist2=dist2+(st_bemv%zcol(i)-st_bemv%zcol(j))**2
    
    HACApK_entry_ij= 1.0d0/(1d-12 + dist2/st_bemv%scale)**0.5 ! made up kernel
    
    
    
  end function HACApK_entry_ij

endmodule m_HACApK_calc_entry_ij
