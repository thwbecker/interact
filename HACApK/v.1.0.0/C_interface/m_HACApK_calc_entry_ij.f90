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
  use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double
  implicit none

  !
  ! this is an example of how to use a C function
  !

  
!*** type :: st_HACApK_calc_entry
  type :: st_HACApK_calc_entry
     ! has to be here
     integer :: nd,lp61
     real*8,pointer :: ao(:)
     !
     ! user defined
     type(c_ptr) :: ckernel_par
     real*8,pointer :: xcol(:),ycol(:),zcol(:)
  
  end type st_HACApK_calc_entry

  ! change here and in C code as well as header
  interface
     function ckernel_func(i,j,kp) bind(c, name='ckernel_func')
       import :: c_ptr, c_int, c_double 
       integer(c_int), value, intent(in) :: i,j !input, pass by value
       type(c_ptr), value, intent(in) :: kp     !input, holds all parameters
       real(c_double) :: ckernel_func           !output is a real*8
     end function ckernel_func
  end interface
  
  public :: HACApK_entry_ij

contains
  !
  !***HACApK_entry_ij
  !
  real*8 function HACApK_entry_ij(i, j, st_bemv)
    implicit none
    type(st_HACApK_calc_entry), intent(in) :: st_bemv
    integer, intent(in) :: i,j
    ! call C function with pointer as part of st_bemv structure
    HACApK_entry_ij = ckernel_func(i-1, j-1, st_bemv%ckernel_par)
    
  end function HACApK_entry_ij

endmodule m_HACApK_calc_entry_ij
