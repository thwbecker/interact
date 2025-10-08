MODULE hacapk_c_interface
  !
  ! C interface for HACApK 
  !
  use, intrinsic :: iso_c_binding, only:  c_int, c_double, c_ptr, c_loc, c_f_pointer

  use m_HACApK_base
  use m_HACApK_solve
  use m_HACApK_use

  implicit none 

  TYPE :: hacapk_chandle_struct
     type(st_HACApK_lcontrol) :: st_ctl
     type(st_HACApK_leafmtxp) :: st_leafmtxp
     type(st_HACApK_calc_entry) :: st_bemv
     integer(c_int) :: ndim
     logical hmat_init, coord_init, init, kernel_init
  END TYPE hacapk_chandle_struct

  ! A pointer to our custom type, used as the opaque handle.
  type(hacapk_chandle_struct), pointer :: hacapk_int_handle => null()


CONTAINS
  !
  ! make the structure and allocate arrays for nd by nd problem
  !
  function cinit_hacapk_struct(nd)  result(c_pointer) BIND(C, name='cinit_hacapk_struct')
    include 'mpif.h'
    type(c_ptr) :: c_pointer   !output is c pointer
    integer(c_int), value, intent(in) :: nd
    integer lrtrn
    
    allocate(hacapk_int_handle)
    
    hacapk_int_handle%hmat_init = .false.
    hacapk_int_handle%coord_init = .false.
    hacapk_int_handle%ndim = nd
    hacapk_int_handle%st_bemv%nd = nd

    allocate(hacapk_int_handle%st_bemv%xcol(nd))
    allocate(hacapk_int_handle%st_bemv%ycol(nd))
    allocate(hacapk_int_handle%st_bemv%zcol(nd))
    !
    lrtrn = HACApK_init(hacapk_int_handle%st_bemv%nd, hacapk_int_handle%st_ctl, &
         hacapk_int_handle%st_bemv, MPI_COMM_WORLD)
    !lrtrn = HACApK_init(st_bemv%nd,st_ctl,st_bemv)
    !
    ! H matrix parameters
    !
    hacapk_int_handle%st_ctl%param(61) = 1            ! matrix normalization
    !
    hacapk_int_handle%init = .true.
    !print *,'initialized hacapk internal structure with n ',hacapk_int_handle%ndim, hacapk_int_handle%init 
    c_pointer = c_loc(hacapk_int_handle)
  END function cinit_hacapk_struct

  ! free arrays
  SUBROUTINE cdeallocate_hacapk_struct(c_pointer) BIND(C, name='cdeallocate_hacapk_struct')
    TYPE(C_PTR), value, INTENT(IN) :: c_pointer
    TYPE(hacapk_chandle_struct), POINTER :: lf_struct
    integer lrtrn
    ! Associate the C handle with a Fortran pointer.
    CALL C_F_POINTER(c_pointer, lf_struct)
    IF (ASSOCIATED(lf_struct)) THEN
       if(lf_struct%coord_init)then
          DEALLOCATE(lf_struct%st_bemv%xcol)
          DEALLOCATE(lf_struct%st_bemv%ycol)
          DEALLOCATE(lf_struct%st_bemv%zcol)
       END IF
       if(lf_struct%hmat_init)then 
          lrtrn = HACApK_free_leafmtxp(lf_struct%st_leafmtxp)
       endif
       lrtrn = HACApK_finalize(lf_struct%st_ctl)
       DEALLOCATE(lf_struct)

       
    END IF
  END SUBROUTINE cdeallocate_hacapk_struct
  !
  ! assign coordinates, x,y,z, have to be nd dimensional vectors
  !
  SUBROUTINE cset_hacapk_struct_coord(c_pointer, x,y,z) BIND(C, name='cset_hacapk_struct_coord')
    TYPE(C_PTR), value, INTENT(in) :: c_pointer
    real(c_double), INTENT(IN) :: x(*),y(*),z(*)
    !
    TYPE(hacapk_chandle_struct), POINTER :: lf_struct
    integer i
    call c_f_pointer(c_pointer, lf_struct) ! Associate the C handle with a Fortran pointer.
    print *,lf_struct%ndim, lf_struct%st_bemv%nd, lf_struct%init
    if(.not.lf_struct%init)then
       print *,'cset_hacapk_struct_coord: structure not initilized'
       stop
    end if

    do i=1,lf_struct%ndim
       !print *,i,lf_struct%st_bemv%nd,x(i),y(i),z(i)
       lf_struct%st_bemv%xcol(i) = x(i)
       lf_struct%st_bemv%ycol(i) = y(i)
       lf_struct%st_bemv%zcol(i) = z(i)
    end do
    lf_struct%coord_init = .true. 
  END SUBROUTINE cset_hacapk_struct_coord

  !
  ! generate the H matrix with kernel provided
  !
  SUBROUTINE cmake_hacapk_struct_hmat(c_pointer, ztol) BIND(C, name='cmake_hacapk_struct_hmat')
    TYPE(C_PTR), value, INTENT(in) :: c_pointer
    real(c_double), value, intent(in) :: ztol
    !
    REAL*8, dimension(:,:),allocatable :: coord
    TYPE(hacapk_chandle_struct), POINTER :: lf_struct
    integer lrtrn,i

    CALL C_F_POINTER(c_pointer, lf_struct)
    if(.not.lf_struct%coord_init)then
       print *,'cmake_hacapk_struct_hmat: coord not initialized'
       stop
    endif
    allocate(coord(lf_struct%st_bemv%nd,3))
    do i=1,lf_struct%st_bemv%nd
       coord(i,1) = lf_struct%st_bemv%xcol(i)
       coord(i,2) = lf_struct%st_bemv%ycol(i)
       coord(i,3) = lf_struct%st_bemv%zcol(i)
       !print *,i,coord(i,:)
    enddo
    if(.not.lf_struct%kernel_init)then
       print *,'cmake_hacapk_struct_hmat: kernel parameters not initialized'
       stop
    endif
    !
    lrtrn=HACApK_generate(lf_struct%st_leafmtxp, lf_struct%st_bemv, lf_struct%st_ctl, coord, ztol)
    if(lrtrn.gt.0)then
       print *,'cmake_hacapk_struct_hmat: HACApK_generate returned code ',lrtrn
       stop
    endif 
    lf_struct%hmat_init = .true.
    print *,'made hmat with ztol ',ztol
    deallocate(coord)
  END SUBROUTINE cmake_hacapk_struct_hmat
  
  !
  ! assemble a dense matrix with same kernel - Fortran style matrix
  !
  SUBROUTINE chacapk_assemble_dense_mat(c_pointer, Ad, n) &
       BIND(C, name='chacapk_assemble_dense_mat')
    TYPE(C_PTR), value, INTENT(in) :: c_pointer
    integer(c_int), value, intent(in) :: n
    real*8, INTENT(out) :: Ad(n,n)
    !
    TYPE(hacapk_chandle_struct), POINTER :: lf_struct
    integer :: i,j
    CALL C_F_POINTER(c_pointer, lf_struct)
    if(.not.lf_struct%kernel_init)then
       print *,'chacapk_assemble_dense_mat: kernel parameters not initialized'
       stop
    endif
    do i=1,n
       do j=1,n
          Ad(i,j) = HACApK_entry_ij(i, j, lf_struct%st_bemv)
          !print *,i,j,lf_struct%st_bemv%xcol(i),lf_struct%st_bemv%ycol(i),lf_struct%st_bemv%zcol(i),&
          !lf_struct%st_bemv%xcol(j),lf_struct%st_bemv%ycol(j),lf_struct%st_bemv%zcol(j),Ad(i,j)
       enddo
    enddo
    
  END SUBROUTINE chacapk_assemble_dense_mat


  subroutine chacapk_set_kernel_par(c_pointer, kscale) BIND(C, name='chacapk_set_kernel_par')
    TYPE(C_PTR), value, INTENT(in) :: c_pointer
    real(c_double), value, INTENT(IN) :: kscale
    !
    TYPE(hacapk_chandle_struct), POINTER :: lf_struct
    call c_f_pointer(c_pointer, lf_struct) ! Associate the C handle with a Fortran pointer.
    !
    ! set kernel parameters
    !
    lf_struct%st_bemv%scale = kscale
    lf_struct%kernel_init = .true. 
  END SUBROUTINE chacapk_set_kernel_par


  !
  ! solve x = A\b for dense matrix, Fortran style matrix
  !
  subROUTINE chacapk_solve_dense(Ad, n, b, x) BIND(C, name='chacapk_solve_dense')
    integer(c_int), value, intent(in) :: n
    REAL*8, dimension(n,n),intent(in) :: Ad
    real*8, dimension(n), intent(in) :: b
    real*8, dimension(n), intent(out) :: x
    !
    integer info,i,j
    REAL*8, dimension(:,:),allocatable :: Ause
    REAL*8, dimension(:),allocatable :: buse
    integer,dimension(:),allocatable :: ipiv
    
    allocate(Ause(n,n),buse(n),ipiv(n))

    Ause=Ad
    buse=b
    call dgesv(n,1,Ause,n,ipiv,buse,n,info)
    x=buse
    deallocate(Ause,buse,ipiv)
    
  END SUBROUTINE chacapk_solve_dense
  
  subroutine chacapk_assign_random_coord(x,y,z,n) BIND(C, name='chacapk_assign_random_coord')
    integer, intent(in) :: n
    real*8, intent(inout) :: x(n),y(n),z(n)
    integer i,itmp
    integer, allocatable :: seed_size(:)
    call random_seed(size=itmp) ! set specific seed
    allocate(seed_size(itmp))
    seed_size=1234
    call random_seed(put=seed_size)

    !print *,'init for random xyz with N ',n
    call random_number(x)
    call random_number(y)
    call random_number(z)
    !do i=1,n
    !print *,i,x(i),y(i),z(i)
    !enddo
    deallocate(seed_size)
  end subroutine chacapk_assign_random_coord
  
  !
  ! b = Ax using H matrix
  !
  SUBROUTINE chacapk_mult_Ax_H(c_pointer, x, b) BIND(C, name='chacapk_mult_Ax_H')
    TYPE(C_PTR), value, INTENT(in) :: c_pointer
    REAL*8, intent(in) :: x(*)
    REAL*8, intent(out) :: b(*)
    !
    TYPE(hacapk_chandle_struct), POINTER :: lf_struct
    integer lrtrn

    CALL C_F_POINTER(c_pointer, lf_struct)
    if(.not.lf_struct%hmat_init)then
       print *,'chacapk_mult_Ax_H: H matrix not initialized'
       stop
    endif

    ! version 1
    !lrtrn = HACApK_adot_pmt_lfmtx_hyp(lf_struct%st_leafmtxp,lf_struct%st_bemv,lf_struct%st_ctl,b,x)
    ! version 2 
    lrtrn = HACApK_adot_pmt_lfmtx_p(lf_struct%st_leafmtxp,lf_struct%st_bemv,lf_struct%st_ctl,b,x)
    if(lrtrn.gt.0)then
       print *,'chacapk_mult_Ax_H:  HACApK_adot_pmt returned code ',lrtrn
       stop
    endif 

  END SUBROUTINE chacapk_mult_Ax_H

  !
  ! x = A\b for H matrix, potentially destructive for b
  !
  SUBROUTINE chacapk_solve_Ab_H(c_pointer, b, x, ztol) BIND(C, name='chacapk_solve_Ab_H')
    TYPE(C_PTR), value, INTENT(in) :: c_pointer
    REAL*8, intent(out) :: x(*)
    REAL*8, intent(in) :: b(*)
    real(c_double), value, intent(in) :: ztol !
    !
    REAL*8, dimension(:,:),allocatable :: coord
    !REAL*8, dimension(:),allocatable :: buse
    TYPE(hacapk_chandle_struct), POINTER :: lf_struct
    integer lrtrn, n, i
    CALL C_F_POINTER(c_pointer, lf_struct)

    if(.not.lf_struct%coord_init)then
       print *,'chacapk_solve_Ab: coord not initialized'
       stop
    endif
    n = lf_struct%st_bemv%nd
    !allocate(buse(n))
    allocate(coord(n,3))
    coord(:,1) = lf_struct%st_bemv%xcol(:)
    coord(:,2) = lf_struct%st_bemv%ycol(:)
    coord(:,2) = lf_struct%st_bemv%zcol(:)
    do i=1,n
       x(i) = 0.d0                    ! starting guess
       !buse(i) = b(i)
    enddo

    
    !lrtrn=HACApK_gensolv(lf_struct%st_leafmtxp,lf_struct%st_bemv,lf_struct%st_ctl,coord,buse,x,ztol)
    lrtrn=HACApK_gensolv(lf_struct%st_leafmtxp,lf_struct%st_bemv,lf_struct%st_ctl,coord,b,x,ztol)
    if(lrtrn.gt.0)then
       print *,'chacapk_solve_Ab:  HACApK_gensolv returned code ',lrtrn
       stop
    endif
    print *,'chacapk_solve_Ab:  called solve with ',ztol
    
    !deallocate(coord,buse)
    deallocate(coord)

  END SUBROUTINE chacapk_solve_Ab_H



  
END MODULE hacapk_c_interface
