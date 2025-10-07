MODULE hacapk_c_interface
  !
  ! C interface for HACApK 
  !

  use m_HACApK_base
  use m_HACApK_solve
  use m_HACApK_use

  use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_f_pointer
  
  implicit none 

  TYPE hacapk_struct
     type(st_HACApK_lcontrol) :: st_ctl
     type(st_HACApK_leafmtxp) :: st_leafmtxp
     type(st_HACApK_calc_entry) :: st_bemv
     logical hmat_init,coord_init
  END TYPE hacapk_struct

CONTAINS

  SUBROUTINE init_hacapk_struct(lc_pointer, nd) BIND(C, name='init_hacapk_struct')
    USE ISO_C_BINDING
    include 'mpif.h'
    TYPE(C_PTR), INTENT(OUT) :: lc_pointer
    INTEGER, INTENT(IN) :: nd
    TYPE(hacapk_struct), POINTER :: lf_struct
    integer lrtrn
    ALLOCATE(lf_struct)

    lf_struct%hmat_init = .false.
    lf_struct%coord_init = .false.
    !
    lf_struct%st_bemv%nd = nd
    !
    allocate(lf_struct%st_bemv%xcol(nd),lf_struct%st_bemv%ycol(nd),lf_struct%st_bemv%zcol(nd))
    !
    lrtrn = HACApK_init(nd,lf_struct%st_ctl,lf_struct%st_bemv,MPI_COMM_WORLD)
    !lrtrn = HACApK_init(st_bemv%nd,st_ctl,st_bemv)

    lf_struct%st_ctl%param(61)=1            ! matrix normalization

    lc_pointer = C_LOC(lf_struct)
    
  END SUBROUTINE init_hacapk_struct

  
  SUBROUTINE deallocate_hacapk_struct(lc_pointer) BIND(C, name='deallocate_hacapk_struct')
    USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN) :: lc_pointer
    TYPE(hacapk_struct), POINTER :: lf_struct
    integer lrtrn
    ! Associate the C handle with a Fortran pointer.
    CALL C_F_POINTER(lc_pointer, lf_struct)
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
  END SUBROUTINE deallocate_hacapk_struct
  !
  ! assign coordinates, x,y,z have to be nd dimensional vectors
  !
  SUBROUTINE set_hacapk_struct_coord(lc_pointer, x,y,z) BIND(C, name='set_hacapk_struct_coord')
    USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN) :: lc_pointer
    REAL*8, INTENT(IN),dimension(:) :: x,y,z
    TYPE(hacapk_struct), POINTER :: lf_struct
    integer i
    
    CALL C_F_POINTER(lc_pointer, lf_struct) ! Associate the C handle with a Fortran pointer.
    !do i=1,lf_struct%nd
    lf_struct%st_bemv%xcol(:) = x(:)
    lf_struct%st_bemv%ycol(:) = y(:)
    lf_struct%st_bemv%zcol(:) = z(:)
    !end do
    lf_struct%coord_init = .true. 
  END SUBROUTINE set_hacapk_struct_coord


  
  SUBROUTINE make_hacapk_struct_hmat(lc_pointer, ztol) BIND(C, name='make_hacapk_struct_hmat')
    USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN) :: lc_pointer
    real*8, intent(in) :: ztol
    !
    REAL*8, dimension(:,:),allocatable :: coord
    TYPE(hacapk_struct), POINTER :: lf_struct
    integer lrtrn

    CALL C_F_POINTER(lc_pointer, lf_struct)
    if(.not.lf_struct%coord_init)then
       print *,'make_hacapk_struct_hmat: coord not initialized'
       stop
    endif
    allocate(coord(lf_struct%st_bemv%nd,3))
    coord(:,1) = lf_struct%st_bemv%xcol(:)
    coord(:,2) = lf_struct%st_bemv%ycol(:)
    coord(:,2) = lf_struct%st_bemv%zcol(:)
    !
    lrtrn=HACApK_generate(lf_struct%st_leafmtxp, lf_struct%st_bemv, lf_struct%st_ctl, coord, ztol)
    if(lrtrn.gt.0)then
       print *,'make_hacapk_struct_hmat: HACApK_generate returned code ',lrtrn
       stop
    endif 
    lf_struct%hmat_init = .true. 
    deallocate(coord)
  END SUBROUTINE make_hacapk_struct_hmat

  SUBROUTINE hacapk_assemble_dense_mat(lc_pointer, Ad, n) BIND(C, name='hacapk_assemble_dense_mat')
    USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN) :: lc_pointer
    integer, intent(in) :: n
    REAL*8, dimension(:,:),intent(out) :: Ad
    TYPE(hacapk_struct), POINTER :: lf_struct
    integer i,j

    CALL C_F_POINTER(lc_pointer, lf_struct)
    do i=1,n
       do j=1,n
          Ad(i,j) = HACApK_entry_ij(i, j, lf_struct%st_bemv)
       enddo
    enddo
    
  END SUBROUTINE hacapk_assemble_dense_mat

  ! x = A\b for dense matrix
  subROUTINE hacapk_solve_dense(Ad, n, b, x) BIND(C, name='hacapk_solve_dense')
    USE ISO_C_BINDING
    integer, intent(in) :: n
    REAL*8, dimension(:,:),intent(in) :: Ad
    real*8, dimension(:), intent(in) :: b
    real*8, dimension(:), intent(out) :: x
    integer info
    REAL*8, dimension(:,:),allocatable :: Ause
    REAL*8, dimension(:),allocatable :: buse
    TYPE(hacapk_struct), POINTER :: lf_struct
    integer,dimension(:),allocatable :: ipiv
    
    allocate(Ause(n,n),buse(n),ipiv(n))

    Ause=Ad
    buse=b
    call dgesv(n,1,Ause,n,ipiv,buse,n,info)
    x=buse
    deallocate(Ause,buse,ipiv)
    
  END SUBROUTINE hacapk_solve_dense
  
  
  
  !
  ! b = Ax
  SUBROUTINE hacapk_mult_Ax(lc_pointer, x, b) BIND(C, name='hacapk_mult_Ax')
    USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN) :: lc_pointer
    REAL*8, intent(in) :: x(:)
    REAL*8, intent(out) :: b(:)

    TYPE(hacapk_struct), POINTER :: lf_struct
    integer lrtrn

    CALL C_F_POINTER(lc_pointer, lf_struct)
    if(.not.lf_struct%hmat_init)then
       print *,'hacapk_mult_Ax: H matrix not initialized'
       stop
    endif

    ! version 1
    !lrtrn = HACApK_adot_pmt_lfmtx_hyp(lf_struct%st_leafmtxp,lf_struct%st_bemv,lf_struct%st_ctl,b,x)
    ! version 2 
    lrtrn = HACApK_adot_pmt_lfmtx_p(lf_struct%st_leafmtxp,lf_struct%st_bemv,lf_struct%st_ctl,b,x)
    if(lrtrn.gt.0)then
       print *,'hacapk_mult_Ax:  HACApK_adot_pmt returned code ',lrtrn
       stop
    endif 

  END SUBROUTINE hacapk_mult_Ax

  !
  ! x = A\b for H matrix
  !
  SUBROUTINE hacapk_solve_Ab(lc_pointer, b, x, ztol) BIND(C, name='hacapk_solve_Ab')
    USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN) :: lc_pointer
    REAL*8, intent(out) :: x(:)
    REAL*8, intent(in) :: b(:)
    real*8, intent(in) :: ztol
    REAL*8, dimension(:,:),allocatable :: coord
    REAL*8, dimension(:),allocatable :: buse
   TYPE(hacapk_struct), POINTER :: lf_struct
    integer lrtrn

    CALL C_F_POINTER(lc_pointer, lf_struct)

    if(.not.lf_struct%coord_init)then
       print *,'hacapk_solve_Ab: coord not initialized'
       stop
    endif
    allocate(buse(lf_struct%st_bemv%nd))
    allocate(coord(lf_struct%st_bemv%nd,3))
    coord(:,1) = lf_struct%st_bemv%xcol(:)
    coord(:,2) = lf_struct%st_bemv%ycol(:)
    coord(:,2) = lf_struct%st_bemv%zcol(:)

    x = 0.d0                    ! starting guess
    buse = b
    lrtrn=HACApK_gensolv(lf_struct%st_leafmtxp,lf_struct%st_bemv,lf_struct%st_ctl,coord,buse,x,ztol)
    if(lrtrn.gt.0)then
       print *,'hacapk_solve_Ab:  HACApK_gensolv returned code ',lrtrn
       stop
    endif
    deallocate(coord,buse)
  END SUBROUTINE hacapk_solve_Ab



  
END MODULE hacapk_c_interface
