MODULE hacapk_testing_routines
  !
  ! shared routines of the C and F90 testing programs for HACApK
  !
  use, intrinsic :: iso_c_binding, only:  c_int, c_double, c_ptr, c_loc, c_f_pointer

  use m_HACApK_base
  use m_HACApK_solve
  use m_HACApK_use

  implicit none 


CONTAINS

  subroutine hacapk_assign_random_coord(x,y,z,n) BIND(C, name='hacapk_assign_random_coord')
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
  end subroutine hacapk_assign_random_coord
  
END MODULE hacapk_testing_routines
