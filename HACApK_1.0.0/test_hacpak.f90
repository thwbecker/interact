program Example_Using_HACApK
  use m_HACApK_base
  use m_HACApK_solve
  use m_HACApK_use
  implicit none
  !include 'mpif.h'
  type(st_HACApK_lcontrol) :: st_ctl
  type(st_HACApK_leafmtxp) :: st_leafmtxp
  type(st_HACApK_calc_entry) :: st_bemv
  
  real*8,dimension(:,:),allocatable :: coord
  real*8,dimension(:),allocatable :: rhs,sol
  real*8 :: ztol
  integer i,j
  integer ierr,lrtrn,icomm


  !call MPI_Init (ierr);icomm = MPI_COMM_WORLD
  call random_seed()

  
  st_bemv%ndim = 3                 ! dimension
  st_bemv%n = 100                        ! number of points

  print *,'initializing'
  !lrtrn = HACApK_init(st_bemv%n,st_ctl,st_bemv,icomm)
  lrtrn = HACApK_init(st_bemv%n,st_ctl,st_bemv)

  print *,'allocating'
  allocate(st_bemv%xcol(st_bemv%n),st_bemv%ycol(st_bemv%n),st_bemv%zcol(st_bemv%n))
  allocate(coord(st_bemv%n,st_bemv%ndim),rhs(st_bemv%n),sol(st_bemv%n))
  
  st_bemv%scale =1.d0


  call random_number(st_bemv%xcol)
  call random_number(st_bemv%ycol)
  call random_number(st_bemv%zcol)
  ! why can't we just use the full coord structure?
  coord(:,1)=st_bemv%xcol(:)
  coord(:,2)=st_bemv%ycol(:)
  coord(:,3)=st_bemv%zcol(:)

  print *,'loading'

  ! RHS and solution guess
  rhs = 1.d0
  sol = 0.d0

  
  ztol=1.0e-5
  print *,'making H '
  lrtrn=HACApK_generate(st_leafmtxp,st_bemv,st_ctl,coord,ztol)
  print *,'solve'
  lrtrn=HACApK_solve(st_leafmtxp,st_bemv,st_ctl,rhs,sol,ztol)
  
  
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp)
  lrtrn=HACApK_finalize(st_ctl)
  
  !call MPI_Finalize (ierr)
end program
