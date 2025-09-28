program Example_Using_HACApK
  use m_HACApK_base
  use m_HACApK_solve
  use m_HACApK_use
  implicit none
  include 'mpif.h'
  type(st_HACApK_lcontrol) :: st_ctl
  type(st_HACApK_leafmtxp) :: st_leafmtxp
  type(st_HACApK_calc_entry) :: st_bemv
  
  real*8,dimension(:,:),allocatable :: coord,mdens,mdens_backup
  real*8,dimension(:),allocatable :: xd,xh,bd,bh,bsave
  real*8 :: ztol
  integer i,j
  integer ierr,lrtrn,icomm,info
  integer,dimension(:),allocatable :: ipiv

  call MPI_Init (ierr);icomm = MPI_COMM_WORLD
  call random_seed()

  
  st_bemv%ndim = 3                 ! dimension
  st_bemv%n = 1000                        ! number of points

  print *,'initializing'
  lrtrn = HACApK_init(st_bemv%n,st_ctl,st_bemv,icomm)
  !lrtrn = HACApK_init(st_bemv%n,st_ctl,st_bemv)

  print *,'allocating'
  allocate(st_bemv%xcol(st_bemv%n),st_bemv%ycol(st_bemv%n),st_bemv%zcol(st_bemv%n))
  allocate(coord(st_bemv%n,st_bemv%ndim),xd(st_bemv%n),xh(st_bemv%n),bh(st_bemv%n),bd(st_bemv%n),bsave(st_bemv%n))
  allocate(mdens(st_bemv%n,st_bemv%n),mdens_backup(st_bemv%n,st_bemv%n),ipiv(st_bemv%n))
  st_bemv%scale = 1.d0


  call random_number(st_bemv%xcol)
  call random_number(st_bemv%ycol)
  call random_number(st_bemv%zcol)
  ! why can't we just use the full coord structure?
  coord(:,1)=st_bemv%xcol(:)
  coord(:,2)=st_bemv%ycol(:)
  coord(:,3)=st_bemv%zcol(:)

  print *,'loading'

  ztol=1.0e-6
  print *,'making H '
  lrtrn=HACApK_generate(st_leafmtxp,st_bemv,st_ctl,coord,ztol)

  stop 
  print *,'making dense matrix'
  do i=1,st_bemv%n
     do j=1,st_bemv%n
        mdens(i,j) = HACApK_entry_ij(i, j, st_bemv)
     enddo
  enddo 

  !
  !
  ! RHS and solution guess
  call random_number(bsave)


  !
  print *,'solve H'
  xh = 0.d0
  bh = bsave
  lrtrn=HACApK_gensolv(st_leafmtxp,st_bemv,st_ctl,coord,bh,xh,ztol)

  print *,'dense solve'
  mdens_backup = mdens
  bd = bsave
  !print *,dsol
  call dgesv(st_bemv%n,1,mdens,st_bemv%n,ipiv,bd,st_bemv%n,info)
  xd = bd
  ! check dens solution
  print *,'dense solve solution error:',norm2(bsave-matmul(mdens_backup,xd))
  !
  !
  print *,'H vs dense solve x error:',norm2(xd-xh)


  ! multiply dense
  call random_number(xd)
  bd = matmul(mdens_backup,xd) ! dense b = Ax 
  !
  xh=xd
  lrtrn = HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp,st_bemv,st_ctl,bh,xh)
  print *,'H vs dense multiplication b error', norm2(bh-bd)
  

  lrtrn=HACApK_free_leafmtxp(st_leafmtxp)
  lrtrn=HACApK_finalize(st_ctl)
  
  call MPI_Finalize (ierr)
end program
