program Example_Using_HACApK
  !
  ! does not work with HBI copied HACApK
  !
  use m_HACApK_base
  use m_HACApK_solve
  use m_HACApK_use
  implicit none
  include 'mpif.h'
  type(st_HACApK_lcontrol) :: st_ctl
  type(st_HACApK_leafmtxp) :: st_leafmtxp
  type(st_HACApK_calc_entry) :: st_bemv
  
  real*8,dimension(:,:),allocatable :: coord,mdens
  real*8,dimension(:),allocatable :: xd,xh,bd,bh,bsave
  real*8 :: ztol
  integer i,j,ni
  integer ierr,lrtrn,info
  integer,dimension(:),allocatable :: ipiv

  call MPI_Init (ierr)
  call random_seed()
  !
  !
  st_bemv%nd = 20                        ! number of points

  print *,'initializing'
  lrtrn = HACApK_init(st_bemv%nd,st_ctl,st_bemv,MPI_COMM_WORLD)
  !lrtrn = HACApK_init(st_bemv%nd,st_ctl,st_bemv)
  
 
  
  print *,'allocating'
  allocate(st_bemv%xcol(st_bemv%nd),st_bemv%ycol(st_bemv%nd),st_bemv%zcol(st_bemv%nd))
  allocate(coord(st_bemv%nd,3))
  allocate(xd(st_bemv%nd),xh(st_bemv%nd),bh(st_bemv%nd),bd(st_bemv%nd),bsave(st_bemv%nd))
  allocate(mdens(st_bemv%nd,st_bemv%nd),ipiv(st_bemv%nd))


  st_bemv%scale = 1.d0
  
  
  call random_number(st_bemv%xcol)
  call random_number(st_bemv%ycol)
  call random_number(st_bemv%zcol)
  ! why can't we just use the full coord structure?
  coord(:,1)=st_bemv%xcol(:)
  coord(:,2)=st_bemv%ycol(:)
  coord(:,3)=st_bemv%zcol(:)

  ztol=1.0e-6

  print *,'making H '
  lrtrn=HACApK_generate(st_leafmtxp,st_bemv,st_ctl,coord,ztol)

  st_ctl%param(61)=1            ! matrix normalization
 
  
  !
  print *,'making dense matrix'
  do i=1,st_bemv%nd
     do j=1,st_bemv%nd
        mdens(i,j) = HACApK_entry_ij(i, j, st_bemv)
        !print *,i,j, mdens(i,j)
     enddo
  enddo 

  
 
  !
  ! multiply b = A x
  !
  ! multiply dense
  !call random_number(xd)
  xd = 1.d0
  bd = matmul(mdens,xd) ! dense b = Ax 
  !
  xh=xd
  ! version 1 - does not work 
  !lrtrn = HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp,st_bemv,st_ctl,bh,xh)
  ! version 2 - does not work 
  lrtrn = HACApK_adot_pmt_lfmtx_p(st_leafmtxp,st_bemv,st_ctl,bh,xh)
  !
  !
  !

  if(st_bemv%nd.lt.30)then
     do i=1,st_bemv%nd
        print *,i,' H ',bh(i),' D ',bd(i), ' diff ',abs(bh(i)-bd(i))
     enddo
  endif
  print *,'H vs dense multiplication b error', norm2(bh-bd)
  

  lrtrn=HACApK_free_leafmtxp(st_leafmtxp)
  lrtrn=HACApK_finalize(st_ctl)
  
  call MPI_Finalize (ierr)
end program
