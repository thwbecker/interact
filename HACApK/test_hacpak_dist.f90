program Example_Using_HACApK
  !
  ! does not work with HBI copied HACApK
  !
  use m_HACApK_base
  use m_HACApK_solve
  use m_HACApK_use
  implicit none
  include 'mpif.h'
  !
  type(st_HACApK_lcontrol) :: st_ctl
  type(st_HACApK_leafmtxp) :: st_leafmtxp
  type(st_HACApK_calc_entry) :: st_bemv
  !
  type(st_HACApK_latticevec) :: st_lat_out,st_lat_in
  type(st_HACApK_LHp) :: st_LHp
  real(8),allocatable::wws(:)
  !
  real*8,dimension(:,:),allocatable :: coord,mdens,mdens_backup
  real*8,dimension(:),allocatable :: xd,xh,bd,bh,bsave
  real*8 :: ztol
  integer i,j,ni,ind,ng,nl
  integer ierr,lrtrn,info,imode,mpi_np,mpi_rank
  integer,dimension(:),allocatable :: ipiv


  imode = 1                     ! 1: old multiplication 2: lattive
  
  call MPI_Init (ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpi_np,ierr )
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpi_rank,ierr )
  
  call random_seed()
  !
  !
  ng = 20
  st_bemv%nd = ng                        ! global number of points

  if(mpi_rank == 0)then
     print *,'initializing'
  endif
  lrtrn = HACApK_init(ng,st_ctl,st_bemv,MPI_COMM_WORLD)

  
  st_ctl%param(8)=20            ! needed?
 
  
  if(mpi_rank == 0)print *,'allocating'
  allocate(st_bemv%xcol(ng),st_bemv%ycol(ng),st_bemv%zcol(ng))
  allocate(coord(ng,3))
  allocate(xd(ng),xh(ng),bh(ng),bd(ng),bsave(ng))
  
  if(mpi_rank == 0)then
     allocate(mdens(ng,ng),mdens_backup(ng,ng),ipiv(ng))
  endif


  st_bemv%scale = 1.d0
  
  
  call random_number(st_bemv%xcol)
  call random_number(st_bemv%ycol)
  call random_number(st_bemv%zcol)
  ! why can't we just use the full coord structure?
  coord(:,1)=st_bemv%xcol(:)
  coord(:,2)=st_bemv%ycol(:)
  coord(:,3)=st_bemv%zcol(:)

  
  if(mpi_rank == 0)print *,'making H '
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  ztol=1.0e-6
  lrtrn=HACApK_generate(st_leafmtxp,st_bemv,st_ctl,coord,ztol)

  if(imode.eq.2)then
     lrtrn=HACApK_construct_LH(st_LHp,st_leafmtxp,st_bemv,st_ctl,coord,ztol)
     allocate(wws(st_leafmtxp%ndlfs))
     lrtrn=HACApK_gen_lattice_vector(st_lat_out,st_leafmtxp,st_ctl)
     lrtrn=HACApK_gen_lattice_vector(st_lat_in, st_leafmtxp,st_ctl)
     nl = st_lat_out%ndc        !local node number
  else
     nl = ng                    ! how to determine local nodes?
  endif
  print *,'nc  ',mpi_np, ' rank ',mpi_rank, ' n ',ng, ' nl ',nl
  !
  if(mpi_rank == 0)then
     print *,'making dense matrix' ! only do on lead node
     do i=1,ng
        do j=1,ng
           mdens_backup(i,j) = HACApK_entry_ij(i, j, st_bemv)
           !print *,i,j, mdens_backup(i,j)
        enddo
     enddo
  endif
  
  !
  !
  ! RHS
  !call random_number(bsave)
  bsave = 1.d3
  !
  ! x = A\b
  !
  print *,'solve H'
  xh = 0.d0                     ! initial guess
  bh = bsave
  lrtrn=HACApK_gensolv(st_leafmtxp,st_bemv,st_ctl,coord,bh,xh,ztol)
  !
  !
  if(mpi_rank == 0)then
     print *,'dense solve'
     mdens = mdens_backup 
     bd = bsave
     !print *,dsol
     ! mdens will get overwritten
     call dgesv(ng,1,mdens,ng,ipiv,bd,ng,info)
     xd = bd
     ! check dens solution
     print *,'dense solve solution error:',norm2(bsave-matmul(mdens_backup,xd))
     ! send to all nodes
     call MPI_Bcast(xd,ng,MPI_REAL8,0,MPI_COMM_WORLD)
  endif
  if(ng.lt.30)then
     ! print solution
     do i=1,ng
        print *,i,' H ',xh(i),' D ',xd(i),' diff ',abs(xh(i)-xd(i))
     enddo
  endif
  if(mpi_rank == 0)print *,'H vs dense solve x error:',norm2(xd-xh)
  
  
  !
  ! multiply b = A x
  !
  xd = 1.d0

  if(mpi_rank == 0)then
     ! multiply dense

     bd = matmul(mdens_backup,xd) ! dense b = Ax
     !call MPI_Bcast(bd,ng,MPI_REAL8,0,MPI_COMM_WORLD)
     print *,'dense bd norm ',norm2(bd)
  endif
  !
  ! H
  xh=xd
  if(imode.eq.1)then
     lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp,st_bemv,st_ctl,bh,xh)
     print *,'H bh norm ',norm2(bh)
  else
     do i=1,nl                  ! assign from global to local
        ind = st_lat_in%lodc(i)
        st_lat_in%vs(i) = xh(ind)
     enddo
     call HACApK_adot_lattice_hyp(st_lat_out,st_LHp,st_ctl,wws,st_lat_in)
     do i=1,nl                  ! assign from local to global
        ind = st_lat_out%lodc(i)
        bh(ind) = st_lat_out%vs(i)
     enddo
     call MPI_Bcast(bh,ng,MPI_REAL8,0,MPI_COMM_WORLD)
  endif
  !
  !
  !
  if(mpi_rank == 0) then 
     if(ng.lt.30)then
        do i=1,ng
           print *,i,' H ',bh(i),' D ',bd(i), ' diff ',abs(bh(i)-bd(i))
        enddo
     endif 
     print *,'mode ',imode, ' H vs dense multiplication b error', norm2(bh-bd)
  endif

  lrtrn=HACApK_free_leafmtxp(st_leafmtxp)
  lrtrn=HACApK_finalize(st_ctl)
  
  call MPI_Finalize (ierr)
end program
