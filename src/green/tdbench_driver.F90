! tdbench_driver.F90
!
! Standalone micro-benchmark + accuracy harness for the triangular-dislocation
! half-space kernels tdstresshs() and tddisphs().
!
! Two passes:
!   1. accuracy: evaluate both kernels once at every test point and write the
!      full-precision result vector (stress 6, strain 6, disp 3) to
!      td_results.dat. The build script diffs this file across optimisation
!      regimes to get max / RMS deviation from the strict reference.
!   2. timing: call each kernel NREP times over the NPT-point set, accumulate a
!      checksum (so nothing is optimised away), and report wall time.
!
! Compiled in double precision (-DUSE_DOUBLE_PRECISION); the kernels' C_PREC is
! then real*8, matched here by real*8.
!
program tdbench
  implicit none
  integer, parameter :: npt = 48, nrep = 200000   ! nrep: timing load, raise/lower for stable wall times
  real*8, dimension(3)   :: p1, p2, p3, loc
  real*8                 :: ss, ds, ts
  real*8, dimension(6)   :: stress, strain
  real*8, dimension(3)   :: u
  real*8, dimension(3,npt) :: pts
  real*8                 :: csum
  real*8                 :: tstress, tdisp
  integer                :: i, rep
  integer*8              :: c0, c1, crate

  ! triangular dislocation below the free surface, and a mixed slip vector
  p1 = (/ 0.0d0, 0.0d0, -2.0d0 /)
  p2 = (/ 1.0d0, 0.2d0, -2.0d0 /)
  p3 = (/ 0.1d0, 1.0d0, -3.0d0 /)
  ss = 1.0d0
  ds = 0.6d0
  ts = 0.25d0

  call build_points(pts, npt, p1, p2, p3)

  ! ---------------- accuracy pass ----------------
  open(unit=10, file='td_results.dat', status='replace')
  do i = 1, npt
     loc = pts(:,i)
     call tdstresshs(loc, p1, p2, p3, ss, ds, ts, stress, strain)
     call tddisphs  (loc, p1, p2, p3, ss, ds, ts, u)
     write(10,'(15es26.17)') stress, strain, u
  end do
  close(10)

  ! ---------------- timing pass: stress ----------------
  csum = 0.0d0
  call system_clock(c0, crate)
  do rep = 1, nrep
     do i = 1, npt
        loc = pts(:,i)
        call tdstresshs(loc, p1, p2, p3, ss, ds, ts, stress, strain)
        csum = csum + stress(1) + stress(4) + stress(6)
     end do
  end do
  call system_clock(c1)
  tstress = dble(c1 - c0) / dble(crate)

  ! ---------------- timing pass: displacement ----------------
  call system_clock(c0, crate)
  do rep = 1, nrep
     do i = 1, npt
        loc = pts(:,i)
        call tddisphs(loc, p1, p2, p3, ss, ds, ts, u)
        csum = csum + u(1) + u(3)
     end do
  end do
  call system_clock(c1)
  tdisp = dble(c1 - c0) / dble(crate)

  write(*,'(a,es24.16)') 'CHECKSUM ', csum
  write(*,'(a,f12.5)')   'STRESS_TIME_S ', tstress
  write(*,'(a,f12.5)')   'DISP_TIME_S ', tdisp
  write(*,'(a,i0)')      'CALLS_EACH ', npt*nrep

contains

  subroutine build_points(pts, n, p1, p2, p3)
    integer, intent(in)  :: n
    real*8, intent(out)  :: pts(3,n)
    real*8, intent(in)   :: p1(3), p2(3), p3(3)
    real*8 :: c(3), e1(3), e2(3), e3(3), r, ang
    integer :: i

    c  = (p1 + p2 + p3) / 3.0d0
    e1 = (p1 + p2) / 2.0d0
    e2 = (p2 + p3) / 2.0d0
    e3 = (p3 + p1) / 2.0d0

    ! points 1..32: representative, scattered on shells around the centroid,
    ! always kept below the free surface (z < 0)
    do i = 1, 32
       r   = 0.3d0 + 9.7d0 * dble(modulo(i*7, 17)) / 16.0d0
       ang = 6.2831853071795862d0 * dble(modulo(i*5, 13)) / 13.0d0
       pts(1,i) = c(1) + r * cos(ang)
       pts(2,i) = c(2) + r * sin(ang) * 0.8d0
       pts(3,i) = c(3) - r * 0.5d0 * abs(sin(ang*1.3d0)) - 0.05d0
    end do

    ! points 33..48: near-singular stress tests (small but nonzero offsets from
    ! vertices, edge midpoints, the free surface, and the TD plane)
    pts(:,33) = p1 + (/  1.0d-2,  1.0d-2, -1.0d-2 /)
    pts(:,34) = p1 + (/  1.0d-4, -1.0d-4, -1.0d-4 /)
    pts(:,35) = p2 + (/ -1.0d-3,  1.0d-3, -1.0d-3 /)
    pts(:,36) = p2 + (/  1.0d-6,  1.0d-6, -1.0d-6 /)
    pts(:,37) = p3 + (/  1.0d-3, -1.0d-3, -1.0d-3 /)
    pts(:,38) = p3 + (/ -1.0d-5,  1.0d-5, -1.0d-5 /)
    pts(:,39) = e1 + (/  0.0d0,   0.0d0,  -1.0d-3 /)
    pts(:,40) = e1 + (/  1.0d-5,  0.0d0,  -1.0d-5 /)
    pts(:,41) = e2 + (/  0.0d0,   1.0d-3, -1.0d-3 /)
    pts(:,42) = e2 + (/  0.0d0,   0.0d0,  -1.0d-6 /)
    pts(:,43) = e3 + (/ -1.0d-3,  0.0d0,  -1.0d-3 /)
    pts(:,44) = e3 + (/  0.0d0,  -1.0d-5, -1.0d-5 /)
    pts(:,45) = (/ c(1)+0.1d0, c(2)+0.1d0, 0.0d0 /)
    pts(:,46) = (/ c(1)-0.3d0, c(2)+0.2d0, 0.0d0 /)
    pts(:,47) = c + (/  1.0d-4,  1.0d-4,  1.0d-3 /)
    pts(:,48) = c + (/ -1.0d-4, -1.0d-4, -1.0d-3 /)
  end subroutine build_points

end program tdbench
