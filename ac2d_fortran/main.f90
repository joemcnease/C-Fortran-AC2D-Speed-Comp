program main
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    integer, parameter  :: nx = 301
    integer, parameter  :: nz = 201
    integer, parameter  :: sx = 150
    integer, parameter  :: sz = 100
    integer, parameter  :: nt = 1000

    real(dp), parameter :: dx = 5.0_dp
    real(dp), parameter :: dz = 5.0_dp
    real(dp), parameter :: dt = 0.0005_dp
    real(dp), parameter :: t0 = 0.01_dp
    real(dp), parameter :: f0 = 100.0_dp
    real(dp), parameter :: c0 = 3000.0_dp
    real(dp), parameter :: amp = 1.0_dp

    real(dp), dimension(nt, nz, nx) :: p(0:nt-1, 0:nz-1, 0:nx-1)
    real(dp), dimension(nz, nx)     :: c(0:nz-1, 0:nx-1)
    real(dp), dimension(nt)         :: stf(0:nt-1)

    real(dp) :: t1, t2


    call system('mkdir -p pressure/')


    p(:, :, :) = 0.0_dp
    c(:, :) = c0
    stf(:) = 0.0_dp


    call gaussian(stf, nt, dt, t0, f0, amp)


    call cpu_time(t1)
    call ac2d(p, c, stf, nx, nz, dx, dz, nt, dt, sx, sz)
    call cpu_time(t2)

    print *, "Time to compute: ", t2-t1


    contains

    subroutine ac2d(p, c, stf, nx, nz, dx, dz, nt, dt, sx, sz)
        implicit none

        real(dp), dimension(nt, nz, nx), intent(inout) :: p(0:nt-1, 0:nz-1, 0:nx-1)
        real(dp), dimension(nz, nx), intent(in)        :: c(0:nz-1, 0:nx-1)
        real(dp), dimension(nt), intent(in)            :: stf(0:nt-1)

        real(dp), intent(in) :: dx
        real(dp), intent(in) :: dz
        real(dp), intent(in) :: dt

        integer, intent(in) :: nx
        integer, intent(in) :: nz
        integer, intent(in) :: nt
        integer, intent(in) :: sx
        integer, intent(in) :: sz

        real(dp), dimension(nz, nx) :: d2px(0:nz-1, 0:nx-1)
        real(dp), dimension(nz, nx) :: d2pz(0:nz-1, 0:nx-1)
        real(dp), dimension(nz, nx) :: pNew(0:nz-1, 0:nx-1)
        real(dp), dimension(nz, nx) :: pOld(0:nz-1, 0:nx-1)

        character(len=50) :: fname
        character(len=10) :: fidx

        integer :: io

        integer :: i, j


        d2px(:, :) = 0.0_dp
        d2pz(:, :) = 0.0_dp
        pNew(:, :) = 0.0_dp
        pOld(:, :) = 0.0_dp


        do i=0, nt-1
            do j=1, nz-2
                d2pz(j, :) = (p(i, j-1, :) - 2*p(i, j, :) + &
                    p(i, j+1, :)) / dz**2
            end do

            do j=1, nx-2
                d2px(:, j) = (p(i, :, j-1) - 2*p(i, :, j) + &
                    p(i, :, j+1)) / dx**2
            end do

            pNew(:, :) = 2.*p(i, :, :) - pOld(:, :) + &
                (dt**2)*(c(:, :)**2)*(d2px(:, :) + d2pz(:, :))

            pNew(sz, sx) = pNew(sz, sx) + stf(i)
            pOld(:, :) = p(i, :, :)
            p(i+1, :, :) = pNew(:, :)

            if (mod(i, 10) == 0) then
                write(fidx, '(I0.3)') i
                fname = "pressure/pressure"//trim(adjustl(fidx))//".csv"
                open(newunit=io, file=fname, action='write', form='unformatted')
                write(io) p(i+1, :, :)
                close(io)
            end if

            print *, "Time step: ", i
        end do


    end subroutine ac2d


    subroutine gaussian(stf, nt, dt, t0, f0, amp)
        implicit none

        real(dp), dimension(nt), intent(inout) :: stf(0:nt-1)

        integer, intent(in) :: nt

        real(dp), intent(in) :: dt
        real(dp), intent(in) :: t0
        real(dp), intent(in) :: f0
        real(dp), intent(in) :: amp

        integer :: i


        do i=0, nt-1
            stf(i) = exp(-(f0**2)*((i*dt - t0)**2))
        end do


    end subroutine gaussian


end program main
