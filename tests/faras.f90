module aras_utils
    implicit none
contains

function l2_distance(x1, x2, n1, n2, q1, q2, qall1, qall2, ksi1, ksi2, pd) result(aa)

    implicit none

    double precision, dimension(:,:), intent(in) :: x1
    double precision, dimension(:,:), intent(in) :: x2

    integer, intent(in) :: n1
    integer, intent(in) :: n2

    integer, intent(in) :: q1
    integer, intent(in) :: q2

    integer, dimension(:), intent(in) :: qall1
    integer, dimension(:), intent(in) :: qall2

    double precision, dimension(:), intent(in) :: ksi1 
    double precision, dimension(:), intent(in) :: ksi2

    double precision, dimension(:,:), intent(in) :: pd

    double precision :: a, b, c, d, e
    integer :: i, j, k, l

    double precision :: aa
    double precision :: ddist

    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)
    double precision, parameter :: t_width = pi / 4.0d0 ! ~0.7853981633974483
    double precision, parameter :: d_width = 0.2d0
    double precision, parameter :: cut_distance = 6.0d0
    double precision, parameter :: r_width = 1.0d0
    double precision, parameter :: c_width = 0.5d0

    double precision, parameter :: inv_cut = pi / (2.0d0 * cut_distance)
    double precision, parameter :: sqrt_pi = sqrt(pi)

    double precision, parameter :: t_fast = -1.0d0 / (4.0d0 * t_width**2)
    double precision, parameter :: d_fast = -1.0d0 / (4.0d0 * d_width**2)
    double precision :: maxgausdist2 = (8.0d0 * d_width)**2

    a = 1.0d0 / (sum(ksi1(:n1)) * sum(ksi2(:n2))) &
        &  * sqrt_pi * d_width * r_width**4 * c_width**4 &
        & * pd(q1, q2)
    
    b = 0.0d0

    do i = 1,  n1
        c = 0.0d0
        if (x1(1,i) > cut_distance) exit

        do j = 1, n2
            d = 0.0d0

            if (x2(1,j) > cut_distance) exit
            ddist = (x1(1,i) - x2(1,j))**2

            if (ddist < maxgausdist2) then

                do k = 1, n1
                    e = 0.0d0

                    do l = 1, n2
                         e = e + ksi2(l) * exp((x1(i+3,k) - x2(j+3,l))**2 * t_fast)

                    enddo
                    d = d + e * ksi1(k)

                enddo
                c = c + d * exp(ddist * d_fast) * ksi2(j) &
                    & * pd(qall1(i), qall2(j))
            endif

        enddo
        b = b + c * ksi1(i)

    enddo

    aa = a * b
    
end function l2_distance

end module aras_utils


! subroutine faras_molecular_distance(x1, x2, q1, q2, n1, n2, nm1, nm2, amax, pd, d)
! 
! 
!     use aras_utils, only: l2_distance
!     implicit none
! 
!     double precision, dimension(:,:,:,:), intent(in) :: x1
!     double precision, dimension(:,:,:,:), intent(in) :: x2
!     integer, dimension(:,:), intent(in) :: q1
!     integer, dimension(:,:), intent(in) :: q2
!     integer, dimension(:), intent(in) :: n1
!     integer, dimension(:), intent(in) :: n2
!     integer, intent(in) :: nm1
!     integer, intent(in) :: nm2
!     integer, intent(in) :: amax
!     double precision, dimension(:,:), intent(in) :: pd
!     double precision, dimension(:,:,:,:), intent(inout) :: d
! 
!     integer :: a, b, i, j, ii, jj
!     integer :: ni, nj
!     
!     double precision, dimension(nm1,amax,amax) :: ksi1
!     double precision, dimension(nm2,amax,amax) :: ksi2
! 
!     double precision, parameter :: pi = 4.0d0 * atan(1.0d0)
!     double precision, parameter :: cut_distance = 6.0d0
! 
!     double precision, parameter :: inv_cut = pi / (2.0d0 * cut_distance)
! 
!     double precision, dimension(nm1,amax) :: aa
!     double precision, dimension(nm2,amax) :: bb
! 
!     ksi1 = 0.0d0
!     ksi2 = 0.0d0
! 
! !$OMP PARALLEL DO PRIVATE(ni)
!     do a = 1, nm1
!         ni = n1(a)
!         do i = 1, ni
!             do ii = 1, ni
!                 ksi1(a,i,ii) = 1.0d0 - sin(x1(a,i,1,ii) * inv_cut)
!             enddo
!         enddo
!     enddo
! !$OMP END PARALLEL DO
! 
! !$OMP PARALLEL DO PRIVATE(nj)
!     do b = 1, nm2
!         nj = n2(b)
!         do j = 1, nj
!             do jj = 1, nj
!                 ksi2(b,j,jj) = 1.0d0 - sin(x2(b,j,1,jj) * inv_cut)
!             enddo
!         enddo
!     enddo
! !$OMP END PARALLEL DO
! 
!     aa = 0.0d0
!     bb = 0.0d0
! 
! !$OMP PARALLEL DO PRIVATE(ni)
!     do a = 1, nm1
!         ni = n1(a)
!         do i = 1, ni
!             aa(a,i) = l2_distance(x1(a,i,:,:), x1(a,i,:,:), n1(a), n1(a), q1(a,i), q1(a,i), q1(a,:), q1(a,:), &
!                                  & ksi1(a,i,:), ksi1(a,i,:), pd)
!         enddo
!     enddo
! !$OMP END PARALLEL DO
! 
! !$OMP PARALLEL DO PRIVATE(nj)
!     do b = 1, nm2
!         nj = n2(b)
!         do j = 1, nj
!             bb(b,j) = l2_distance(x2(b,j,:,:), x2(b,j,:,:), n2(b), n2(b), q2(b,j), q2(b,j), q2(b,:), q2(b,:), &
!                                  & ksi2(b,j,:), ksi2(b,j,:), pd)
!         enddo
!     enddo
! !$OMP END PARALLEL DO
! 
! !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ni,nj) SHARED(d)
!     do a = 1, nm1
!     ni = n1(a)
!     do b = 1, nm2
!     nj = n2(b)
!         do i = 1, ni
!         do j = 1, nj
! 
!             if (1>1) then
!                 ii = 1
!             endif        
!              d(a,b,i,j) = aa(a,i) + bb(b,j) &
!                 & - 2.0d0 * l2_distance(x1(a,i,:,:), x2(b,j,:,:), n1(a), n2(b), q1(a,i), q2(b,j), q1(a,:), q2(b,:), &
!                                  & ksi1(a,i,:), ksi2(b,j,:), pd)
!         enddo
!         enddo
! 
!     enddo
!     enddo
! !$OMP END PARALLEL DO
! 
! 
! end subroutine faras_molecular_distance

subroutine faras_molecular_distance(x1, x2, q1, q2, n1, n2, nm1, nm2, amax, pd, d)

    use aras_utils, only: l2_distance
    implicit none

    double precision, dimension(:,:,:,:), intent(in) :: x1
    double precision, dimension(:,:,:,:), intent(in) :: x2
    integer, dimension(:,:), intent(in) :: q1
    integer, dimension(:,:), intent(in) :: q2
    integer, dimension(:), intent(in) :: n1
    integer, dimension(:), intent(in) :: n2
    integer, intent(in) :: nm1
    integer, intent(in) :: nm2
    integer, intent(in) :: amax
    double precision, dimension(:,:), intent(in) :: pd
    double precision, dimension(nm1,nm2,amax,amax), intent(out) :: d

    integer :: a, b, i, j, ii, jj
    integer :: ni, nj
    
    double precision, allocatable, dimension(:,:,:) :: ksi1
    double precision, allocatable, dimension(:,:,:) :: ksi2

    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)
    double precision, parameter :: cut_distance = 6.0d0

    double precision, parameter :: inv_cut = pi / (2.0d0 * cut_distance)

    double precision, allocatable, dimension(:,:) :: aa
    double precision, allocatable, dimension(:,:) :: bb

    allocate(ksi1(nm1, maxval(n1), maxval(n1)))
    allocate(ksi2(nm2, maxval(n2), maxval(n2)))
    allocate(aa(nm1,maxval(n1)))
    allocate(bb(nm2,maxval(n2)))

    ksi1(:,:,:) = 0.0d0
    ksi2(:,:,:) = 0.0d0

!$OMP PARALLEL DO PRIVATE(ni)
    do a = 1, nm1
        ni = n1(a)
        do i = 1, ni
            do ii = 1, ni
                ksi1(a,i,ii) = 1.0d0 - sin(x1(a,i,1,ii) * inv_cut)
            enddo
        enddo
    enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(nj)
    do b = 1, nm2
        nj = n2(b)
        do j = 1, nj
            do jj = 1, nj
                ksi2(b,j,jj) = 1.0d0 - sin(x2(b,j,1,jj) * inv_cut)
            enddo
        enddo
    enddo
!$OMP END PARALLEL DO

    aa(:,:) = 0.0d0
    bb(:,:) = 0.0d0

!$OMP PARALLEL DO PRIVATE(ni)
    do a = 1, nm1
        ni = n1(a)
        do i = 1, ni
            aa(a,i) = l2_distance(x1(a,i,:,:), x1(a,i,:,:), n1(a), n1(a), q1(a,i), q1(a,i), q1(a,:), q1(a,:), &
                                 & ksi1(a,i,:), ksi1(a,i,:), pd)
        enddo
    enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(nj)
    do b = 1, nm2
        nj = n2(b)
        do j = 1, nj
            bb(b,j) = l2_distance(x2(b,j,:,:), x2(b,j,:,:), n2(b), n2(b), q2(b,j), q2(b,j), q2(b,:), q2(b,:), &
                                 & ksi2(b,j,:), ksi2(b,j,:), pd)
        enddo
    enddo
!$OMP END PARALLEL DO

    d(:,:,:,:) = 0.0d0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ni,nj) SHARED(d)
    do a = 1, nm1
    ni = n1(a)
    do b = 1, nm2
    nj = n2(b)
        do i = 1, ni
        do j = 1, nj

             d(a,b,i,j) = aa(a,i) + bb(b,j) &
                & - 2.0d0 * l2_distance(x1(a,i,:,:), x2(b,j,:,:), n1(a), n2(b), q1(a,i), q2(b,j), q1(a,:), q2(b,:), &
                                 & ksi1(a,i,:), ksi2(b,j,:), pd)
        enddo
        enddo

    enddo
    enddo
!$OMP END PARALLEL DO

    deallocate(ksi1)
    deallocate(ksi2)
    deallocate(aa)
    deallocate(bb)

end subroutine faras_molecular_distance
