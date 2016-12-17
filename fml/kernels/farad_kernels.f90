! MIT License
!
! Copyright (c) 2016 Anders Steen Christensen
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module arad

    implicit none

contains

function atomic_distl2(X1, X2, N1, N2, sin1, sin2, width, cut_distance, r_width, c_width) result(aadist)

    implicit none

    double precision, dimension(:,:), intent(in) :: X1
    double precision, dimension(:,:), intent(in) :: X2

    integer, intent(in) :: N1
    integer, intent(in) :: N2

    double precision, dimension(:), intent(in) :: sin1
    double precision, dimension(:), intent(in) :: sin2

    double precision, intent(in) :: width
    double precision, intent(in) :: cut_distance
    double precision, intent(in) :: r_width
    double precision, intent(in) :: c_width

    double precision :: aadist

    double precision :: d

    integer :: m_1, m_2

    double precision :: maxgausdist2

    double precision :: inv_width
    double precision :: c_width2, r_width2, r2

    inv_width = -1.0d0 / (4.0d0 * width**2)

    maxgausdist2 = (8.0d0 * width)**2
    r_width2 = r_width**2
    c_width2 = c_width**2

    aadist = 0.0d0

    do m_1 = 1, N1

        if (X1(1, m_1) > cut_distance) exit

        do m_2 = 1, N2

            if (X2(1, m_2) > cut_distance) exit

            r2 = (X2(1,m_2) - X1(1,m_1))**2

            if (r2 < maxgausdist2) then

                d = exp(r2 * inv_width )  * sin1(m_1) * sin2(m_2)

                d = d * (r_width2/(r_width2 + (x1(2,m_1) - x2(2,m_2))**2) * &
                    & c_width2/(c_width2 + (x1(3,m_1) - x2(3,m_2))**2))

                aadist = aadist + d * (1.0d0 + x1(4,m_1)*x2(4,m_2) + x1(5,m_1)*x2(5,m_2))

            end if
        end do
    end do

end function atomic_distl2

end module arad


subroutine fget_kernel_arad(q, z, n, sigma, nm, kernel)

    use arad, only: atomic_distl2

    implicit none

    double precision, parameter :: width = 0.2d0
    double precision, parameter :: cut_distance = 6.0d0
    double precision, parameter :: r_width = 1.0d0
    double precision, parameter :: c_width = 0.5d0

    ! ARAD descriptors for the training set, format (i,j_1,5,m_1)
    double precision, dimension(:,:,:,:), intent(in) :: q

    ! ARAD atom-types for each atom in each molecule, format (i, j_1, 2)
    double precision, dimension(:,:,:), intent(in) :: z

    ! List of numbers of atoms in each molecule
    integer, dimension(:), intent(in) :: n

    ! Sigma in the Gaussian kernel
    double precision, intent(in) :: sigma

    ! Resulting alpha vector
    double precision, dimension(nm,nm), intent(out) :: kernel

    ! Number of atomic representations (i.e. the total number of atoms)
    integer :: na

    ! Number of molecules
    integer :: nm

    ! Internal counters
    integer :: i, j, ni,nj
    integer :: m_1, i_1, j_1

    ! -1.0 / sigma^2 for use in the kernel
    double precision :: inv_sigma2

    double precision :: r_width2
    double precision :: c_width2
    double precision :: l2dist, inv_cut
    double precision, allocatable, dimension(:,:) :: atomic_kernel

    double precision, allocatable, dimension(:,:) :: norm

    ! Pre-computed sine terms
    double precision, allocatable, dimension(:,:,:) :: sin1

    ! Value of PI at full FORTRAN precision.
    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)

    ! Small number // to test for numerical stability
    double precision, parameter :: eps = 1.0e-11

    r_width2 = r_width**2
    c_width2 = c_width**2

    inv_cut = pi / (2.0d0 * cut_distance)

    ! nm = size(n, dim=1) ! Number of molecules
    na = sum(n(:nm)) ! Number of atomic representations

    inv_sigma2 = -1.0d0 / sigma**2

    allocate(sin1(nm, maxval(n), maxval(n)))
    allocate(norm(nm, maxval(n)))

!$OMP PARALLEL DO PRIVATE(ni)
    do i = 1, nm
        ni = n(i)
        do m_1 = 1, ni
            do i_1 = 1, ni
                sin1(i, i_1, m_1) = 1.0d0 - sin(q(i,i_1,1,m_1) * inv_cut)
            enddo
        enddo
    enddo
!$OMP END PARALLEL DO



!$OMP PARALLEL DO PRIVATE(ni)
    do i = 1, nm
        ni = n(i)
        do i_1 = 1, ni
            norm(i,i_1) = atomic_distl2(q(i,i_1,:,:), q(i,i_1,:,:), n(i), n(i), &
                & sin1(i,i_1,:), sin1(i,i_1,:), width, cut_distance, r_width, c_width)
        enddo
    enddo
!$OMP END PARALLEL DO

    allocate(atomic_kernel(maxval(n), maxval(n)))

    kernel(:,:) = 0.0d0
    atomic_kernel(:,:) = 0.0d0

!$OMP PARALLEL DO PRIVATE(l2dist,atomic_kernel,ni,nj)
    do j = 1, nm
        nj = n(j)
        do i = 1, j ! Only upper triangle is necessary for Cholesky decomposition
            ni = n(i)

            atomic_kernel(:,:) = 0.0d0

            do i_1 = 1, ni
                do j_1 = 1, nj

                    l2dist = atomic_distl2(q(i,i_1,:,:), q(j,j_1,:,:), n(i), n(j), &
                        & sin1(i,i_1,:), sin1(j,j_1,:), width, cut_distance, r_width, c_width)

                    l2dist = norm(i,i_1) + norm(j,j_1) - 2.0d0 * l2dist &
                        & * (r_width2/(r_width2 + (z(i,i_1,1) - z(j,j_1,1))**2) * &
                        & c_width2/(c_width2 + (z(i,i_1,2) - z(j,j_1,2))**2))

                    if (abs(l2dist) < eps) l2dist = 0.0d0

                    atomic_kernel(i_1,j_1) =  exp(inv_sigma2 * l2dist)

                enddo
            enddo

            kernel(i,j) =  sum(atomic_kernel(:,:))
            kernel(j,i) =  sum(atomic_kernel(:,:))

        enddo
    enddo
!$OMP END PARALLEL DO

    deallocate(atomic_kernel)
    deallocate(norm)
    deallocate(sin1)

end subroutine fget_kernel_arad



subroutine fget_alpha_arad(q, z, n, y, sigma, lambda, alpha)

    use arad, only: atomic_distl2

    implicit none

    double precision, parameter :: width = 0.2d0
    double precision, parameter :: cut_distance = 6.0d0
    double precision, parameter :: r_width = 1.0d0
    double precision, parameter :: c_width = 0.5d0

    ! ARAD descriptors for the training set, format (i,j_1,5,m_1)
    double precision, dimension(:,:,:,:), intent(in) :: q

    ! ARAD atom-types for each atom in each molecule, format (i, j_1, 2)
    double precision, dimension(:,:,:), intent(in) :: z

    ! List of numbers of atoms in each molecule
    integer, dimension(:), intent(in) :: n

    ! List of learning property for each molecule
    double precision, dimension(:), intent(in) :: y

    ! Lambda for kernel diagonal
    double precision, intent(in) :: lambda

    ! Sigma in the Gaussian kernel
    double precision, intent(in) :: sigma

    ! Resulting alpha vector
    double precision, dimension(:), intent(inout) :: alpha

    ! Intermediary array
    ! double precision, dimension(:), allocatable :: ltc

    ! Intermediary array
    double precision, dimension(:,:), allocatable :: ltcl

    ! Number of atomic representations (i.e. the total number of atoms)
    integer :: na

    ! Number of molecules
    integer :: nm

    ! Internal counters
    integer :: i, j, k, ni, nj
    integer :: k_start, k_end
    integer :: m_1, i_1, j_1

    ! Info for LAPACK
    integer :: info

    ! Intermediary array
    double precision, allocatable, dimension(:) :: ltcly

    ! List of indexes
    integer, allocatable, dimension(:) :: k_starts

    ! List of indexes
    integer, allocatable, dimension(:) :: k_ends

    ! -1.0 / sigma^2 for use in the kernel
    double precision :: inv_sigma2

    double precision :: r_width2
    double precision :: c_width2
    double precision :: l2dist, inv_cut


    double precision, allocatable, dimension(:,:) :: atomic_kernel
    double precision, allocatable, dimension(:,:) :: norm

    ! Pre-computed sine terms
    double precision, allocatable, dimension(:,:,:) :: sin1

    ! Value of PI at full FORTRAN precision.
    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)

    ! Small number // to test for numerical stability
    double precision, parameter :: eps = 1.0e-11

    r_width2 = r_width**2
    c_width2 = c_width**2

    inv_cut = pi / (2.0d0 * cut_distance)

    nm = size(n, dim=1) ! Number of molecules
    na = sum(n(:nm)) ! Number of atomic representations

    inv_sigma2 = -1.0d0 / sigma**2

    allocate(sin1(nm, maxval(n), maxval(n)))
    allocate(norm(nm, maxval(n)))

!$OMP PARALLEL DO PRIVATE(ni)
    do i = 1, nm
        ni = n(i)
        do m_1 = 1, ni
            do i_1 = 1, ni
                sin1(i, i_1, m_1) = 1.0d0 - sin(q(i,i_1,1,m_1) * inv_cut)
            enddo
        enddo
    enddo
!$OMP END PARALLEL DO



!$OMP PARALLEL DO PRIVATE(ni)
    do i = 1, nm
        ni = n(i)
        do i_1 = 1, ni
            norm(i,i_1) = atomic_distl2(q(i,i_1,:,:), q(i,i_1,:,:), n(i), n(i), &
                & sin1(i,i_1,:), sin1(i,i_1,:), width, cut_distance, r_width, c_width)
        enddo
    enddo
!$OMP END PARALLEL DO

    allocate(atomic_kernel(maxval(n), maxval(n)))
    allocate(ltcl(nm,nm))

    ltcl(:,:) = 0.0d0

!$OMP PARALLEL DO PRIVATE(l2dist,atomic_kernel,ni,nj)
    do j = 1, nm
        nj = n(j)
        do i = 1, j ! Only upper triangle is necessary for Cholesky decomposition
            ni = n(i)

            atomic_kernel = 0.0d0

            do i_1 = 1, ni
                do j_1 = 1, nj

                    l2dist = atomic_distl2(q(i,i_1,:,:), q(j,j_1,:,:), n(i), n(j), &
                        & sin1(i,i_1,:), sin1(j,j_1,:), width, cut_distance, r_width, c_width)

                    l2dist = norm(i,i_1) + norm(j,j_1) - 2.0d0 * l2dist &
                        & * (r_width2/(r_width2 + (z(i,i_1,1) - z(j,j_1,1))**2) * &
                        & c_width2/(c_width2 + (z(i,i_1,2) - z(j,j_1,2))**2))

                    if (abs(l2dist) < eps) l2dist = 0.0d0

                    atomic_kernel(i_1,j_1) =  exp(inv_sigma2 * l2dist)

                enddo
            enddo

            ltcl(i,j) =  sum(atomic_kernel(:,:))
            ltcl(j,i) =  sum(atomic_kernel(:,:))

        enddo
    enddo
!$OMP END PARALLEL DO

    deallocate(atomic_kernel)
    deallocate(norm)
    deallocate(sin1)

! !$OMP PARALLEL DO
    do i = 1, nm
        ltcl(i,i) = ltcl(i,i) + lambda
    enddo
! !$OMP END PARALLEL DO

    call dpotrf("U", nm, ltcl, nm, info)
    if (info > 0) then
        write (*,*) "WARNING: Cholesky decomposition DPOTRF() exited with error code:", info
    endif
    call dpotri("U", nm, ltcl, nm, info )
    if (info > 0) then
        write (*,*) "WARNING: Cholesky inversion DPOTRI() exited with error code:", info
    endif

    allocate(ltcly(nm))
    ltcly(:) = 0.0d0

    call dsymv("U", nm, 1.0d0, ltcl, nm, y, 1, 0.0d0, ltcly, 1)

    deallocate(ltcl)

    allocate(k_starts(nm))
    allocate(k_ends(nm))

! !$OMP PARALLEL DO
    do i = 1, nm
        k_ends(i) = sum(n(:i))
        k_starts(i) = k_ends(i) - n(i) + 1
    enddo
! !$OMP END PARALLEL DO

    alpha = 0.0d0
! !$OMP PARALLEL DO PRIVATE(k_start,k_end,k)
    do j = 1, nm
        k_start = k_starts(j)
        k_end = k_ends(j)
        do k = k_start, k_end
            alpha(k) = ltcly(j)
        enddo
    enddo
! !$OMP END PARALLEL DO

    deallocate(k_starts)
    deallocate(k_ends)
    deallocate(ltcly)

end subroutine fget_alpha_arad
