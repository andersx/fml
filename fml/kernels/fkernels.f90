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

    double precision :: maxgausdist
    double precision :: d

    double precision :: r_dist
    double precision :: c_dist

    integer :: m_1, m_2

    double precision :: maxgausdist2

    double precision :: inv_cut, inv_width, killer
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
    double precision, dimension(:), allocatable :: ltc

    ! Intermediary array
    double precision, dimension(:,:), allocatable :: ltcl

    ! Number of atomic representations (i.e. the total number of atoms)
    integer :: na

    ! Number of molecules
    integer :: nm

    ! Internal counters
    integer :: i, j, k, l
    integer :: k_start, k_end, l_start, l_end, idx
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

    double precision :: rdist, r_width2
    double precision :: cdist, c_width2
    double precision :: l2dist, inv_cut


    double precision, allocatable, dimension(:,:) :: norm

    ! Pre-computed sine terms
    double precision, allocatable, dimension(:,:,:) :: sin1

    ! Value of PI at full FORTRAN precision.
    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)

    ! Small number // to test for numerical stability
    double precision, parameter :: eps = 1.0e-12

    r_width2 = r_width**2
    c_width2 = c_width**2

    inv_cut = pi / (2.0d0 * cut_distance)

    nm = size(n, dim=1) ! Number of molecules
    na = sum(n(:nm)) ! Number of atomic representations

    inv_sigma2 = -1.0d0 / sigma**2

    allocate(sin1(nm, maxval(n), maxval(n)))
    allocate(norm(nm, maxval(n)))

!$OMP PARALLEL DO
    do i = 1, nm
        do m_1 = 1, n(i)
            do j_1 = 1, n(i)
                sin1(i, j_1, m_1) = 1.0d0 - sin(q(i,j_1,1,m_1) * inv_cut)
            enddo
        enddo
    enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
    do i = 1, nm
        do j_1 = 1, N(i)

            norm(i,j_1) = atomic_distl2(q(i,i_1,:,:), q(i,i_1,:,:), n(i), n(i), &
                & sin1(i,i_1,:), sin1(i,i_1,:), width, cut_distance, r_width, c_width)

        enddo
    enddo
!$OMP END PARALLEL DO

    allocate(k_starts(nm))
    allocate(k_ends(nm))

!$OMP PARALLEL DO
    do i = 1, nm
        k_ends(i) = sum(n(:i))
        k_starts(i) = k_ends(i) - n(i) + 1
    enddo
!$OMP END PARALLEL DO

    allocate(ltc(na))
    allocate(ltcl(nm,nm))

!$OMP PARALLEL DO PRIVATE(ltc,l_start,l_end,j_1,l2dist)

    ! Loop over MOL1
    do j = 1, nm

        ltc(:) = 0.0d0
        ! k_start = k_starts(j)
        ! k_end = k_ends(j)

        do i = 1, nm
            do i_1 = 1, n(i)
                idx = k_starts(i) + i_1 - 1

                do j_1 = 1, n(j)

                    l2dist = atomic_distl2(q(i,i_1,:,:), q(j,j_1,:,:), n(i), n(j), &
                        & sin1(i,i_1,:), sin1(j,j_1,:), width, cut_distance, r_width, c_width)

                    l2dist = norm(i,i_1) + norm(j,j_1) - 2.0d0 * l2dist &
                        & * (r_width2/(r_width2 + (z(i,i_1,1) - z(j,j_1,1))**2) * &
                        & c_width2/(c_width2 + (z(i,i_1,2) - z(j,j_1,2))**2))

                    if (abs(l2dist) < eps) l2dist = 0.0d0

                    ltc(idx) = ltc(idx) + exp(inv_sigma2 * l2dist)

                enddo
            enddo
        enddo

        do i = 1, j ! Only upper triangle is necessary for Cholesky decomposition
            ltcl(i,j) = sum(ltc(k_starts(i):k_ends(i)))
        enddo

    enddo
!$OMP END PARALLEL DO

    deallocate(norm)
    deallocate(sin1)

    deallocate(ltc)

!$OMP PARALLEL DO
    do i = 1, nm
        ltcl(i,i) = ltcl(i,i) + lambda
    enddo
!$OMP END PARALLEL DO

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

    alpha = 0.0d0
!$OMP PARALLEL DO PRIVATE(k_start,k_end,k)
    do j = 1, nm
        k_start = k_starts(j)
        k_end = k_ends(j)
        do k = k_start, k_end
            alpha(k) = ltcly(j)
        enddo
    enddo
!$OMP END PARALLEL DO

    deallocate(k_starts)
    deallocate(k_ends)
    deallocate(ltcly)

end subroutine fget_alpha_arad


subroutine fgaussian_kernel(a, na, b, nb, k, sigma)

    implicit none

    double precision, dimension(:,:), intent(in) :: a
    double precision, dimension(:,:), intent(in) :: b

    integer, intent(in) :: na, nb

    double precision, dimension(:,:), intent(inout) :: k
    double precision, intent(in) :: sigma

    double precision, allocatable, dimension(:) :: temp

    double precision :: inv_sigma
    integer :: i, j

    inv_sigma = -0.5d0 / (sigma*sigma)

    allocate(temp(size(a, dim=1)))

!$OMP PARALLEL DO PRIVATE(temp)
    do i = 1, nb
        do j = 1, na
            temp(:) = a(:,j) - b(:,i)
            k(j,i) = exp(inv_sigma * sqrt(sum(temp*temp)))
        enddo
    enddo
!$OMP END PARALLEL DO

    deallocate(temp)

end subroutine fgaussian_kernel

subroutine flaplacian_kernel(a, na, b, nb, k, sigma)

    implicit none

    double precision, dimension(:,:), intent(in) :: a
    double precision, dimension(:,:), intent(in) :: b

    integer, intent(in) :: na, nb

    double precision, dimension(:,:), intent(inout) :: k
    double precision, intent(in) :: sigma

    double precision :: inv_sigma

    integer :: i, j

    inv_sigma = -1.0d0 / sigma

!$OMP PARALLEL DO
    do i = 1, nb
        do j = 1, na
            k(j,i) = exp(inv_sigma * sum(abs(a(:,j) - b(:,i))))
        enddo
    enddo
!$OMP END PARALLEL DO

end subroutine flaplacian_kernel


subroutine fget_alpha(q, n, y, sigma, lambda, alpha)

    implicit none

    double precision, dimension(:,:), intent(in) :: q
    integer, dimension(:), intent(in) :: n
    double precision, dimension(:), intent(in) :: y
    double precision, intent(in) :: lambda
    double precision, intent(in) :: sigma

    double precision, dimension(:), intent(inout) :: alpha
    double precision, dimension(:), allocatable :: ltc
    double precision, dimension(:,:), allocatable :: ltcl

    integer :: na, nm

    integer :: i, j, k, k_start, k_end
    integer :: info
    double precision, allocatable, dimension(:) :: ltcly
    integer, allocatable, dimension(:) :: k_starts
    integer, allocatable, dimension(:) :: k_ends

    double precision :: inv_sigma

    na = size(q, dim=2) ! Number of atomic representations
    nm = size(n, dim=1) ! Number of molecules

    inv_sigma = -1.0d0 / sigma

    allocate(ltc(na))
    allocate(ltcl(nm,nm))

    allocate(k_starts(nm))
    allocate(k_ends(nm))

!$OMP PARALLEL DO
    do i = 1, nm
        k_ends(i) = sum(n(:i))
        k_starts(i) = k_ends(i) - n(i) + 1
    enddo
!$OMP END PARALLEL DO


!$OMP PARALLEL DO PRIVATE(ltc,k_start,k_end)
    do j = 1, nm
        ltc(:na) = 0.0d0
        k_start = k_starts(j)
        k_end = k_ends(j)
        do i = 1, na
            do k = k_start, k_end
                ltc(i) = ltc(i) + exp(inv_sigma * sum(abs(q(:,k) - q(:,i))))
            enddo
        enddo
        do i = 1, j ! Only upper triangle is necessary for Cholesky decomposition
            ltcl(i,j) = sum(ltc(k_starts(i):k_ends(i)))
        enddo
    enddo
!$OMP END PARALLEL DO

    deallocate(ltc)

!$OMP PARALLEL DO
    do i = 1, nm
        ltcl(i,i) = ltcl(i,i) + lambda
    enddo
!$OMP END PARALLEL DO

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

    alpha = 0.0d0
!$OMP PARALLEL DO PRIVATE(k_start,k_end,k)
    do j = 1, nm
        k_start = k_starts(j)
        k_end = k_ends(j)
        do k = k_start, k_end
            alpha(k) = ltcly(j)
        enddo
    enddo
!$OMP END PARALLEL DO

    deallocate(k_starts)
    deallocate(k_ends)
    deallocate(ltcly)

end subroutine fget_alpha

subroutine fget_prediction(Q, Q2, N2, alpha, sigma, Y)

    implicit none

    double precision, dimension(:,:), intent(in) :: Q
    double precision, dimension(:,:), intent(in) :: Q2
    integer, dimension(:), intent(in) :: N2
    double precision, dimension(:), intent(in) :: alpha
    double precision, intent(in) :: sigma

    double precision, dimension(:), intent(inout) :: Y

    double precision, allocatable, dimension(:) :: csa
    double precision :: inv_sigma

    integer, allocatable, dimension(:) :: k_starts
    integer, allocatable, dimension(:) :: k_ends

    integer :: i, j, na, na2, nm2

    na = size(Q, dim=2)
    na2 = size(Q2, dim=2)
    nm2 = size(N2, dim=1)

    inv_sigma = -1.0d0 / sigma

    allocate(k_starts(nm2))
    allocate(k_ends(nm2))

!$OMP PARALLEL DO
    do i = 1, nm2
        k_ends(i) = sum(n2(:i))
        k_starts(i) = k_ends(i) - n2(i) + 1
    enddo
!$OMP END PARALLEL DO

    allocate(csa(na2))
    csa(:) = 0.0d0

!$OMP PARALLEL DO
    do i = 1, na2
        do j = 1, na
            csa(i) = csa(i) + exp(inv_sigma * sum(abs(q(:,j) - q2(:,i))))*alpha(j)
        enddo
    enddo
!$OMP END PARALLEL DO


!$OMP PARALLEL DO
    do j = 1, nm2
        Y(j) = sum(csa(k_starts(j):k_ends(j)))
    enddo
!$OMP END PARALLEL DO

    deallocate(csa)
    deallocate(k_starts)
    deallocate(k_ends)


end subroutine fget_prediction

subroutine fget_alpha_from_distance(d, n, y, dgamma, lambda, alpha)

    implicit none

    double precision, dimension(:,:), intent(in) :: d
    integer, dimension(:), intent(in) :: n
    double precision, dimension(:), intent(in) :: y
    double precision, intent(in) :: dgamma
    double precision, intent(in) :: lambda

    double precision, dimension(:), intent(inout) :: alpha
    double precision, dimension(:,:), allocatable :: ltc
    double precision, dimension(:,:), allocatable :: ltcl

    integer :: na, nm

    integer :: i, j, k
    integer :: info
    double precision, allocatable, dimension(:) :: ltcly
    integer, allocatable, dimension(:) :: k_starts
    integer, allocatable, dimension(:) :: k_ends

    double precision :: inv_gamma

    na = size(d, dim=2) ! Number of atomic representations
    nm = size(n, dim=1) ! Number of molecules

    inv_gamma = -1.0d0 / dgamma

    allocate(k_starts(nm))
    allocate(k_ends(nm))

!$OMP PARALLEL DO
    do i = 1, nm
        k_ends(i) = sum(n(:i))
        k_starts(i) = k_ends(i) - n(i) + 1
    enddo
!$OMP END PARALLEL DO

    allocate(ltc(nm,na))
    allocate(ltcl(nm,nm))

!$OMP PARALLEL DO
    do i = 1, na

        ! READ IN SLICES HERE
        do j = 1, nm
            ltc(j,i) = 0.0d0
            do k = k_starts(j), k_ends(j)
                ltc(j,i) = ltc(j,i) + exp(inv_gamma * d(k,i))
            enddo
        enddo
    enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
    do j = 1, nm
        do i = 1, nm
            ltcl(i,j) = sum(ltc(j,k_starts(i):k_ends(i)))
        enddo
    enddo
!$OMP END PARALLEL DO

    deallocate(ltc)

!$OMP PARALLEL DO
    do i = 1, nm
        ltcl(i,i) = ltcl(i,i) + lambda
    enddo
!$OMP END PARALLEL DO

    call dpotrf("U", nm, ltcl, nm, info)
    if (info > 0) then
        write (*,*) "WARNING: Cholesky decomposition DPOTRF() exited with error code:", info
    endif
    call dpotri("U", nm, ltcl, nm, info )
    if (info > 0) then
        write (*,*) "WARNING: Cholesky inversion DPOTRI() exited with error code:", info
    endif

    allocate(ltcly(nm))
    ltcly = 0.0d0

    call dsymv("U", nm, 1.0d0, ltcl, nm, y, 1, 0.0d0, ltcly, 1)

    deallocate(ltcl)

    alpha = 0.0d0
!$OMP PARALLEL DO
    do j = 1, nm
        do k = k_starts(j), k_ends(j)
            alpha(k) = ltcly(j)
        enddo
    enddo
!$OMP END PARALLEL DO

    deallocate(k_starts)
    deallocate(k_ends)
    deallocate(ltcly)

end subroutine fget_alpha_from_distance



subroutine fget_prediction_from_distance(D2, N2, alpha, dgamma, Y)

    implicit none

    double precision, dimension(:,:), intent(in) :: D2
    nteger, dimension(:), intent(in) :: N2
    double precision, dimension(:), intent(in) :: alpha
    double precision, intent(in) :: dgamma

    double precision, dimension(:), intent(inout) :: Y

    double precision, allocatable, dimension(:) :: csa
    double precision :: inv_gamma

    integer, allocatable, dimension(:) :: k_starts
    integer, allocatable, dimension(:) :: k_ends

    integer :: i, j, na, na2, nm2

    na = size(D2, dim=1)
    na2 = size(D2, dim=2)
    nm2 = size(N2, dim=1)

    inv_gamma = -1.0d0 / dgamma

    allocate(k_starts(nm2))
    allocate(k_ends(nm2))

!$OMP PARALLEL DO
    do i = 1, nm2
        k_ends(i) = sum(n2(:i))
        k_starts(i) = k_ends(i) - n2(i) + 1
    enddo
!$OMP END PARALLEL DO

    allocate(csa(na2))

!$OMP PARALLEL DO
    do i = 1, na2
        do j = 1, na
            csa(i) = csa(i) + exp(inv_gamma * D2(j,i))*alpha(j)
        enddo
    enddo
!$OMP END PARALLEL DO


!$OMP PARALLEL DO
    do j = 1, nm2
        Y(j) = sum(csa(k_starts(j):k_ends(j)))
    enddo
!$OMP END PARALLEL DO

    deallocate(csa)
    deallocate(k_starts)
    deallocate(k_ends)


end subroutine fget_prediction_from_distance

subroutine fmanhattan_distance(a, na, b, nb, k)

    implicit none

    double precision, dimension(:,:), intent(in) :: a
    double precision, dimension(:,:), intent(in) :: b

    integer, intent(in) :: na, nb

    double precision, dimension(:,:), intent(inout) :: k

    integer :: i, j

!$OMP PARALLEL DO
    do i = 1, nb
        do j = 1, na
            k(j,i) = sum(abs(a(:,j) - b(:,i)))
        enddo
    enddo
!$OMP END PARALLEL DO

end subroutine fmanhattan_distance
