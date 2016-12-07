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
    double precision, dimension(:), intent(in) :: n
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
    double precision, dimension(:), intent(in) :: N2
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
    double precision, dimension(:), intent(in) :: n
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
    double precision, dimension(:), intent(in) :: N2
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
