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
    double precision, allocatable, dimension(:) :: temp

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


subroutine fget_alpha(q, n, y, sigma, alpha)

    implicit none

    double precision, dimension(:,:), intent(in) :: q
    double precision, dimension(:), intent(in) :: n
    double precision, dimension(:), intent(in) :: y
    double precision, intent(in) :: sigma

    double precision, dimension(:), intent(inout) :: alpha
    double precision, dimension(:,:), allocatable :: ltc
    double precision, dimension(:,:), allocatable :: ltcl

    integer :: na, nm

    integer :: i, j, k, k_start, k_end
    integer :: lwork, info
    integer, allocatable, dimension(:) :: ipiv
    integer, allocatable, dimension(:) :: work
    double precision, allocatable, dimension(:) :: ltcly

    double precision :: inv_sigma
    double precision :: temp

    na = size(q, dim=2) ! Number of atomic representations
    nm = size(n, dim=1) ! Number of molecules

    inv_sigma = -1.0d0 / sigma

    allocate(ltc(nm,na))
    allocate(ltcl(nm,nm))

!$OMP PARALLEL DO PRIVATE(temp,k_start,k_end,k)
    do j = 1, nm
        do i = 1, na    
            k_start = sum(n(:j)) - n(j) + 1
            k_end = sum(n(:j))
            ! write (*,*) "start, end", k_start, k_end
            temp = 0.0d0
            do k = k_start, k_end
                temp = temp + exp(inv_sigma * sum(abs(q(:,k) - q(:,i))))
            enddo
            ltc(j,i) = temp
            ! write (*,*) i,j, ltc(j,i)
        enddo
    enddo
!$OMP END PARALLEL DO
    
!$OMP PARALLEL DO PRIVATE(temp,k_start,k_end,k)
    do i = 1, nm
        do j = 1, nm
            k_start = sum(n(:j)) - n(j) + 1
            k_end = sum(n(:j))
            temp = 0.0d0
            do k = k_start, k_end
                temp = temp + ltc(i,k)
            enddo
            ltcl(j,i) = temp
        enddo
    enddo
!$OMP END PARALLEL DO

    deallocate(ltc)

!$OMP PARALLEL DO
    do i = 1, nm
        ltcl(i,i) = ltcl(i,i) + 1.0d-10
    enddo 
!$OMP END PARALLEL DO

    allocate(ipiv(nm))

    call dgetrf(nm, nm, ltcl, nm, ipiv, info)
    if (info > 0) then
        write (*,*) "WARNING: L -decomposition DGETRF() exited with error code:", info
    endif

    lwork = nm * nm
    allocate(work(lwork))
    call dgetri(nm, ltcl, nm, ipiv, work, lwork, info )
    if (info > 0) then
        write (*,*) "WARNING: Cholesky inversion DPOTRI() exited with error code:", info
    endif

    deallocate(ipiv)
    deallocate(work)

    allocate(ltcly(nm))
    ltcly = 0.0d0
    call dgemv("N", nm, nm, 1.0d0, ltcl, nm, y, 1, 0.0d0, ltcly, 1)

    deallocate(ltcl)

    alpha = 0.0d0 
!$OMP PARALLEL DO PRIVATE(k_start,k_end,k)
    do j = 1, nm
        k_start = sum(n(:j)) - n(j) + 1
        k_end = sum(n(:j))
        do k = k_start, k_end
            alpha(k) = ltcly(j)
        enddo
    enddo
!$OMP END PARALLEL DO

    deallocate(ltcly)
end subroutine fget_alpha
