subroutine fgaussian_kernel(a, na, b, nb, k, sigma)

    double precision, dimension(:,:), intent(in) :: a
    double precision, dimension(:,:), intent(in) :: b

    integer, intent(in) :: na, nb

    double precision, dimension(:,:), intent(inout) :: k
    double precision, intent(in) :: sigma

    double precision, allocatable, dimension(:) :: temp

    double precision :: inv_sigma
    integer :: l

    inv_sigma = -0.5d0 / (sigma*sigma)

    allocate(temp(size(a, dim=1)))

!$OMP PARALLEL DO PRIVATE(temp)
    do i = 1, nb 
        do j = 1, na
            temp(:) = a(:,j) - b(:,i)
            K(j,i) = exp(inv_sigma * sqrt(sum(temp*temp)))
        enddo
    enddo
!$OMP END PARALLEL DO

    deallocate(temp)

end subroutine fgaussian_kernel

subroutine flaplacian_kernel(a, na, b, nb, k, sigma)

    double precision, dimension(:,:), intent(in) :: a
    double precision, dimension(:,:), intent(in) :: b

    integer, intent(in) :: na, nb

    double precision, dimension(:,:), intent(inout) :: k
    double precision, intent(in) :: sigma

    double precision :: inv_sigma

    inv_sigma = -1.0d0 / sigma

!$OMP PARALLEL DO
    do i = 1, nb 
        do j = 1, na
        k(j,i) = exp(inv_sigma * sum(abs(a(:,j) - b(:,i))))
        enddo
    enddo
!$OMP END PARALLEL DO

end subroutine flaplacian_kernel
