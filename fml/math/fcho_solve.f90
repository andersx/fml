subroutine fcho_solve(A,y,x)

    implicit none

    double precision, dimension(:,:), intent(in) :: A
    double precision, dimension(:), intent(in) :: y
    double precision, dimension(:), intent(inout) :: x
    
    integer :: info, na

    na = size(A, dim=1)

    call dpotrf("U", na, A, na, info)
    call dpotrs("U", na, 1, A, na, y, na, info)

    x(:na) = y(:na)

end subroutine fcho_solve
