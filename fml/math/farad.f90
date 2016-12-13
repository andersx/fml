module funcs

contains

function stoch_dist(D1, D2, W1, W2) result(a)

    implicit none

    double precision, intent(in) :: D1
    double precision, intent(in) :: D2
    double precision, intent(in) :: W1
    double precision, intent(in) :: W2

    double precision :: a

    a = W1**2/(W1**2 + D1**2) * W2**2/(W2**2 + D2**2)

end function stoch_dist


function dist(x1, x2, w1, w2) result(a)

    implicit none

    double precision, intent(in) :: x1
    double precision, intent(in) :: x2
    double precision, intent(in) :: w1
    double precision, intent(in) :: w2

    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)

    double precision :: a

    a = exp(-((x1-x2)**2)/(4.0d0*w1**2)) *  (1.0d0 - sin(pi * x1/(2.0d0 * w2))) * &
        & (1.0d0 - sin(pi * x2/(2.0d0 * w2)))

end function dist


function m_dist(X1, X2, N1, N2, width, cut_distance, r_width, c_width) result(aadist)

    implicit none

    double precision, dimension(:,:), intent(in) :: X1
    double precision, dimension(:,:), intent(in) :: X2

    integer, intent(in) :: N1
    integer, intent(in) :: N2

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

    maxgausdist = 8.0d0*width
    aadist = 0.0d0

    do m_1 = 1, N1

        if (X1(1, m_1) > cut_distance) exit

        do m_2 = 1, N2

            if (X2(1, m_2) > cut_distance) exit

            if (abs(X2(1,m_2) - X1(1,m_1)) < maxgausdist) then

                r_dist = abs(x1(2,m_1) - x2(2,m_2))
                c_dist = abs(x1(3,m_1) - x2(3,m_2))

                d = dist(x1(1,m_1), x2(1,m_2),width,cut_distance)
                d = d * stoch_dist(r_dist,c_dist,r_width,c_width)
                aadist = aadist + d * (1.0d0 + x1(4,m_1)*x2(4,m_2) + x1(5,m_1)*x2(5,m_2))

            end if
        end do
    end do

end function m_dist

end module funcs

subroutine molecular_arad_l2_distance(X1, X2, Z1, Z2, N1, N2, width, &
    & cut_distance, r_width, c_width, distance)

    use funcs, only: m_dist, stoch_dist

    implicit none

    double precision, dimension(:,:,:), intent(in) :: X1
    double precision, dimension(:,:,:), intent(in) :: X2

    integer, dimension(:,:), intent(in) :: Z1
    integer, dimension(:,:), intent(in) :: Z2

    integer, intent(in) :: N1
    integer, intent(in) :: N2

    double precision, intent(in) :: width
    double precision, intent(in) :: cut_distance
    double precision, intent(in) :: r_width
    double precision, intent(in) :: c_width

    double precision, intent(out) :: distance

    integer :: j_1, j_2
    double precision :: D1, D2

    double precision :: dd
    double precision :: pair

    double precision :: rdist
    double precision :: cdist

    D1 = 0.0d0

    do j_1 = 1, N1
        do j_2 = 1, N1

            rdist = abs(z1(j_1,1) - z1(j_2,1))
            CDist = abs(Z1(j_1,2) - Z1(j_2,2))

            dd = M_Dist(X1(j_1,:,:),X1(j_2,:,:),n1,n1,width,cut_distance,R_Width,C_width)
            D1 = D1 + dd * stoch_dist(RDist,CDist,R_Width,C_width)

        enddo
    enddo

    D2 = 0.0d0

    do j_1 = 1, N2
        do j_2 = 1, N2

            rdist = abs(z2(j_1,1) - z2(j_2,1))
            CDist = abs(Z2(j_1,2) - Z2(j_2,2))

            dd = M_Dist(X2(j_1,:,:),X2(j_2,:,:),n2,n2,width,cut_distance,R_Width,C_width)
            D2 = D2 + dd * stoch_dist(RDist,CDist,R_Width,C_width)

        enddo
    enddo

    pair = 0.0d0

    do j_1 = 1, N1
        do j_2 = 1, N2

            rdist = abs(z1(j_1,1) - z2(j_2,1))
            CDist = abs(Z1(j_1,2) - Z2(j_2,2))

            dd = M_Dist(X1(j_1,:,:),X2(j_2,:,:),n1,n2,width,cut_distance,R_Width,C_width)
            dd = dd * stoch_dist(RDist,CDist,R_Width,C_width)
            pair = pair + dd

        enddo
    enddo

    distance = D1 + D2 - 2.0d0 * pair

end subroutine molecular_arad_l2_distance

subroutine atomic_arad_l2_distance(X1, X2, Z1, Z2, N1, N2, width, &
    & cut_distance, r_width, c_width, distance)

    use funcs, only: m_dist, stoch_dist

    implicit none

    double precision, dimension(:,:,:), intent(in) :: X1
    double precision, dimension(:,:,:), intent(in) :: X2

    integer, dimension(:,:), intent(in) :: Z1
    integer, dimension(:,:), intent(in) :: Z2

    integer, intent(in) :: N1
    integer, intent(in) :: N2

    double precision, intent(in) :: width
    double precision, intent(in) :: cut_distance
    double precision, intent(in) :: r_width
    double precision, intent(in) :: c_width

    double precision, intent(out) :: distance

    integer :: j_1, j_2
    double precision :: D1, D2

    double precision :: dd
    double precision :: pair

    double precision :: rdist
    double precision :: cdist

    D1 = 0.0d0

    do j_1 = 1, N1
        do j_2 = 1, N1

            rdist = abs(z1(j_1,1) - z1(j_2,1))
            CDist = abs(Z1(j_1,2) - Z1(j_2,2))

            dd = M_Dist(X1(j_1,:,:),X1(j_2,:,:),n1,n1,width,cut_distance,R_Width,C_width)
            D1 = D1 + dd * stoch_dist(RDist,CDist,R_Width,C_width)

        enddo
    enddo

    D2 = 0.0d0

    do j_1 = 1, N2
        do j_2 = 1, N2

            rdist = abs(z2(j_1,1) - z2(j_2,1))
            CDist = abs(Z2(j_1,2) - Z2(j_2,2))

            dd = M_Dist(X2(j_1,:,:),X2(j_2,:,:),n2,n2,width,cut_distance,R_Width,C_width)
            D2 = D2 + dd * stoch_dist(RDist,CDist,R_Width,C_width)

        enddo
    enddo

    pair = 0.0d0

    do j_1 = 1, N1
        do j_2 = 1, N2

            rdist = abs(z1(j_1,1) - z2(j_2,1))
            CDist = abs(Z1(j_1,2) - Z2(j_2,2))

            dd = M_Dist(X1(j_1,:,:),X2(j_2,:,:),n1,n2,width,cut_distance,R_Width,C_width)
            dd = dd * stoch_dist(RDist,CDist,R_Width,C_width)
            pair = pair + dd

        enddo
    enddo

    distance = D1 + D2 - 2.0d0 * pair

end subroutine atomic_arad_l2_distance
