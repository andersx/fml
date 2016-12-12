
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

! def _dist(x1,x2,w1,w2):
!     return math.exp(-((x1-x2)**2)/(4*w1**2))*(1 - math.sin(np.pi * x1/(2 * w2)))*(1 - math.sin(np.pi * x2/(2 * w2)))
! 
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

    ! double precision :: stoch_dist
    ! double precision :: dist

    double precision :: aadist

    double precision :: maxgausdist
    double precision :: d

    double precision :: r_dist
    double precision :: c_dist

    integer :: m_1, m_2

    maxgausdist = 8.0d0*width
    aadist = 0.0d0

    do m_1 = 1, N1
    ! for m_1 in range(len(X1[0])):
   
        if (X1(1, m_1) > cut_distance) exit
        ! if X1[0,m_1] > cutDist:
        !    break

        do m_2 = 1, N2
        ! for m_2 in range(len(X2[0])):

            if (X1(1, m_1) > cut_distance) exit

            ! if X2[0,m_2] > cutDist:
            !     break

            if (abs(X2(1,m_2) - X1(1,m_1)) < maxgausdist) then
            ! if  abs(X2[0,m_2] - X1[0,m_1]) < maxGausDist:

                r_dist = abs(x1(2,m_1) - x2(2,m_2))
                c_dist = abs(x1(3,m_1) - x2(3,m_2))
            !     RDist = abs(X1[1,m_1] - X2[1,m_2])
            !     CDist = abs(X1[2,m_1] - X2[2,m_2])

                d = dist(x1(1,m_1), x2(1,m_2),width,cut_distance)
                d = d * stoch_dist(r_dist,c_dist,r_width,c_width)
                aadist = aadist + d * (1.0d0 + x1(4,m_1)*x2(4,m_2) + x1(5,m_1)*x2(5,m_2))
            !     d = _dist(X1[0,m_1], X2[0,m_2],width,cutDist)
            !     d = d *_StochDist(RDist,CDist,RWidth,Cwidth)
            !     AAdist += d * (1 + X1[3,m_1]*X2[3,m_2] + X1[4,m_1]*X2[4,m_2])

            end if
        end do
    end do


end function m_dist

end module funcs

subroutine molecular_arad_l2_distance(X1, X2, Z1, Z2, N1, N2, width, &
    & cut_distance, r_width, c_width, distance)

    use funcs

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

    ! double precision :: stoch_dist
    ! double precision :: dist
    ! double precision :: m_dist

    integer :: j_1, j_2
    double precision :: D1, D2
    
    double precision :: dd

    double precision :: rdist
    double precision :: cdist

    distance = 0.0d0
    ! DistMatrix  = zeros((len(X1), len(X1)))

    do j_1 = 1, N1
    ! for j_1 in  range(len(X1)):

        ! if Z1[j_1,0] == -1:
        !    break

        rdist = abs(z1(j_1,1) - z1(j_2,1))
        ! RDist = abs(Z1[j_1,0] - Z1[j_2,0])

        CDist = abs(Z1(j_1,2) - Z1(j_2,2))
        ! CDist = abs(Z1[j_1,1] - Z1[j_2,1])

        dd = M_Dist(X1(j_1,:,:),X1(j_2,:,:),n1,n2,width,cut_distance,R_Width,C_width)

        print *, "fdd", dd
        dd = dd * stoch_dist(RDist,CDist,R_Width,C_width)

        print *, "fdd", dd
        distance = distance + dd
        ! DistMatrix[j_1,j_2] = dd

    enddo

end subroutine molecular_arad_l2_distance
