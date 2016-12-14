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

function distl2(X1, X2, Z1, Z2, N1, N2, width, &
    & cut_distance, r_width, c_width) result(D12)

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

    double precision :: D12 

    integer :: j_1, j_2
    
    integer :: m_1, m_2

    double precision :: r_dist, c_dist, aadist, d, maxgausdist

    double precision :: inv_cut, inv_width, sin1
    double precision :: c_width2, r_width2

    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)

    inv_cut = pi / (2.0d0 * cut_distance)
    inv_width = 1.0d0 / (4.0d0 * width**2)
    maxgausdist = 8.0d0 * width
    r_width2 = r_width**2
    c_width2 = c_width**2

    D12 = 0.0d0

    do j_1 = 1, N1
        do j_2 = 1, N2

            aadist = 0.0d0

            do m_1 = 1, N1

                ! if (X1(1, m_1, j_1) > cut_distance) exit
                if (X1(j_1, 1,m_1) > cut_distance) exit

                ! sin1 = 1.0d0 - sin(x1(1, m_1, j_1) * inv_cut)
                sin1 = 1.0d0 - sin(x1(j_1,1,m_1) * inv_cut)

                do m_2 = 1, N2

                    ! if (X2(1, m_2, j_2) > cut_distance) exit
                    if (X2(j_2,1, m_2) > cut_distance) exit

                    ! if (abs(X2(1, m_2, j_2) - X1(1, m_1, j_1)) < maxgausdist) then
                    if (abs(X2(j_2,1,m_2) - X1(j_1,1,m_1)) < maxgausdist) then

                        ! d = exp(-((x1(1, m_1 ,j_1)-x2(1, m_2, j_2))**2) * inv_width ) *  & 
                        !     & sin1 * (1.0d0 - sin(x2(1, m_2, j_2) * inv_cut))
                        d = exp(-((x1(j_1,1,m_1)-x2(j_2,1,m_2))**2) * inv_width ) *  & 
                            & sin1 * (1.0d0 - sin(x2(j_2,1,m_2) * inv_cut))

                        ! r_dist = abs(x1(2, m_1, j_1) - x2(2,m_2, j_2))
                        ! c_dist = abs(x1(3, m_1, j_1) - x2(3,m_2, j_2))
                        r_dist = abs(x1(j_1,2,m_1) - x2(j_2,2,m_2))
                        c_dist = abs(x1(j_1,3,m_1) - x2(j_2,3,m_2))
    
                        ! d = d * (r_width2/(r_width2 + r_dist**2) * c_width2/(c_width2 + c_dist**2))
                        d = d * (r_width2/(r_width2 + r_dist**2) * c_width2/(c_width2 + c_dist**2))


                        ! aadist = aadist + d * (1.0d0 + x1(4, m_1, j_1)*x2(4, m_2, j_2) + & 
                        !     & x1(5, m_1, j_1)*x2(5, m_2, j_2))
                        aadist = aadist + d * (1.0d0 + x1(j_1,4,m_1)*x2(j_2,4,m_2) + & 
                            & x1(j_1,5,m_1)*x2(j_2,5,m_2))

                    end if
                end do
            end do

            r_dist = abs(z1(j_1,1) - z2(j_2,1))
            c_dist = abs(z1(j_1,2) - z2(j_2,2))

            D12 = D12 + aadist * (r_width2/(r_width2 + r_dist**2) * c_width2/(c_width2 + c_dist**2))

        enddo
    enddo

end function distl2

end module funcs

subroutine molecular_arad_l2_distance_all(X1, X2, Z1, Z2, N1, N2, nmol1, nmol2, width, &
    & cut_distance, r_width, c_width, D12)

    use funcs, only: distl2, m_dist, stoch_dist

    implicit none

    double precision, dimension(:,:,:,:), intent(in) :: X1
    double precision, dimension(:,:,:,:), intent(in) :: X2

    integer, dimension(:,:,:), intent(in) :: Z1
    integer, dimension(:,:,:), intent(in) :: Z2

    integer, dimension(:), intent(in) :: N1
    integer, dimension(:), intent(in) :: N2

    integer, intent(in) :: nmol1
    integer, intent(in) :: nmol2

    double precision, intent(in) :: width
    double precision, intent(in) :: cut_distance
    double precision, intent(in) :: r_width
    double precision, intent(in) :: c_width

    double precision, dimension(nmol1, nmol2), intent(out) :: D12 

    double precision, dimension(nmol1) :: D11
    double precision, dimension(nmol2) :: D22

    integer :: i, j

    do i = 1, nmol1
        ! D11(i) = distl2(X1(:,:,:,i), X1(:,:,:,i), Z1(i,:,:), Z1(i,:,:), N1(i), N1(i), &
        D11(i) = distl2(X1(i,:,:,:), X1(i,:,:,:), Z1(i,:,:), Z1(i,:,:), N1(i), N1(i), &
            & width, cut_distance, r_width, c_width)
    enddo 
    do i = 1, nmol2
        ! D22(i) = distl2(X2(:,:,:,i), X2(:,:,:,i), Z2(i,:,:), Z2(i,:,:), N2(i), N2(i), &
        D22(i) = distl2(X2(i,:,:,:), X2(i,:,:,:), Z2(i,:,:), Z2(i,:,:), N2(i), N2(i), &
            & width, cut_distance, r_width, c_width)
    enddo

    do j = 1, nmol2
        do i = 1, nmol1
            ! D12(i,j) = distl2(X1(:,:,:,i), X2(:,:,:,j), Z1(i,:,:), Z2(j,:,:), N1(i), N2(j), &
            D12(i,j) = distl2(X1(i,:,:,:), X2(j,:,:,:), Z1(i,:,:), Z2(j,:,:), N1(i), N2(j), &
                & width, cut_distance, r_width, c_width)

            D12(i,j) = D11(i) + D22(j) - 2.0d0 * D12(i,j)
        enddo
    enddo

end subroutine molecular_arad_l2_distance_all

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
    double precision :: D11, D22, D12

    double precision :: dd

    double precision :: rdist
    double precision :: cdist

    D11 = 0.0d0
    write (*,*) "N1:", N1
    write (*,*) "N2:", N2

    do j_1 = 1, N1
        write(*,*) "Z1:", "i", j_1, Z1(j_1,1), Z1(j_1,2)
    enddo
        

    do j_1 = 1, N1
        do j_2 = 1, N1

            rdist = abs(z1(j_1,1) - z1(j_2,1))
            CDist = abs(Z1(j_1,2) - Z1(j_2,2))

            dd = M_Dist(X1(j_1,:,:),X1(j_2,:,:),n1,n1,width,cut_distance,R_Width,C_width)
            D11 = D11 + dd * stoch_dist(RDist,CDist,R_Width,C_width)

        enddo
    enddo

    write(*,*) "D11:", "i", D11
    D22 = 0.0d0

    do j_1 = 1, N2
        do j_2 = 1, N2

            rdist = abs(z2(j_1,1) - z2(j_2,1))
            CDist = abs(Z2(j_1,2) - Z2(j_2,2))

            dd = M_Dist(X2(j_1,:,:),X2(j_2,:,:),n2,n2,width,cut_distance,R_Width,C_width)
            D22 = D22 + dd * stoch_dist(RDist,CDist,R_Width,C_width)

        enddo
    enddo

    D12 = 0.0d0

    do j_1 = 1, N1
        do j_2 = 1, N2

            rdist = abs(z1(j_1,1) - z2(j_2,1))
            CDist = abs(Z1(j_1,2) - Z2(j_2,2))

            dd = M_Dist(X1(j_1,:,:),X2(j_2,:,:),n1,n2,width,cut_distance,R_Width,C_width)
            D12 = D12 + dd * stoch_dist(RDist,CDist,R_Width,C_width)

        enddo
    enddo

    distance = D11 + D22 - 2.0d0 * D12

end subroutine molecular_arad_l2_distance

subroutine atomic_arad_l2_distance(X1, X2, Z1, Z2, N1, N2, width, &
    & cut_distance, r_width, c_width, D12)

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

    double precision, dimension(N1,N2), intent(out) :: D12

    integer :: j_1, j_2

    double precision, dimension(N1) :: D1
    double precision, dimension(N2) :: D2

    double precision :: dd

    double precision :: rdist
    double precision :: cdist


    do j_1 = 1, N1
            D1(j_1) = M_Dist(X1(j_1,:,:),X1(j_1,:,:),n1,n1,width,cut_distance,R_Width,C_width)
    enddo

    do j_1 = 1, N2
            D2(j_1) = M_Dist(X2(j_1,:,:),X2(j_1,:,:),n2,n2,width,cut_distance,R_Width,C_width)
    enddo

    do j_1 = 1, N1
        do j_2 = 1, N2

            rdist = abs(z1(j_1,1) - z2(j_2,1))
            CDist = abs(Z1(j_1,2) - Z2(j_2,2))

            dd = M_Dist(X1(j_1,:,:),X2(j_2,:,:),n1,n2,width,cut_distance,R_Width,C_width)
            D12(j_1, j_2) = D1(j_1) + D2(j_2) - 2.0d0 * dd * stoch_dist(RDist,CDist,R_Width,C_width)

        enddo
    enddo


end subroutine atomic_arad_l2_distance
