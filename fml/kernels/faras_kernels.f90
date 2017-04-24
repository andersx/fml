module aras_utils

    implicit none

contains

function atomic_distl2(X1, X2, N1, N2, ksi1, ksi2, sin1, sin2, cos1, cos2, &
    & t_width, r_width, c_width, d_width, cut_distance, order, pd, ang_norm2) result(aadist)

    implicit none

    double precision, dimension(:,:), intent(in) :: X1
    double precision, dimension(:,:), intent(in) :: X2

    integer, intent(in) :: N1
    integer, intent(in) :: N2

    double precision, dimension(:), intent(in) :: ksi1
    double precision, dimension(:), intent(in) :: ksi2

    double precision, dimension(:,:,:), intent(in) :: sin1
    double precision, dimension(:,:,:), intent(in) :: sin2
    double precision, dimension(:,:,:), intent(in) :: cos1
    double precision, dimension(:,:,:), intent(in) :: cos2

    double precision, intent(in) :: t_width
    double precision, intent(in) :: r_width
    double precision, intent(in) :: c_width 
    double precision, intent(in) :: d_width 
    double precision, intent(in) :: cut_distance
    integer, intent(in) :: order
    double precision, dimension(:,:), intent(in) :: pd

    double precision :: aadist

    double precision :: d

    integer :: m_1, m_2

    integer :: i, m, p1, p2, n

    double precision :: angular 
    ! double precision :: angular2

    ! double precision :: dcos1, dcos2, dsin1, dsin2
    double precision :: maxgausdist2

    integer :: pmax1
    integer :: pmax2

    double precision :: inv_width !, diff, dist
    double precision :: c_width2, r_width2, r2

    double precision, dimension(order) :: s
    
    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)

    ! double precision, parameter :: g1 = 2.0d0
    ! double precision, parameter :: a0 = g1 / 2.0d0
    double precision :: g1 
    double precision :: a0 

    logical, allocatable, dimension(:) :: mask1
    logical, allocatable, dimension(:) :: mask2

    double precision :: temp, sin1_temp, cos1_temp, ang_norm1

    double precision, intent(in):: ang_norm2
    
    double precision :: cos_sum1, cos_sum2


    ! write(*,*) "ANG_NORM", t_width, ang_norm2

    pmax1 = int(maxval(x1(2,:)))
    pmax2 = int(maxval(x2(2,:)))

    allocate(mask1(pmax1))
    allocate(mask2(pmax2))

    mask1 = .true.
    mask2 = .true.

    do i = 1, n1
        mask1(int(x1(2,i))) = .false.
    enddo
        
    do i = 1, n2
        mask2(int(x2(2,i))) = .false.
    enddo

    a0 = 0.0d0
    g1 = sqrt(2.0d0 * pi)/ang_norm2
    ! g1 = 1.0d0/ang_norm2
    ! g1 = pi/ang_norm2

    do m = 1, order
        s(m) = g1 * exp(-(t_width * m)**2 / 2.0d0)
        ! write (*,*) "S(M)", m, s(m)
    enddo 

    inv_width = -1.0d0 / (4.0d0 * d_width**2)

    maxgausdist2 = (8.0d0 * d_width)**2
    r_width2 = r_width**2
    c_width2 = c_width**2

    aadist = 0.0d0

    ! ang_norm1 = 1.0d0! / (sum(ksi1(:N1)) * sum(ksi2(:N2)))


    do m_1 = 1, N1

        if (X1(1, m_1) > cut_distance) exit

        do m_2 = 1, N2

            if (X2(1, m_2) > cut_distance) exit

            r2 = (X2(1,m_2) - X1(1,m_1))**2

            if (r2 < maxgausdist2) then

                d = exp(r2 * inv_width )  * ksi1(m_1) * ksi2(m_2) * &
                    & pd(int(x1(2,m_1)), int(x2(2,m_2)))


                angular = a0 * a0

                ! cos_sum1 = 0.0d0
                ! cos_sum2 = 0.0d0
                
                do m = 1, order

                    temp = 0.0d0

                    do p1 = 1, pmax1
                        if (mask1(p1)) cycle
                        cos1_temp = cos1(m,m_1,p1)
                        sin1_temp = sin1(m,m_1,p1)

                        do p2 = 1, pmax2
                            if (mask2(p2)) cycle

                        ! angular = angular + s(m) * (cos1(m,m_1,p1) * cos2(m,m_2,p2) + sin1(m,m_1,p1) * sin2(m,m_2,p2))!* pd(p2,p1)
                        temp = temp + (cos1_temp * cos2(m,m_2,p2) + sin1_temp * sin2(m,m_2,p2)) * pd(p2,p1)

                        enddo 
                    enddo 

                    angular = angular + temp * s(m)
                    ! angular = angular + s(m) * (sum(cos1(m,m_1,:)) * sum(cos2(m,m_2,:)) + sum(sin1(m,m_1,:)) * sum(sin2(m,m_2,:))) 

                    ! write (*,*) m_1, m_2, sum(cos1(m,m_1,:)), sum(cos2(m,m_2,:))
                enddo 
                

                ! write (*,*) cos_sum1, cos_sum2
                ! Note: this 0.5 is an empirical hyperparameter
                aadist = aadist + d * (1.0d0 + angular)

            end if
        end do
    end do
    

    deallocate(mask1)
    deallocate(mask2)
end function atomic_distl2


end module aras_utils


subroutine fget_kernels_aras(x1, x2, n1, n2, sigmas, nm1, nm2, nsigmas, &
       & t_width, r_width, c_width, d_width, cut_distance, order, pd, kernels)

    use aras_utils, only: atomic_distl2

    implicit none

    ! ARAD descriptors for the training set, format (i,j_1,5,m_1)
    double precision, dimension(:,:,:,:), intent(in) :: x1
    double precision, dimension(:,:,:,:), intent(in) :: x2

    ! List of numbers of atoms in each molecule
    integer, dimension(:), intent(in) :: n1
    integer, dimension(:), intent(in) :: n2

    ! Sigma in the Gaussian kernel
    double precision, dimension(:), intent(in) :: sigmas

    ! Number of molecules
    integer, intent(in) :: nm1
    integer, intent(in) :: nm2

    ! Number of sigmas
    integer, intent(in) :: nsigmas

    double precision, intent(in) :: t_width
    double precision, intent(in) :: r_width
    double precision, intent(in) :: c_width 
    double precision, intent(in) :: d_width 
    double precision, intent(in) :: cut_distance
    integer, intent(in) :: order

    ! -1.0 / sigma^2 for use in the kernel
    double precision, dimension(nsigmas) :: inv_sigma2

    double precision, dimension(:,:), intent(in) :: pd

    ! Resulting alpha vector
    double precision, dimension(nsigmas,nm1,nm2), intent(out) :: kernels

    ! Internal counters
    integer :: i, j, k, ni, nj
    integer :: m_1, i_1, j_1, a, m, n


    ! Pre-computed constants
    double precision :: r_width2
    double precision :: c_width2
    double precision :: inv_cut

    ! Temporary variables necessary for parallelization
    double precision :: l2dist
    double precision, allocatable, dimension(:,:) :: atomic_distance

    ! Pre-computed terms in the full distance matrix
    double precision, allocatable, dimension(:,:) :: selfl21
    double precision, allocatable, dimension(:,:) :: selfl22

    ! Pre-computed sine terms
    double precision, allocatable, dimension(:,:,:) :: ksi1
    double precision, allocatable, dimension(:,:,:) :: ksi2

    double precision, allocatable, dimension(:,:,:,:,:) :: sinp1
    double precision, allocatable, dimension(:,:,:,:,:) :: sinp2
    double precision, allocatable, dimension(:,:,:,:,:) :: cosp1
    double precision, allocatable, dimension(:,:,:,:,:) :: cosp2

    ! Value of PI at full FORTRAN precision.
    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)

    ! Small number // to test for numerical stability
    double precision, parameter :: eps = 5.0d-12

    ! counter for periodic distance
    integer :: p
    integer :: pmax1
    integer :: pmax2
    double precision :: theta

    double precision :: ang_norm2
    ang_norm2 = 0.0d0

    do n = -10000, 10000
        ! ang_norm2 = ang_norm2 + t_width**2 * 2.0d0 * pi &
        ! & * exp(-((t_width * n)**2)) * (2.0d0 - 2.0d0 * cos(n * pi))
        ang_norm2 = ang_norm2 + exp(-((t_width * n)**2)) * (2.0d0 - 2.0d0 * cos(n * pi))
    end do
    ! ang_norm2 = sqrt(ang_norm2)
    ang_norm2 = sqrt(ang_norm2 * pi) * 2.0d0
    ! write(*,*) "ANGNORM", ang_norm2

    pmax1 = 0
    pmax2 = 0

    do a = 1, nm1
        pmax1 = max(pmax1, int(maxval(x1(a,1,2,:n1(a)))))
    enddo
    do a = 1, nm2
        pmax2 = max(pmax2, int(maxval(x2(a,1,2,:n1(a)))))
    enddo


    r_width2 = r_width**2
    c_width2 = c_width**2

    inv_cut = pi / (2.0d0 * cut_distance)
    inv_sigma2(:) = -1.0d0 / (sigmas(:))**2

    allocate(ksi1(nm1, maxval(n1), maxval(n1)))
    allocate(ksi2(nm2, maxval(n2), maxval(n2)))

    ksi1 = 0.0d0
    ksi2 = 0.0d0

    ! write (*,*) "INV_CUT", inv_cut

    !$OMP PARALLEL DO PRIVATE(ni)
    do i = 1, nm1
        ni = n1(i)
        do m_1 = 1, ni
            do i_1 = 1, ni
                if (x1(i,i_1,1,m_1) < cut_distance) then
                    ksi1(i, i_1, m_1) = 1.0d0 - sin(x1(i,i_1,1,m_1) * inv_cut)
                endif
            enddo
        enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(ni)
    do i = 1, nm2
        ni = n2(i)
        do m_1 = 1, ni
            do i_1 = 1, ni
                if (x2(i,i_1,1,m_1) < cut_distance) then
                    ksi2(i, i_1, m_1) = 1.0d0 - sin(x2(i,i_1,1,m_1) * inv_cut)
                endif
            enddo
        enddo
    enddo
    !$OMP END PARALLEL DO

    allocate(cosp1(nm1, maxval(n1), order, maxval(n1), pmax1))
    allocate(sinp1(nm1, maxval(n1), order, maxval(n1), pmax1))

    cosp1 = 0.0d0
    sinp1 = 0.0d0

    do a = 1, nm1
        ni = n1(a)
        do m = 1, order
            do i = 1, ni
                do j = 1, ni
                    do k = 1, ni
                        if ((j /= k) .and. (j /= 1) .and. (k /= 1)) then
                        p =  int(x1(a,i,2,k))
                        ! cosp1(a, i, m, j, p) = cosp1(a, i, m, j, p) + cos(m * x1(a,i,j+3,k))
                        ! sinp1(a, i, m, j, p) = sinp1(a, i, m, j, p) + sin(m * x1(a,i,j+3,k))
                        theta = x1(a,i,j+3,k)
                        cosp1(a, i, m, j, p) = cosp1(a, i, m, j, p) + &
                            & (cos(m * theta) - cos((theta + pi) * m))*ksi1(a,i,k)
    
            !write (*,*) "COS1", j,k, (cos(m * theta) - cos((theta + pi) * m))*ksi1(a,i,k)

                        sinp1(a, i, m, j, p) = sinp1(a, i, m, j, p) + &
                            & (sin(m * theta) - sin((theta + pi) * m))*ksi1(a,i,k)

                        endif
                    enddo
                    ! write (*,*) "COS1", cosp1(a, i, m, j, p)
                enddo
            enddo
        enddo
    enddo

    allocate(cosp2(nm2, maxval(n2), order, maxval(n2), pmax2))
    allocate(sinp2(nm2, maxval(n2), order, maxval(n2), pmax2))

    cosp2 = 0.0d0
    sinp2 = 0.0d0

    do a = 1, nm2
        ni = n2(a)
        do m = 1, order
            do i = 1, ni
                do j = 1, ni
                    do k = 1, ni
                        if ((j /= k) .and. (j /= 1) .and. (k /= 1)) then
                        p =  int(x2(a,i,2,k))
                        ! cosp2(a, i, m, j, p) = cosp2(a, i, m, j, p) + cos(m * x2(a,i,j+3,k))
                        ! sinp2(a, i, m, j, p) = sinp2(a, i, m, j, p) + sin(m * x2(a,i,j+3,k))
                        theta = x2(a,i,j+3,k)
                        cosp2(a, i, m, j, p) = cosp2(a, i, m, j, p) + &
                            & (cos(m * theta) - cos((theta + pi) * m)) * ksi2(a,i,k)
                        sinp2(a, i, m, j, p) = sinp2(a, i, m, j, p) + &
                            & (sin(m * theta) - sin((theta + pi) * m)) * ksi2(a,i,k)
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo

    allocate(selfl21(nm1, maxval(n1)))
    allocate(selfl22(nm2, maxval(n2)))

    !$OMP PARALLEL DO PRIVATE(ni)
    do i = 1, nm1
        ni = n1(i)
        do i_1 = 1, ni
            selfl21(i,i_1) = atomic_distl2(x1(i,i_1,:,:), x1(i,i_1,:,:), n1(i), n1(i), &
                & ksi1(i,i_1,:), ksi1(i,i_1,:), sinp1(i,i_1,:,:,:), sinp1(i,i_1,:,:,:), cosp1(i,i_1,:,:,:), cosp1(i,i_1,:,:,:), &
            & t_width, r_width, c_width, d_width, cut_distance, order, pd, ang_norm2)
            !    & t_width, r_width, c_width, d_width, cut_distance, order, pd)
        enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(ni)
    do i = 1, nm2
        ni = n2(i)
        do i_1 = 1, ni
            selfl22(i,i_1) = atomic_distl2(x2(i,i_1,:,:), x2(i,i_1,:,:), n2(i), n2(i), &
                & ksi2(i,i_1,:), ksi2(i,i_1,:), sinp2(i,i_1,:,:,:), sinp2(i,i_1,:,:,:), cosp2(i,i_1,:,:,:), cosp2(i,i_1,:,:,:), &
            & t_width, r_width, c_width, d_width, cut_distance, order, pd, ang_norm2)
        enddo
    enddo
    !$OMP END PARALLEL DO


    allocate(atomic_distance(maxval(n1), maxval(n2)))

    kernels(:,:,:) = 0.0d0
    atomic_distance(:,:) = 0.0d0

    !$OMP PARALLEL DO PRIVATE(l2dist,atomic_distance,ni,nj)
    do j = 1, nm2
        nj = n2(j)
        do i = 1, nm1
            ni = n1(i)

            atomic_distance(:,:) = 0.0d0

            do i_1 = 1, ni
                do j_1 = 1, nj

                    l2dist = atomic_distl2(x1(i,i_1,:,:), x2(j,j_1,:,:), n1(i), n2(j), &
                        & ksi1(i,i_1,:), ksi2(j,j_1,:), sinp1(i,i_1,:,:,:), sinp2(j,j_1,:,:,:), &
                        & cosp1(i,i_1,:,:,:), cosp2(j,j_1,:,:,:), &
            & t_width, r_width, c_width, d_width, cut_distance, order, pd, ang_norm2)

                    l2dist = selfl21(i,i_1) + selfl22(j,j_1) - 2.0d0 * l2dist * pd(int(x1(i,i_1,2,1)),int(x2(j,j_1,2,1)))

                    if (abs(l2dist) < eps) l2dist = 0.0d0

                    atomic_distance(i_1,j_1) = l2dist
                    
                        ! write (*,*) i_1, j_1, l2dist
                enddo
            enddo

            do k = 1, nsigmas
                kernels(k, i, j) =  sum(exp(atomic_distance(:ni,:nj) * inv_sigma2(k)))
            enddo

        enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(atomic_distance)
    deallocate(selfl21)
    deallocate(selfl22)
    deallocate(ksi1)
    deallocate(ksi2)
    deallocate(cosp1)
    deallocate(cosp2)
    deallocate(sinp1)
    deallocate(sinp2)

end subroutine fget_kernels_aras