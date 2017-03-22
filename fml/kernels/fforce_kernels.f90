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

module force_tools

    contains

function outer3d(ri, rj)

    implicit none

    double precision, dimension(3), intent(in) :: ri, rj
    double precision, dimension(3,3) :: outer3d

    !write (*,*) "RI", ri(:3)
    ! write (*,*) "RJ", rj(:3)

    outer3d(:3,:3) = spread(ri(:3), dim=2, ncopies=3) * spread(rj(:3), dim=1, ncopies=3)
    ! write (*,*) "outer[0,:3] = ", outer3d(1,1), outer3d(1,2), outer3d(1,3)
    ! write (*,*) "outer[1,:3] = ", outer3d(2,1), outer3d(2,2), outer3d(2,3)
    ! write (*,*) "outer[2,:3] = ", outer3d(3,1), outer3d(3,2), outer3d(3,3)

end function outer3d

function periodic_distance(zi, zj, r_width, c_width)

    implicit none

    integer, dimension(2), intent(in) :: zi
    integer, dimension(2), intent(in) :: zj
    double precision, intent(in) :: c_width
    double precision, intent(in) :: r_width
    double precision :: periodic_distance

    double precision  :: dr, dc

    ! Row-distance
    dr = exp(-dble(abs(zi(1) - zj(1)))/r_width)

    ! Column-distance
    dc = exp(-dble(abs(zi(2) - zj(2)))/c_width)

    periodic_distance = dr * dc

end function periodic_distance


subroutine fsingle_force_kernel_inout_all(xi, xj, ni, nj, qi, qj, sigma_space, r_width, c_width, K)

    ! use force_tools, only: outer3d, periodic_distance

    implicit none

    ! Coordinates, molecule i
    double precision, dimension(:,:), intent(in) :: xi

    ! Coordinates, molecule j
    double precision, dimension(:,:), intent(in) :: xj

    ! Number of atoms, molecule i
    integer, intent(in) :: ni

    ! Number of atoms, molecule j
    integer, intent(in) :: nj

    ! Periodic table,  molecule i
    integer, dimension(:,:), intent(in) :: qi

    ! Periodic table, molecule j
    integer, dimension(:,:), intent(in) :: qj

    ! Kernel width for real-space and periodic table
    double precision, intent(in) :: sigma_space
    double precision, intent(in) :: r_width
    double precision, intent(in) :: c_width

    ! The actual kernel
    ! double precision, dimension(3, 3, nj, ni), intent(inout) :: K
    double precision, dimension(:,:,:,:), intent(inout) :: K

    ! Value of PI at full FORTRAN precision.
    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)

    ! Intermediate distance vectors
    double precision, dimension(3) :: ri, rj

    ! Norm of intermediate distance vectors
    double precision :: ri_norm, rj_norm

    ! Loop counter variables
    integer :: i, j, ii, jj

    ! Intermediate results
    double precision :: inv_L, gamma_ij, phi_ij, alpha_ij
    double precision :: quarter_inv_sigma_space2
    double precision :: half_inv_sigma_space2

    quarter_inv_sigma_space2 = 0.25d0 / (sigma_space**2)
    half_inv_sigma_space2 = 0.5d0 / (sigma_space**2)

    K = 0.0d0

    do i = 1, ni
        do j = 1, nj
            do ii = 1, ni
                do jj = 1, nj

                    if ((i == ii).or.(j == jj)) cycle

                    ri = xi(i,:) - xi(ii,:)
                    rj = xj(j,:) - xj(jj,:)

                    ri_norm = norm2(ri)
                    rj_norm = norm2(rj)

                    alpha_ij = (ri_norm**2 + rj_norm**2) * quarter_inv_sigma_space2
                    gamma_ij = (ri_norm * rj_norm) * half_inv_sigma_space2
                    phi_ij = exp(-alpha_ij) / (gamma_ij**2) * (gamma_ij * cosh(gamma_ij) - sinh(gamma_ij))
                    K(:3,:3,j,i) = K(:3,:3,j,i) + outer3d(ri, rj) * phi_ij &
                        & * periodic_distance(qi(ii,:), qj(jj,:), r_width, c_width)

                end do
            end do

            K(:3,:3,j,i) = K(:3,:3,j,i) * periodic_distance(qi(i,:), qj(j,:), r_width, c_width)
        end do
    end do

    inv_L = 1.0d0 / (2.0d0 * sqrt(pi * sigma_space**2))**3
    K = K * inv_L

end subroutine fsingle_force_kernel_inout_all


subroutine fsingle_force_kernel_inout(xi, xj, ni, nj, qi, qj, &
        & acti, actj, nacti, nactj, sigma_space, r_width, c_width, K)

    ! use force_tools, only: outer3d, periodic_distance

    implicit none

    ! Coordinates, molecule i
    double precision, dimension(:,:), intent(in) :: xi

    ! Coordinates, molecule j
    double precision, dimension(:,:), intent(in) :: xj

    ! Number of atoms, molecule i
    integer, intent(in) :: ni

    ! Number of atoms, molecule j
    integer, intent(in) :: nj

    ! Periodic table,  molecule i
    integer, dimension(:,:), intent(in) :: qi

    ! Periodic table, molecule j
    integer, dimension(:,:), intent(in) :: qj

    ! Index of active atoms, molecule i
    integer, dimension(:), intent(in) :: acti

    ! Index of active atoms, molecule j
    integer, dimension(:), intent(in) :: actj

    ! Number of active atoms, molecule i
    integer, intent(in) :: nacti

    ! Number of active atoms, molecule j
    integer, intent(in) :: nactj

    ! Kernel width for real-space and periodic table
    double precision, intent(in) :: sigma_space
    double precision, intent(in) :: r_width
    double precision, intent(in) :: c_width

    ! The actual kernel
    ! double precision, dimension(3, 3, nj, ni), intent(inout) :: K
    double precision, dimension(:,:,:,:), intent(inout) :: K

    ! Value of PI at full FORTRAN precision.
    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)

    ! Intermediate distance vectors
    double precision, dimension(3) :: ri, rj

    ! Norm of intermediate distance vectors
    double precision :: ri_norm, rj_norm

    ! Loop counter variables
    integer :: i, j, ii, jj

    ! Intermediate results
    double precision :: inv_L, gamma_ij, phi_ij, alpha_ij
    double precision :: quarter_inv_sigma_space2
    double precision :: half_inv_sigma_space2

    quarter_inv_sigma_space2 = 0.25d0 / (sigma_space**2)
    half_inv_sigma_space2 = 0.5d0 / (sigma_space**2)

    K = 0.0d0

    ! do i = 1, ni
    do i = 1, nacti
        ! do j = 1, nj
        do j = 1, nactj
            do ii = 1, ni
                do jj = 1, nj

                    if ((acti(i) == ii).or.(actj(j) == jj)) cycle

                    ri = xi(acti(i),:) - xi(ii,:)
                    rj = xj(actj(j),:) - xj(jj,:)

                    ri_norm = norm2(ri)
                    rj_norm = norm2(rj)

                    alpha_ij = (ri_norm**2 + rj_norm**2) * quarter_inv_sigma_space2
                    gamma_ij = (ri_norm * rj_norm) * half_inv_sigma_space2
                    phi_ij = exp(-alpha_ij) / (gamma_ij**2) * (gamma_ij * cosh(gamma_ij) - sinh(gamma_ij))
                    K(:3,:3,j,i) = K(:3,:3,j,i) + outer3d(ri, rj) * phi_ij &
                        & * periodic_distance(qi(ii,:), qj(jj,:), r_width, c_width)

                end do
            end do

            K(:3,:3,j,i) = K(:3,:3,j,i) * periodic_distance(qi(acti(i),:), qj(actj(j),:), &
                & r_width, c_width)
        end do
    end do

    inv_L = 1.0d0 / (2.0d0 * sqrt(pi * sigma_space**2))**3
    K = K * inv_L

end subroutine fsingle_force_kernel_inout

end module force_tools

! subroutine fsingle_force_kernel(xi, xj, ni, nj, qi, qj, sigma_space, sigma_periodic, K)
subroutine fsingle_force_kernel(xi, xj, ni, nj, qi, qj, sigma_space, r_width, c_width, K)

    use force_tools, only: outer3d, periodic_distance

    implicit none

    ! Coordinates, molecule i
    double precision, dimension(:,:), intent(in) :: xi

    ! Coordinates, molecule j
    double precision, dimension(:,:), intent(in) :: xj

    ! Number of atoms, molecule i
    integer, intent(in) :: ni

    ! Number of atoms, molecule j
    integer, intent(in) :: nj

    ! Periodic table,  molecule i
    integer, dimension(:,:), intent(in) :: qi

    ! Periodic table, molecule j
    integer, dimension(:,:), intent(in) :: qj

    ! Kernel width for real-space and periodic table
    double precision, intent(in) :: sigma_space
    double precision, intent(in) :: r_width
    double precision, intent(in) :: c_width

    ! The actual kernel
    double precision, dimension(3, 3, nj, ni), intent(out) :: K

    ! Value of PI at full FORTRAN precision.
    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)

    ! Intermediate distance vectors
    double precision, dimension(3) :: ri, rj

    ! Norm of intermediate distance vectors
    double precision :: ri_norm, rj_norm

    ! Loop counter variables
    integer :: i, j, ii, jj

    ! Intermediate results
    double precision :: inv_L, gamma_ij, phi_ij, alpha_ij
    double precision :: quarter_inv_sigma_space2
    double precision :: half_inv_sigma_space2

    quarter_inv_sigma_space2 = 0.25d0 / (sigma_space**2)
    half_inv_sigma_space2 = 0.5d0 / (sigma_space**2)

    K = 0.0d0

    do i = 1, ni
        do j = 1, nj
            do ii = 1, ni
                do jj = 1, nj

                    if ((i == ii).or.(j == jj)) cycle

                    ri = xi(:,i) - xi(:,ii)
                    rj = xj(:,j) - xj(:,jj)

                    ri_norm = norm2(ri)
                    rj_norm = norm2(rj)

                    alpha_ij = (ri_norm**2 + rj_norm**2) * quarter_inv_sigma_space2
                    gamma_ij = (ri_norm * rj_norm) * half_inv_sigma_space2
                    ! phi_ij = exp(-alpha_ij) / (gamma_ij**2) * (gamma_ij * cosh(gamma_ij) - sinh(gamma_ij))

                    phi_ij = exp(-alpha_ij) / (gamma_ij**2) * ((gamma_ij - 1.0d0) &
                        & * exp(gamma_ij) + (gamma_ij + 1.0d0) * exp(-gamma_ij))

                    K(:3,:3,j,i) = K(:3,:3,j,i) + outer3d(ri, rj) * phi_ij

                end do
            end do

            ! write (*,*), i, j, periodic_distance(qi(i,:), qj(j,:), r_width, c_width)

            K(:,:,j,i) = K(:,:,j,i) * periodic_distance(qi(i,:), qj(j,:), r_width, c_width)
        end do
    end do

    inv_L = 1.0d0 / (2.0d0 * sqrt(pi * sigma_space**2))**3
    K = K * inv_L

end subroutine fsingle_force_kernel


subroutine fcovariant_force_kernel(xi, xj, ni, nj, qi, qj, zi, zj, nmi, nmj, amax, &
    & acti, actj, nacti, nactj, actmax, sigma_space, periodic_distance_matrix, r_width, c_width, K)

    use force_tools, only: outer3d, periodic_distance

    implicit none

    ! Coordinates, molecule i
    double precision, dimension(:,:,:), intent(in) :: xi

    ! Coordinates, molecule j
    double precision, dimension(:,:,:), intent(in) :: xj

    ! Number of atoms, molecule i
    integer, dimension(:), intent(in) :: ni

    ! Number of atoms, molecule j
    integer, dimension(:), intent(in) :: nj

    ! Periodic table,  molecule i
    integer, dimension(:,:,:), intent(in) :: qi

    ! Periodic table, molecule j
    integer, dimension(:,:,:), intent(in) :: qj

    ! Atomic number,  molecule i
    integer, dimension(:,:), intent(in) :: zi

    ! Atomic number, molecule j
    integer, dimension(:,:), intent(in) :: zj

    ! Z x Z matrix of distance in the periodic table with the given r_width and c_width
    double precision, dimension(:,:), intent(in) :: periodic_distance_matrix

    ! Number of molecules i
    integer, intent(in) :: nmi

    ! Number of molecules j
    integer, intent(in) :: nmj

    ! Max atoms in one molecule
    integer, intent(in) :: amax

    ! Indexes of active atoms, molecule i
    integer, dimension(:,:), intent(in) :: acti

    ! Indexes of active atoms, molecule j
    integer, dimension(:,:), intent(in) :: actj

    ! Number of active atoms, molecule i
    integer, dimension(:), intent(in) :: nacti

    ! Number of active atoms, molecule j
    integer, dimension(:), intent(in) :: nactj

    ! Max number of active atoms in one molecule
    integer, intent(in) :: actmax

    ! Kernel width for real-space and periodic table
    double precision, intent(in) :: sigma_space
    double precision, intent(in) :: r_width
    double precision, intent(in) :: c_width

    ! The actual kernel
    ! double precision, dimension(3, 3, amax, amax, nmj, nmi), intent(out) :: K
    ! double precision, dimension(3, 3, actmax, actmax, nmj, nmi), intent(out) :: K
    double precision, dimension(nmi, nmj, actmax, actmax, 3, 3), intent(out) :: K

    ! Value of PI at full FORTRAN precision.
    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)

    ! Intermediate distance vectors
    double precision, dimension(3) :: ri, rj

    ! Norm of intermediate distance vectors
    double precision :: ri_norm, rj_norm

    ! Loop counter variables
    integer :: i, j, ii, jj

    ! Intermediate results
    double precision :: inv_L, gamma_ij, phi_ij, alpha_ij
    double precision :: quarter_inv_sigma_space2
    double precision :: half_inv_sigma_space2

    ! Intermediary kernel
    ! double precision, dimension(3, 3, amax, amax) :: Kij
    ! double precision, dimension(3, 3, actmax, actmax) :: Kij
    double precision, dimension(actmax, actmax, 3, 3) :: Kij

    integer :: iii, jjj

    double precision, dimension(nmi,actmax,amax) :: dmati
    double precision, dimension(nmj,actmax,amax) :: dmatj

    double precision :: lola, lolb

    quarter_inv_sigma_space2 = 0.25d0 / (sigma_space**2)
    half_inv_sigma_space2 = 0.5d0 / (sigma_space**2)
    inv_L = 1.0d0 / (2.0d0 * (pi * sigma_space**2))**3

    !$OMP PARALLEL DO PRIVATE(i,ii,ri,ri_norm)
    do iii = 1, nmi
        do i = 1, nacti(iii)
            do ii = 1, ni(iii)
                ri = xi(iii,acti(iii,i),:) - xi(iii,ii,:)
                ri_norm = norm2(ri)
                dmati(iii,i,ii) = ri_norm
            end do
        end do
    end do

    !$OMP PARALLEL DO PRIVATE(j,jj,rj,rj_norm)
    do jjj = 1, nmj
        do j = 1, nactj(jjj)
            do jj = 1, nj(jjj)
                rj = xj(jjj,actj(jjj,j),:) - xj(jjj,jj,:)
                rj_norm = norm2(rj)
                dmatj(jjj,j,jj) = rj_norm
            end do
        end do
    end do
    !$OMP END PARALLEL DO


    !$OMP PARALLEL DO PRIVATE(i,j,ii,jj,ri,rj,ri_norm,rj_norm,alpha_ij,gamma_ij,phi_ij,Kij)
    do jjj = 1, nmj
        do iii = 1, nmi

            Kij = 0.0d0

            do i = 1, nacti(iii)
                do j = 1, nactj(jjj)

                    do ii = 1, ni(iii)
                        do jj = 1, nj(jjj)

                            if ((acti(iii,i) == ii).or.(actj(jjj,j) == jj)) cycle

                            ri = xi(iii,acti(iii,i),:) - xi(iii,ii,:)
                            rj = xj(jjj,actj(jjj,j),:) - xj(jjj,jj,:)

                            ! ri_norm = dmati(iii,i,ii)
                            ri_norm = norm2(ri)
                            ri(:) = ri(:) / ri_norm

                            ! rj_norm = dmatj(jjj,j,jj)
                            rj_norm = norm2(rj)
                            rj(:) = rj(:) / rj_norm

                            alpha_ij = (ri_norm**2 + rj_norm**2) * quarter_inv_sigma_space2
                            gamma_ij = (ri_norm * rj_norm) * half_inv_sigma_space2

                            phi_ij = exp(-alpha_ij) / (gamma_ij**2) &
                                & * (gamma_ij *  cosh(gamma_ij) - sinh(gamma_ij))

                            ! Kij(:3,:3,j,i) = Kij(:3,:3,j,i) + outer3d(ri, rj) * phi_ij &
                            !   &  * periodic_distance_matrix(zi(iii,ii),zj(jjj,jj))
                            Kij(i,j,:3,:3) = Kij(i,j,:3,:3) + outer3d(ri, rj) * phi_ij &
                                & * periodic_distance_matrix(zi(iii,ii),zj(jjj,jj))

                        end do
                    end do

                ! Kij(:3,:3,j,i) = Kij(:3,:3,j,i) * periodic_distance_matrix(zi(iii,i),zj(jjj,j))
                    Kij(i,j,:3,:3) = Kij(i,j,:3,:3) * periodic_distance_matrix(zi(iii,i),zj(jjj,j))

                end do
            end do

            !K(:3, :3, :actmax, :actmax, jjj, iii) = K(:3, :3, :actmax, :actmax, jjj, iii) + Kij(:3, :3, :actmax, :actmax)
            K(iii, jjj, :actmax, :actmax, :3, :3) =  K(iii, jjj, :actmax, :actmax, :3, :3) + Kij(:actmax, :actmax, :3, :3)

        end do
    end do
    !$OMP END PARALLEL DO

    K = K * inv_L

end subroutine fcovariant_force_kernel
