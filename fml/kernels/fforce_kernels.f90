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

    outer3d(:3,:3) = spread(ri(:3), dim=1, ncopies=3) * spread(rj(:3), dim=2, ncopies=3)

end function outer3d

function periodic_distance(zi, zj, c_width, r_width)

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


subroutine fforce_kernel(xi, xj, ni, nj, qi, qj, nmi, nmj, amax, &
        & acti, actj, nacti, nactj, actmax, sigma_space, r_width, c_width, K)

    use force_tools, only: fsingle_force_kernel_inout

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
    double precision, dimension(3, 3, actmax, actmax, nmj, nmi), intent(out) :: K

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
    double precision, dimension(3, 3, actmax, actmax) :: Kij

    do j = 1, nmj
        do i = 1, nmi

            Kij = 0.0d0
            call fsingle_force_kernel_inout(xi(i,:,:), xj(j,:,:), ni(i), nj(j), qi(i,:,:), qj(j,:,:), &
                & acti(i,:), actj(j,:), nacti(i), nactj(j), sigma_space, r_width, c_width, Kij)

            K(:3, :3, :actmax, :actmax, j, i) = K(:3, :3, :actmax, :actmax, j, i) + Kij(:3, :3, :actmax, :actmax)

        end do
    end do

    ! do j = 1, nmj
    !     do i = 1, nmi

    !         Kij = 0.0d0
    !         call fsingle_force_kernel_inout(xi(i,:,:), xj(j,:,:), ni(i), nj(j), &
    !             & qi(i,:,:), qj(j,:,:), sigma_space, r_width, c_width, Kij)

    !         K(:3, :3, :amax, :amax, j, i) = K(:3, :3, :amax, :amax, j, i) + Kij(:3, :3, :amax, :amax)

    !     end do
    ! end do

end subroutine fforce_kernel
