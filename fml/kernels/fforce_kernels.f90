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

module tools

    contains
function outer3d(ri, rj)

    double precision, dimension(3), intent(in) :: ri, rj
    double precision, dimension(3,3) :: outer3d

    outer3d(:3,:3) = spread(ri(:3), dim=2, ncopies=3) * spread(rj(:3), dim=1, ncopies=3)

end function outer3d
end module tools

! subroutine fsingle_force_kernel(xi, xj, ni, nj, qi, qj, sigma_space, sigma_periodic, K)
subroutine fsingle_force_kernel(xi, xj, ni, nj, sigma_space, sigma_periodic, K)

    use tools, only: outer3d 

    ! Coordinates, molecule i
    double precision, dimension(:,:), intent(in) :: xi

    ! Coordinates, molecule j
    double precision, dimension(:,:), intent(in) :: xj

    ! Number of atoms, molecule i
    integer, intent(in) :: ni

    ! Number of atoms, molecule j
    integer, intent(in) :: nj

    ! ! Periodic table,  molecule i
    ! double precision, dimension(:,:), intent(in) :: qi

    ! ! Periodic table, molecule j
    ! double precision, dimension(:,:), intent(in) :: qj

    ! Kernel width for real-space and periodic table
    double precision, intent(in) :: sigma_space
    double precision, intent(in) :: sigma_periodic

    ! The actual kernel
    double precision, dimension(ni, nj, 3, 3), intent(out) :: K

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

    K(:,:,:,:) = 0.0d0

    inv_L = 1.0d0 / (2.0d0 * sqrt(pi * sigma_space**2)**3)

    do i = 1, ni
        do j = 1, nj
            do ii = 1, ni
                do jj = 1, nj

                    if ((i == ii).or.(j == jj)) cycle

                    ri = xi(i,:) - xi(ii,:)
                    rj = xj(j,:) - xj(jj,:)

                    ri_norm = norm2(ri)
                    rj_norm = norm2(rj)

                    alpha_ij = (ri_norm**2 + rj_norm**2) / (4.0d0 * sigma_space**2)
                    gamma_ij = (ri_norm * rj_norm) / (2.0d0 * sigma_space**2)

                    phi_ij = exp(-alpha_ij) / (gamma_ij**2) * (gamma_ij * cosh(gamma_ij) - sinh(gamma_ij))

                    K(i,j,:3,:3) = K(i,j,:3,:3) + outer3d(ri, rj) * phi_ij

                end do
            end do
        end do
    end do

    K(:,:,:,:) = K(:,:,:,:) * inv_L

end subroutine fsingle_force_kernel
