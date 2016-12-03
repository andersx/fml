subroutine fgenerate_coulomb_matrix(nuclear_charges, coordinates, natoms, nmax, cm)
    double precision, dimension(:), intent(in) :: nuclear_charges
    double precision, dimension(:,:), intent(in) :: coordinates
    integer, intent(in) :: natoms
    integer, intent(in) :: nmax

    double precision, dimension((nmax + 1) * nmax / 2), intent(out):: cm
    ! double precision, dimension(nmax,nmax) :: cm
    integer :: idx

    double precision, allocatable, dimension(:) :: temp

    double precision, dimension(natoms) :: row_norms
    double precision :: pair_norm

    integer, dimension(natoms) :: sorted_atoms

    double precision  :: inv_dist
    double precision  :: row_norm

    double precision, allocatable, dimension(:,:) :: inv_distance_matrix
    
    allocate(inv_distance_matrix(natoms,natoms))
    allocate(temp(size(coordinates, dim=2)))

    ! Fill out distance matrix
    do i = 1, natoms
        do j = i+1, natoms
            temp(:) = coordinates(j,:) - coordinates(i,:)
            inv_dist = 1.0d0  / sqrt(sum(temp*temp))
            inv_distance_matrix(j, i) = inv_dist
            inv_distance_matrix(i, j) = inv_dist
        enddo
    enddo

    ! Calculate row-norms and store pair-interaction in inv_distance_matrix
    row_norms = 0.0d0
    do i = 1, natoms
        do j = i, natoms
            if (i == j) then
                pair_norm = 0.5d0 * nuclear_charges(i) ** 2.4d0
                row_norms(i) = row_norms(i) + pair_norm
                inv_distance_matrix(i, j) = pair_norm
            else
                pair_norm = nuclear_charges(i) * nuclear_charges(j) * inv_distance_matrix(i, j)
                row_norms(j) = row_norms(j) + pair_norm
                row_norms(i) = row_norms(i) + pair_norm
                inv_distance_matrix(i, j) = pair_norm
                inv_distance_matrix(j, i) = pair_norm
            endif
        enddo
    enddo

    do i = 1, natoms
        j = minloc(row_norms, dim=1)
        write(*,*) "MIN", j, row_norms(j)
        sorted_atoms(natoms - i + 1) = j
        row_norms(j) = huge(row_norms(i))
    enddo
   

    cm = 0.0d0 
    idx = 1
    do m = 1, natoms
        i = sorted_atoms(m)
        do n = 1, m
            j = sorted_atoms(n)
            cm(idx) = inv_distance_matrix(i, j)
            idx = idx + 1
        enddo
    enddo


    deallocate(temp)
    deallocate(inv_distance_matrix)
end subroutine fgenerate_coulomb_matrix
