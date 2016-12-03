subroutine fgenerate_coulomb_matrix(nuclear_charges, coordinates, nmax, cm)

    implicit none

    double precision, dimension(:), intent(in) :: nuclear_charges
    double precision, dimension(:,:), intent(in) :: coordinates
    integer :: natoms
    integer, intent(in) :: nmax

    double precision, dimension((nmax + 1) * nmax / 2), intent(out):: cm
    ! double precision, dimension(nmax,nmax) :: cm
    integer :: idx

    double precision, allocatable, dimension(:) :: row_norms
    double precision :: pair_norm

    integer, allocatable, dimension(:) :: sorted_atoms

    double precision, allocatable, dimension(:,:) :: pair_distance_matrix

    integer i, j, m, n

    if (size(coordinates, dim=1) /= size(nuclear_charges, dim=1)) then
        write(*,*) "ERROR: Coulomb matrix generation"
        write(*,*) size(coordinates, dim=2), "coordinates, but", &
            & size(nuclear_charges, dim=1), "atom_types!"
        stop
    else 
        natoms = size(nuclear_charges, dim=1)
    endif

    ! Allocate temporary     
    allocate(pair_distance_matrix(natoms,natoms))
    allocate(row_norms(natoms))
    allocate(sorted_atoms(natoms))


    ! Calculate row-norms and store pair-distances in pair_distance_matrix
    row_norms = 0.0d0

!$OMP PARALLEL DO PRIVATE(pair_norm) REDUCTION(+:row_norms)
    do i = 1, natoms
        pair_norm = 0.5d0 * nuclear_charges(i) ** 2.4d0
        row_norms(i) = row_norms(i) + pair_norm
        pair_distance_matrix(i, i) = pair_norm
    enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(pair_norm) REDUCTION(+:row_norms)
    do i = 1, natoms
        do j = i+1, natoms
            pair_norm = nuclear_charges(i) * nuclear_charges(j) &
                & / sqrt(sum((coordinates(j,:) - coordinates(i,:))**2))

            row_norms(j) = row_norms(j) + pair_norm
            row_norms(i) = row_norms(i) + pair_norm
            pair_distance_matrix(i, j) = pair_norm
            pair_distance_matrix(j, i) = pair_norm
        enddo
    enddo
!$OMP END PARALLEL DO

    !Generate sorted list of atom ids by row_norms - not really (easily) parallelizable
    do i = 1, natoms
        j = minloc(row_norms, dim=1)
        sorted_atoms(natoms - i + 1) = j
        row_norms(j) = huge(row_norms(i))
    enddo
   
    ! Fill coulomb matrix according to sorted row-norms
    cm = 0.0d0 
!$OMP PARALLEL DO PRIVATE(idx, i, j)
    do m = 1, natoms
        i = sorted_atoms(m)
        idx = (m*m+m)/2 - m
        do n = 1, m
            j = sorted_atoms(n)
            cm(idx+n) = pair_distance_matrix(i, j)
        enddo
    enddo
!$OMP END PARALLEL DO

    ! Clean up
    deallocate(pair_distance_matrix)
    deallocate(row_norms)
    deallocate(sorted_atoms)
end subroutine fgenerate_coulomb_matrix


subroutine fgenerate_local_coulomb_matrix(nuclear_charges, coordinates, natoms, nmax, cm)

    implicit none

    double precision, dimension(:), intent(in) :: nuclear_charges
    double precision, dimension(:,:), intent(in) :: coordinates
    integer,intent(in) :: natoms
    integer, intent(in) :: nmax

    double precision, dimension(natoms,(nmax + 1) * nmax / 2), intent(out):: cm
    ! double precision, dimension(nmax,nmax) :: cm
    integer :: idx

    double precision, allocatable, dimension(:) :: row_norms
    double precision :: pair_norm
    double precision :: huge_double

    integer, allocatable, dimension(:) :: sorted_atoms
    integer, allocatable, dimension(:,:) :: sorted_atoms_all

    double precision, allocatable, dimension(:,:) :: pair_distance_matrix

    integer i, j, m, n, k

    if (size(coordinates, dim=1) /= size(nuclear_charges, dim=1)) then
        write(*,*) "ERROR: Coulomb matrix generation"
        write(*,*) size(coordinates, dim=2), "coordinates, but", &
            & size(nuclear_charges, dim=1), "atom_types!"
        stop
    ! else 
    !    natoms = size(nuclear_charges, dim=1)
    endif

    huge_double = huge(1.0d0)

    ! Allocate temporary     
    allocate(pair_distance_matrix(natoms,natoms))
    allocate(row_norms(natoms))
    allocate(sorted_atoms(natoms))
    allocate(sorted_atoms_all(natoms,natoms))

    ! Calculate row-norms and store pair-distances in pair_distance_matrix
    row_norms = 0.0d0

!$OMP PARALLEL DO PRIVATE(pair_norm) REDUCTION(+:row_norms)
    do i = 1, natoms
        pair_norm = 0.5d0 * nuclear_charges(i) ** 2.4d0
        row_norms(i) = row_norms(i) + pair_norm
        pair_distance_matrix(i, i) = pair_norm
    enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(pair_norm) REDUCTION(+:row_norms)
    do i = 1, natoms
        do j = i+1, natoms
            pair_norm = nuclear_charges(i) * nuclear_charges(j) &
                & / sqrt(sum((coordinates(j,:) - coordinates(i,:))**2))

            row_norms(j) = row_norms(j) + pair_norm
            row_norms(i) = row_norms(i) + pair_norm
            pair_distance_matrix(i, j) = pair_norm
            pair_distance_matrix(j, i) = pair_norm
        enddo
    enddo
!$OMP END PARALLEL DO

    !Generate sorted list of atom ids by row_norms - not really (easily) parallelizable
    do i = 1, natoms
        j = minloc(row_norms, dim=1)
        sorted_atoms(natoms - i + 1) = j
        row_norms(j) = huge(row_norms(i))
    enddo
  
    do k = 1, natoms
        sorted_atoms_all(1, k)  = k
        m = 0
        do i = 1, natoms - 1
            if (sorted_atoms(i) == k) then
                m = 1
            endif
            sorted_atoms_all(i+1, k) = sorted_atoms(i + m)
        enddo
    enddo

    ! Fill coulomb matrix according to sorted row-norms
    cm = 0.0d0 
!$OMP PARALLEL DO PRIVATE(idx, i, j)
    do k = 1, natoms
        do m = 1, natoms
            i = sorted_atoms_all(m, k)
            idx = (m*m+m)/2 - m
            do n = 1, m
                j = sorted_atoms_all(n, k)
                cm(k, idx+n) = pair_distance_matrix(i, j)
            enddo
        enddo
    enddo
!$OMP END PARALLEL DO

    ! Clean up
    deallocate(pair_distance_matrix)
    deallocate(row_norms)
    deallocate(sorted_atoms)
    deallocate(sorted_atoms_all)

end subroutine fgenerate_local_coulomb_matrix


subroutine fgenerate_atomic_coulomb_matrix(nuclear_charges, coordinates, natoms, nmax, cm)

    implicit none

    double precision, dimension(:), intent(in) :: nuclear_charges
    double precision, dimension(:,:), intent(in) :: coordinates
    integer,intent(in) :: natoms
    integer, intent(in) :: nmax

    double precision, dimension(natoms,(nmax + 1) * nmax / 2), intent(out):: cm
    ! double precision, dimension(nmax,nmax) :: cm
    integer :: idx

    double precision :: pair_norm
    double precision :: norm
    double precision :: huge_double

    integer, allocatable, dimension(:) :: sorted_atoms
    integer, allocatable, dimension(:,:) :: sorted_atoms_all

    double precision, allocatable, dimension(:,:) :: pair_distance_matrix
    double precision, allocatable, dimension(:,:) :: distance_matrix

    integer i, j, m, n, k

    if (size(coordinates, dim=1) /= size(nuclear_charges, dim=1)) then
        write(*,*) "ERROR: Coulomb matrix generation"
        write(*,*) size(coordinates, dim=2), "coordinates, but", &
            & size(nuclear_charges, dim=1), "atom_types!"
        stop
    ! else 
    !    natoms = size(nuclear_charges, dim=1)
    endif

    huge_double = huge(1.0d0)

    ! Allocate temporary     
    allocate(distance_matrix(natoms,natoms))
    allocate(pair_distance_matrix(natoms,natoms))
    allocate(sorted_atoms(natoms))
    allocate(sorted_atoms_all(natoms,natoms))

!$OMP PARALLEL DO
    do i = 1, natoms
        distance_matrix(i, i) = 0.0d0
    enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(norm)
    do i = 1, natoms
        do j = i+1, natoms
            norm = sqrt(sum((coordinates(j,:) - coordinates(i,:))**2))
            distance_matrix(i, j) = norm
            distance_matrix(j, i) = norm
        enddo
    enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(pair_norm)
    do i = 1, natoms
        do j = i+1, natoms
            pair_norm = nuclear_charges(i) * nuclear_charges(j) &
                & / distance_matrix(j, i)

            pair_distance_matrix(i, j) = pair_norm
            pair_distance_matrix(j, i) = pair_norm
        enddo
    enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
    do i = 1, natoms
        pair_distance_matrix(i, i) = 0.5d0 * nuclear_charges(i) ** 2.4d0
    enddo
!$OMP END PARALLEL DO

    !Generate sorted list of atom ids by row_norms - not really (easily) parallelizable
    do k = 1, natoms
        do i = 1, natoms
            j = minloc(distance_matrix(:,k), dim=1)
            sorted_atoms_all(i, k) = j
            distance_matrix(j, k) = huge_double
        enddo
    enddo

    ! Fill coulomb matrix according to sorted row-norms
    cm = 0.0d0 

!$OMP PARALLEL DO PRIVATE(idx, i, j)
    do k = 1, natoms
        do m = 1, natoms
            i = sorted_atoms_all(m, k)
            idx = (m*m+m)/2 - m
            do n = 1, m
                j = sorted_atoms_all(n, k)
                cm(k, idx+n) = pair_distance_matrix(i, j)
            enddo
        enddo
    enddo
!$OMP END PARALLEL DO

    ! Clean up
    deallocate(pair_distance_matrix)
    deallocate(distance_matrix)
    deallocate(sorted_atoms)
    deallocate(sorted_atoms_all)

end subroutine fgenerate_atomic_coulomb_matrix
