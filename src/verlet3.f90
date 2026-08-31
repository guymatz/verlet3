module verlet3

    use kinds, ONLY: wp => dp
    use jpca15
    use stdlib_logger
    implicit none
    !  CAN I HAVE PUBLIC / PRIVATE with tests?!
    ! private :: eudist, eudist_with_delta, get_ser, get_delta_ser
    ! public :: verlet3_init, compute_force

    character(len=100) :: log_msg

contains

    ! Mostly for configuring logging
    subroutine verlet3_init(verbose)
        implicit none
        logical, intent(in) :: verbose
        call global_logger%configure(indent=.true., max_width=100)
        if (verbose) then
            call global_logger%configure(level = ALL_LEVEL)
        else
            call global_logger%configure(level = NONE_LEVEL)
        end if
        call global_logger%log_information("Verlet has been Initialized")
    end subroutine verlet3_init

    ! Calculates Euclidean distance between two 3D points
    real(KIND=wp) function eudist(p1, p2) result(dist)
        use kinds, ONLY: wp => dp
        implicit none
        ! Two points
        real(KIND=wp), DIMENSION(3), intent(in) :: p1, p2
        ! Distance to return
        !real(KIND=wp) :: dist
        dist = SQRT( (p1(1) - p2(1))**2 + (p1(2) - p2(2))**2 + (p1(3) - p2(3))**2)
    end function eudist

    !  Returns the distances between pairs of 3 points
    subroutine get_ser(p1, p2, p3, ser)
        use kinds, ONLY: wp => dp
        implicit none
        ! Three points
        real(KIND=wp), DIMENSION(3), intent(in) :: p1, p2, p3
        real(KIND=wp), DIMENSION(3), intent(out) :: ser

        write(log_msg, '(A)') "==== In get_ser . . . "
        call global_logger%log_warning(log_msg)

        write(log_msg, '(A, 3F12.5)') "P1 in get_ser:  ", p1
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F12.5)') "P2 in get_ser:  ", p2
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F12.5)') "P3 in get_ser:  ", p3
        call global_logger%log_warning(log_msg)

        ser(1) = eudist(p1, p2)
        ser(2) = eudist(p1, p3)
        ser(3) = eudist(p2, p3)
        write(log_msg, '(A, 3F12.5)') "SER in get_ser: ", ser
        call global_logger%log_warning(log_msg)
    end subroutine get_ser

    ! Returns the distances between pairs of 3 points, with a delta added in each dimension
    subroutine get_delta_ser(pts, delta, ser)
        use kinds, ONLY: wp => dp
        implicit none
        ! points
        real(KIND=wp), DIMENSION(3, 3), intent(in) :: pts
        real(KIND=wp), intent(in) :: delta
        ! The atom in the list to perturb
        ! Return n x d (3) vector
        real(KIND=wp), DIMENSION(3, 3), intent(out) :: ser
        ! for looping
        integer :: dimn

        do dimn = 1, 3, 1
            ser(1, dimn) = eudist_with_delta(pts(1, :), pts(2, :), delta, dimn)
            ser(2, dimn) = eudist_with_delta(pts(1, :), pts(3, :), delta, dimn)
            ser(3, dimn) = eudist_with_delta(pts(2, :), pts(3, :), delta, dimn)
        end do
    end subroutine get_delta_ser

    ! Returns the distance between a pair of points, with a delta added in each dimension
    real(KIND=wp) function eudist_with_delta(p1, p2, delta, dimn) result(dist)
        ! delta_d_(1, dimn) = eudist_with_delta(points(1, :), points(2, :), delta, dimn)
        use kinds, ONLY: wp => dp
        implicit none
        ! Two points
        real(KIND=wp), DIMENSION(3), intent(in) :: p1, p2
        ! delta to add to dimension, and distance to return
        real(KIND=wp), intent(in) :: delta
        ! Dimension which to add delta
        integer, intent(in) :: dimn
        ! New point with delta added to proper dimension
        real(KIND=wp), DIMENSION(3) :: p1_delta
        p1_delta = p1
        p1_delta(dimn) = p1_delta(dimn) + delta
        ! Distance to return
        !!!!!   NOT Adding delta to second atom.  Only the first
        dist = eudist(p1_delta, p2)
    end function eudist_with_delta

    function fact(n) result(f)
        integer :: i, n, f
        do i = 1, n
            f = f + i
        end do
    end function fact

    subroutine get_all_eudists_with_delta(points, delta, dists)
        integer :: atom_pair_ctr, atom_1, atom_2, dimn
        real(KIND=wp), intent(in) :: delta
        real(KIND=wp), dimension(3, 3), intent(in) :: points
        real(KIND=wp), dimension(6, 3), intent(out) :: dists
        real(KIND=wp) :: tmp

        write(log_msg, '(A)') "IN get_all_eudists_with_delta"
        call global_logger%log_warning(log_msg)
        
        atom_pair_ctr = 0
        do atom_1 = 1, size(points, 1)
            do atom_2 = 1, size(points, 1)
                if (atom_1 == atom_2) then
                    ! only compute distances for pairs
                    cycle
                end if
                atom_pair_ctr = atom_pair_ctr + 1
                do dimn = 1, size(points, 2)
                    dists(atom_pair_ctr, dimn) = eudist_with_delta(points(atom_1, :), points(atom_2, :), delta, dimn)
                end do
            end do
        end do
    end subroutine get_all_eudists_with_delta

    ! Uses Velocity Verlet to compute the force on each atom and in each x,y,z direction
    subroutine compute_force(points, delta, forces)
      ! vars
        use kinds, ONLY: wp => dp
        implicit none
        ! IN
        ! points = x,y,z coords of each point
        real(KIND=wp), DIMENSION(3, 3), intent(in) :: points
        real(KIND=wp), intent(in) :: delta
        ! OUT
        real(KIND=wp), DIMENSION(3, 3), intent(out) :: forces
        ! We will compute `ser` for input to comp_pe, and stores the results in der
        real(KIND=wp), DIMENSION(3) :: ser, der
        ! comp_pe also returns er
        real(KIND=wp) :: er
        ! We will compute `delta_ser` for input to comp_pe, and stores the results in delta_der
        real(KIND=wp), DIMENSION(3) :: delta_ser, delta_der
        ! comp_pe also returns delta_er
        real(KIND=wp) :: delta_er
        ! d_ is used to store the euclidean distances between atoms
        real(KIND=wp), DIMENSION(3) :: d_
        ! delta_d_ is used to store the euclidean distances + delta between atoms
        ! It's 2D so that we can store distances between atoms in all dimensions
        real(KIND=wp), DIMENSION(6, 3) :: delta_d_
        ! for looping
        integer :: dimn, atom_i, atom_pair_ctr, atom_1, atom_2

      ! meat
        write(log_msg, '(A, F5.1)') "----- START compute_forces"
        call global_logger%log_warning(log_msg)

      ! Step 1 of calculating forces: - get euclidean distances between points
        call get_ser(points(1, :), points(2, :), points(3, :), d_)
        write(log_msg, '(A, 3F5.1)') "    DIST: ", d_
        call global_logger%log_warning(log_msg)

      ! Step 2 of calculating forces: - pass distances to comp_pe as `ser`
        ser = (/d_(1), d_(2), d_(3)/)
        write(log_msg, '(A, 3F5.1)') "    SER: ", ser
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F5.1)') "    distances: ", d_
        call global_logger%log_warning(log_msg)
        call comp_pe(ser, er, der)
        write(log_msg, '(A, F15.9)') "    er:", er
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F15.9)') "    der:", der
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.9)') "    sum der:", sum(der)
        call global_logger%log_warning(log_msg)

      ! Step 3 of calculating forces: - get euclidean distances between points + delta
        ! Loop through axes of pairs of points to get all Euclidean distances between pairs,
        ! with a delta added in each dimension

        do atom_pair_ctr = 1, size(delta_d_, 1)
            delta_d_(atom_pair_ctr, :) = (/0.0, 0.0, 0.0/)
        end do

        call get_all_eudists_with_delta(points, delta, delta_d_)

        write(log_msg, '(A)') "JUST worked on delta_d_"
        call global_logger%log_warning(log_msg)

        atom_pair_ctr = 0
        do atom_1 = 1, size(points, 1)
            atom_pair_ctr = atom_pair_ctr + 1
            do atom_2 = 1, size(points, 1)
                if (atom_1 == atom_2) then
                    cycle
                end if
                do dimn = 1, size(points, 2)
                    write(log_msg, '(A, I0, A, I0, A, I0, A, F20.10)') "Atom 1 & 2 - dimn: ", atom_1, " - ", atom_2, " - ", dimn, ": ", delta_d_(atom_pair_ctr, dimn) 
                    call global_logger%log_warning(log_msg)
                end do
            end do
        end do

        do atom_i = 1, 6
            write(log_msg, '(A, I0, 3F10.5)') "Delta D 2: ", atom_i, delta_d_(atom_i, :)
            call global_logger%log_warning(log_msg)
        end do

        do atom_i = 1, 3
            write(log_msg, '(A, I0, 3F5.1)') "      D: ", atom_i, d_(atom_i)
            call global_logger%log_warning(log_msg)
        end do

        ! Step 3.5: I don't think I need to do this, but let's initialize the forces to 0
        do dimn = 1, 3, 1
            do atom_i = 1, 3, 1
                forces(atom_i, dimn) = 0.0_wp
            end do
        end do

  ! Step 4: Update the nx3 `forces` array with the force on each atom in each dimension
      ! Atom 1 jiggled
        write(log_msg, '(A, A40)') "Jiggle atom in dimn", "A-B      A-C      B-C"
        call global_logger%log_warning(log_msg)
    ! x
        delta_ser = (/delta_d_(1, 1), delta_d_(2, 1), d_(3)/)
        call comp_pe(delta_ser, delta_er, delta_der)
        forces(1, 1) = - (delta_er - er) / delta

        write(log_msg, '(A, T25, F15.9, F15.9, F15.9)') "Ax", delta_ser(1), delta_ser(2), delta_ser(3)
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F15.10)') "  delta_ser:", delta_ser
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F15.10)') "  delta_er:", delta_er
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F15.10)') "  er:", er
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F15.10)') " delta_der:", delta_der
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') " sum delta_der:", sum(delta_der)
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  forces(1, 1):", forces(1, 1)
        call global_logger%log_warning(log_msg)

    ! y
        delta_ser = (/delta_d_(1, 2), delta_d_(2, 2), d_(3)/)
        call comp_pe(delta_ser, delta_er, delta_der)
        forces(1, 2) = - (delta_er - er) / delta

        write(log_msg, '(A, T35, F7.2, F10.2, F8.2)') "Ay", delta_ser(1), delta_ser(2), delta_ser(3)
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F15.10)') "  delta_ser:", delta_ser
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  delta_er:", delta_er
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  er:", er
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F15.10)') " delta_der:", delta_der
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') " sum delta_der:", sum(delta_der)
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  forces(1, 2):", forces(1, 2)
        call global_logger%log_warning(log_msg)

    ! z
        delta_ser = (/delta_d_(1, 3), delta_d_(2, 3), d_(3)/)
        call comp_pe(delta_ser, delta_er, delta_der)
        forces(1, 3) = - (delta_er - er) / delta

        write(log_msg, '(A, T35, F7.2, F10.2, F8.2)') "Az", delta_ser(1), delta_ser(2), delta_ser(3)
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F15.10)') "  delta_ser:", delta_ser
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  delta_er:", delta_er
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  er:", er
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F15.10)') " delta_der:", delta_der
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') " sum delta_der:", sum(delta_der)
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  forces(1, 3):", forces(1, 3)
        call global_logger%log_warning(log_msg)

    ! Atom 2 jiggled
        delta_ser = (/delta_d_(3, 1), d_(2), delta_d_(4, 1)/)
        call comp_pe(delta_ser, delta_er, delta_der)
        forces(2, 1) = - (delta_er - er) / delta

        !write(log_msg, '(A, T35, F7.2, F10.2, F8.2)') "Bx", delta_ser(1), delta_ser(2), delta_ser(3)
        write(log_msg, '(A)') "Bx"
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F15.10)') "  delta_ser:", delta_ser
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F15.10)') "  delta_der:", delta_der
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  delta_er:", delta_er
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  er:", er
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F15.10)') " delta_der:", delta_der
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') " sum delta_der:", sum(delta_der)
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  forces(2, 1):", forces(2, 1)
        call global_logger%log_warning(log_msg)

        delta_ser = (/delta_d_(3, 2), d_(2), delta_d_(4, 2)/)
        call comp_pe(delta_ser, delta_er, delta_der)
        forces(2, 2) = - (delta_er - er) / delta

        write(log_msg, '(A, T35, F7.2, F10.2, F8.2)') "By", delta_ser(1), delta_ser(2), delta_ser(3)
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F15.10)') "  delta_ser:", delta_ser
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  delta_er:", delta_er
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  er:", er
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F15.10)') " delta_der:", delta_der
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') " sum delta_der:", sum(delta_der)
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  forces(2, 2):", forces(2, 2)
        call global_logger%log_warning(log_msg)

        delta_ser = (/delta_d_(3, 3), d_(2), delta_d_(4, 3)/)
        call comp_pe(delta_ser, delta_er, delta_der)
        forces(2, 3) = - (delta_er - er) / delta

        write(log_msg, '(A, T35, F7.2, F10.2, F8.2)') "Bz", delta_ser(1), delta_ser(2), delta_ser(3)
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F15.10)') "  delta_ser:", delta_ser
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  delta_er:", delta_er
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  er:", er
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F15.10)') " delta_der:", delta_der
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') " sum delta_der:", sum(delta_der)
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  forces(2, 3):", forces(2, 3)
        call global_logger%log_warning(log_msg)

    ! Atom 3 jiggled
      ! x
        delta_ser = (/d_(1), delta_d_(5, 1), delta_d_(6, 1)/)
        call comp_pe(delta_ser, delta_er, delta_der)
        forces(3, 1) = - (delta_er - er) / delta

        write(log_msg, '(A, T35, F7.2, F10.2, F8.2)') "Cx", delta_ser(1), delta_ser(2), delta_ser(3)
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F15.10)') "  delta_ser:", delta_ser
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  delta_er:", delta_er
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  er:", er
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F15.10)') " delta_der:", delta_der
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') " sum delta_der:", sum(delta_der)
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  forces(3, 1):", forces(3, 1)
        call global_logger%log_warning(log_msg)

      ! y
        delta_ser = (/d_(1), delta_d_(5, 2), delta_d_(6, 2)/)
        call comp_pe(delta_ser, delta_er, delta_der)
        forces(3, 2) = - (delta_er - er) / delta

        write(log_msg, '(A, T35, F7.2, F10.2, F8.2)') "Cy", delta_ser(1), delta_ser(2), delta_ser(3)
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F15.10)') "  delta_ser:", delta_ser
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  delta_er:", delta_er
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  er:", er
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F15.10)') " delta_der:", delta_der
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') " sum delta_der:", sum(delta_der)
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  forces(3, 2):", forces(3, 2)
        call global_logger%log_warning(log_msg)

      ! z
        delta_ser = (/d_(1), delta_d_(5, 3), delta_d_(6, 3)/)
        call comp_pe(delta_ser, delta_er, delta_der)
        forces(3, 3) = - (delta_er - er) / delta

        write(log_msg, '(A, T35, F7.2, F10.2, F8.2)') "Cz", delta_ser(1), delta_ser(2), delta_ser(3)
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F15.10)') "  delta_ser:", delta_ser
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  delta_er:", delta_er
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  er:", er
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, 3F15.10)') " delta_der:", delta_der
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') " sum delta_der:", sum(delta_der)
        call global_logger%log_warning(log_msg)
        write(log_msg, '(A, F15.10)') "  forces(3, 3):", forces(3, 3)
        call global_logger%log_warning(log_msg)

        do atom_i = 1, 3, 1
            write(log_msg, '(A, I0, 3F15.10)') "Force Atom:", atom_i, forces(atom_i, :)
            call global_logger%log_warning(log_msg)
        end do

        write(log_msg, '(A, 3F15.10)') "Total Force: ", sum(forces, dim=1)
        call global_logger%log_warning(log_msg)

        ! DONE!
        write(log_msg, '(A, F15.10)') "----- END compute_forces"
        call global_logger%log_warning(log_msg)
    end subroutine compute_force

end module verlet3
