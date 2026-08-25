module verlet3

    ! use kinds, ONLY: wp => dp
    use jpca15
    ! implicit none
    ! private
    ! public :: eudist, eudist_with_delta, compute_force, get_ser, get_delta_ser
    use,intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT

contains

    real(KIND=wp) function eudist(p1, p2) result(dist)
        use kinds, ONLY: wp => dp
        implicit none
        ! Two points
        real(KIND=wp), DIMENSION(3), intent(in) :: p1, p2
        ! Distance to return
        !real(KIND=wp) :: dist
        dist = SQRT( (p1(1) - p2(1))**2 + (p1(2) - p2(2))**2 + (p1(3) - p2(3))**2)
        !write(stderr,*) "poop"
    end function eudist

    subroutine get_ser(p1, p2, p3, ser)
        !  Returns the distances between pairs of points
        use kinds, ONLY: wp => dp
        implicit none
        ! Three points
        real(KIND=wp), DIMENSION(3), intent(in) :: p1, p2, p3
        real(KIND=wp), DIMENSION(3), intent(out) :: ser

        ! write(stderr,*) "P1 in get_ser: ", p1
        ! write(stderr,*) "P2 in get_ser: ", p2
        ! write(stderr,*) "P3 in get_ser: ", p3

        !ser(1) = eudist(p1, p2)
        ser(1) = eudist(p1, p2)
        ser(2) = eudist(p1, p3)
        ser(3) = eudist(p2, p3)
        ! write(stderr,*) "SER in get_ser: ", ser
    end subroutine get_ser

    subroutine get_delta_ser(pts, delta, ser)
        !  Returns the distances between pairs of points, with a delta added in each dimension
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

        write(stderr,*) "In get_delta_ser: Points"
        do dimn = 1, 3, 1
            write(stderr,*) "Atom:", dimn, ": ", pts(dimn, :)
        end do
        do dimn = 1, 3, 1
            ser(1, dimn) = eudist_with_delta(pts(1, :), pts(2, :), delta, dimn)
            ! write(stderr,*) "ser 1:", dimn, ": ", ser(1, dimn)
            ser(2, dimn) = eudist_with_delta(pts(1, :), pts(3, :), delta, dimn)
            ! write(stderr,*) "ser 2:", dimn, ": ", ser(2, dimn)
            ser(3, dimn) = eudist_with_delta(pts(2, :), pts(3, :), delta, dimn)
            ! write(stderr,*) "ser 3:", dimn, ": ", ser(3, dimn)
        end do
    end subroutine get_delta_ser

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

        ! write(stderr,*) "   -------------"
        ! write(stderr,*) "   ******** p1+:", p1_delta
        ! write(stderr,*) "   ******** p2:", p2
        ! write(stderr,*) "   ******** dist:", dist
        ! write(stderr,*) "   -------------"
        ! write(stderr,*) "eudist_with_delta: Dim:", "p1:", p1(dimn), "p1_d:", p1_delta(dimn)

    end function eudist_with_delta

    ! points = x,y,z coords of each point
    subroutine compute_force(points, delta, forces)
      ! vars
        use kinds, ONLY: wp => dp
        implicit none
        ! IN
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
        integer :: dimn, atom_i

      ! meat
      ! write(stderr,*) "----- START compute_forces"
        ! Step 1 of calculating forces: - get euclidean distances between points
      ! write(stderr,*) "points: "
        do atom_i = 1, 3, 1
          ! write(stderr,*) "    ", atom_i, points(atom_i, :)
        end do
        call get_ser(points(1, :), points(2, :), points(3, :), d_)
        ! write(stderr,*) "    DIST: ", d_
        ! Step 2 of calculating forces: - pass distances to comp_pe as `ser`
        ser = (/d_(1), d_(2), d_(3)/)
        ! write(stderr,*) "    SER: ", ser
        ! write(stderr,*) "    distances: ", d_
        call comp_pe(ser, er, der)
        ! write(stderr,*) "    er:", er

        ! ! Step 3 of calculating forces: - get euclidean distances between points + delta
        ! do dimn = 1, 3, 1
        !     delta_d_(1, dimn) = eudist_with_delta(points(1, :), points(2, :), delta, dimn)
        !   ! write(stderr,*) "Distance between A & B in", dimn,":", delta_d_(1, dimn)
        !     delta_d_(2, dimn) = eudist_with_delta(points(1, :), points(3, :), delta, dimn)
        !   ! write(stderr,*) "Distance between A & C in", dimn,":", delta_d_(2, dimn)
        !     delta_d_(3, dimn) = eudist_with_delta(points(2, :), points(3, :), delta, dimn)
        !   ! write(stderr,*) "Distance between B & C in", dimn,":", delta_d_(3, dimn)
        ! end do

        ! Atom A-B xyz
        delta_d_(1, 1) = eudist_with_delta(points(1, :), points(2, :), delta, 1)
        delta_d_(1, 2) = eudist_with_delta(points(1, :), points(2, :), delta, 2)
        delta_d_(1, 3) = eudist_with_delta(points(1, :), points(2, :), delta, 3)
        ! Atom A-C xyz
        delta_d_(2, 1) = eudist_with_delta(points(1, :), points(3, :), delta, 1)
        delta_d_(2, 2) = eudist_with_delta(points(1, :), points(3, :), delta, 2)
        delta_d_(2, 3) = eudist_with_delta(points(1, :), points(3, :), delta, 3)
        ! Atom B-A xyz
        delta_d_(3, 1) = eudist_with_delta(points(2, :), points(1, :), delta, 1)
        delta_d_(3, 2) = eudist_with_delta(points(2, :), points(1, :), delta, 2)
        delta_d_(3, 3) = eudist_with_delta(points(2, :), points(1, :), delta, 3)
        ! Atom B-C xyz
        delta_d_(4, 1) = eudist_with_delta(points(2, :), points(3, :), delta, 1)
        delta_d_(4, 2) = eudist_with_delta(points(2, :), points(3, :), delta, 2)
        delta_d_(4, 3) = eudist_with_delta(points(2, :), points(3, :), delta, 3)
        ! Atom C-A xyz
        delta_d_(5, 1) = eudist_with_delta(points(3, :), points(1, :), delta, 1)
        delta_d_(5, 2) = eudist_with_delta(points(3, :), points(1, :), delta, 2)
        delta_d_(5, 3) = eudist_with_delta(points(3, :), points(1, :), delta, 3)
        ! Atom C-B xyz
        delta_d_(6, 1) = eudist_with_delta(points(3, :), points(2, :), delta, 1)
        delta_d_(6, 2) = eudist_with_delta(points(3, :), points(2, :), delta, 2)
        delta_d_(6, 3) = eudist_with_delta(points(3, :), points(2, :), delta, 3)

        do atom_i = 1, 6, 1
         ! write(stderr,*) "Delta D: ", atom_i, delta_d_(atom_i, :)
        end do

        do atom_i = 1, 3, 1
         ! write(stderr,*) "      D: ", atom_i, d_(atom_i)
        end do

        ! Step 3.5: I don't think I need to do this, but let's initialize the forces to 0
        do dimn = 1, 3, 1
            do atom_i = 1, 3, 1
                forces(atom_i, dimn) = 0.0_wp
            end do
        end do

        ! Step 4: Update the nx3 `forces` array with the force on each atom in each dimension
    ! Atom 1 jiggled
       !write(stderr,*) "Jiggle atom in dimn", "                   Dist A-B                   Dist A-C              Dist B-C"
      ! x
        delta_ser = (/delta_d_(1, 1), delta_d_(2, 1), d_(3)/)
        call comp_pe(delta_ser, delta_er, delta_der)
        forces(1, 1) = (delta_er - er) / delta

     ! write(stderr,*) "A, x                          ", delta_ser
     ! write(stderr,*) "  delta_ser:", delta_ser
     ! write(stderr,*) "  delta_er:", delta_er
     ! write(stderr,*) "  er:", er
     ! write(stderr,*) "  forces(1, 1):", forces(1, 1)

      ! y
        delta_ser = (/delta_d_(1, 2), delta_d_(2, 2), d_(3)/)
        call comp_pe(delta_ser, delta_er, delta_der)
        forces(1, 2) = (delta_er - er) / delta

     ! write(stderr,*) "A, y                          ", delta_ser
     ! write(stderr,*) "  delta_ser:", delta_ser
     ! write(stderr,*) "  delta_er:", delta_er
     ! write(stderr,*) "  er:", er
     ! write(stderr,*) "  forces(1, 2):", forces(1, 2)

      ! z
        delta_ser = (/delta_d_(1, 3), delta_d_(2, 3), d_(3)/)
        call comp_pe(delta_ser, delta_er, delta_der)
        forces(1, 3) = (delta_er - er) / delta

     ! write(stderr,*) "A, z                          ", delta_ser
     ! write(stderr,*) "  delta_ser:", delta_ser
     ! write(stderr,*) "  delta_er:", delta_er
     ! write(stderr,*) "  er:", er
     ! write(stderr,*) "  forces(1, 3):", forces(1, 3)

    ! Atom 2 jiggled
        delta_ser = (/delta_d_(3, 1), d_(2), delta_d_(4, 1)/)
        call comp_pe(delta_ser, delta_er, delta_der)
        forces(2, 1) = (delta_er - er) / delta

     ! write(stderr,*) "B, x                          ", delta_ser
     ! write(stderr,*) "  delta_ser:", delta_ser
     ! write(stderr,*) "  delta_er:", delta_er
     ! write(stderr,*) "  er:", er
     ! write(stderr,*) "  forces(2, 1):", forces(2, 1)

        delta_ser = (/delta_d_(3, 2), d_(2), delta_d_(4, 2)/)
        call comp_pe(delta_ser, delta_er, delta_der)
        forces(2, 2) = (delta_er - er) / delta

     ! write(stderr,*) "B, y                          ", delta_ser
     ! write(stderr,*) "  delta_ser:", delta_ser
     ! write(stderr,*) "  delta_er:", delta_er
     ! write(stderr,*) "  er:", er
     ! write(stderr,*) "  forces(2, 2):", forces(2, 2)

        delta_ser = (/delta_d_(3, 3), d_(2), delta_d_(4, 3)/)
        call comp_pe(delta_ser, delta_er, delta_der)
        forces(2, 3) = (delta_er - er) / delta

     ! write(stderr,*) "B, z                          ", delta_ser
     ! write(stderr,*) "  delta_ser:", delta_ser
     ! write(stderr,*) "  delta_er:", delta_er
     ! write(stderr,*) "  er:", er
     ! write(stderr,*) "  forces(2, 3):", forces(2, 3)

    ! Atom 3 jiggled
      ! x
        delta_ser = (/d_(1), delta_d_(5, 1), delta_d_(6, 1)/)
        call comp_pe(delta_ser, delta_er, delta_der)
        forces(3, 1) = (delta_er - er) / delta

     ! write(stderr,*) "C, x                          ", delta_ser
     ! write(stderr,*) "  delta_ser:", delta_ser
     ! write(stderr,*) "  delta_er:", delta_er
     ! write(stderr,*) "  er:", er
     ! write(stderr,*) "  forces(3, 1):", forces(3, 1)

      ! y
        delta_ser = (/d_(1), delta_d_(5, 2), delta_d_(6, 2)/)
        call comp_pe(delta_ser, delta_er, delta_der)
        forces(3, 2) = (delta_er - er) / delta

     ! write(stderr,*) "C, y                          ", delta_ser
     ! write(stderr,*) "  delta_ser:", delta_ser
     ! write(stderr,*) "  delta_er:", delta_er
     ! write(stderr,*) "  er:", er
     ! write(stderr,*) "  forces(3, 2):", forces(3, 2)

      ! z
        delta_ser = (/d_(1), delta_d_(5, 3), delta_d_(6, 3)/)
        call comp_pe(delta_ser, delta_er, delta_der)
        forces(3, 3) = (delta_er - er) / delta

     ! write(stderr,*) "C, z                          ", delta_ser
     ! write(stderr,*) "  delta_ser:", delta_ser
     ! write(stderr,*) "  delta_er:", delta_er
     ! write(stderr,*) "  er:", er
     ! write(stderr,*) "  forces(3, 3):", forces(3, 3)

        ! do atom_i = 1, 3, 1
       !  !write(stderr,*) "Atom:", atom_i, forces(atom_i, :)
        ! end do

        ! DONE!
      ! write(stderr,*) "----- END compute_forces"
    end subroutine compute_force

end module verlet3
