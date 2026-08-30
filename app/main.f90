program main
! 1verbose set fdm?

! Var definitions, etc
    ! For printing to log_msg
    use,intrinsic :: iso_fortran_env, only: output_unit
    !use stdlib_logger, only: global => global_logger
    use stdlib_logger
    use kinds, ONLY: wp => dp
    use verlet3, ONLY: verlet3_init, compute_force
    implicit none

! vars
    !  for logging
    character(len=100) :: log_msg
    ! for looping by atoms & dimension
    integer :: atom_num, dimn

    LOGICAL :: OK = .FALSE.
    ! For requesting xyz file output
    LOGICAL :: XYZ = .TRUE.
    CHARACTER(len=32) :: arg
    CHARACTER(len=32) :: file_name = "data/reactive.dat"
    !integer :: arg_len
    !integer :: status
    real(KIND=wp), DIMENSION(:, :), ALLOCATABLE :: x, v, f, fnext
    real(KIND=wp), DIMENSION(:), ALLOCATABLE :: mass
    ! ser, er & der are for parameters to jpca15 function
    ! INPUT
    !   ser: a vector with the three interatomic distances (AB, AC, and BC)
    ! OUTPUT
    !   er: potential energy (in eV)
    !   der: vector of the derivatives of the potential with
    !        respect to the interatomic distances (AB, AC, and BC) (in bohr)
    !real(KIND=wp), DIMENSION(3) :: d_AB_ser, d_AB_er, d_AB_der
    !real(KIND=wp), DIMENSION(3) :: d_AC_ser, d_AC_er, d_AC_der
    !real(KIND=wp), DIMENSION(3) :: d_BC_ser, d_BC_er, d_BC_der
    ! Force vectors (in the x-direction only)
    !real(KIND=wp), DIMENSION(3) :: f_AB, f_AC, f_BC, force
    ! for reading in atomic data from file
    integer :: nk, nk_cli=0
    !real :: sigma, epsilon
    real(KIND=wp) :: tau, tau_cli=0.0_wp
    !real :: tmp
    ! interatomic distances
    !real(KIND=wp) :: d_AB, d_AC, d_BC
    ! Delta for interatomic distances
    !real(KIND=wp), DIMENSION(3) :: dd_AB, dd_AC, dd_BC
    real(KIND=wp) :: delta = 0.1, delta_cli = 0
    ! For computing jPCA with delta
    !real(KIND=wp), DIMENSION(3) :: dd_ser, dd_er, dd_der
    !integer :: num_rows
    integer :: num_atoms
    ! locations and velocites for atom
    real(KIND=wp) :: ax, ay, az, vx, vy, vz
    ! for storing intermediate values & looping
    integer :: i, k
    !real (KIND = wp), DIMENSION(7) :: p_a, p_b ! our two particles

! pre-processing config
    call global_logger%configure(time_stamp=.false.)
    call global_logger%configure(indent=.true., max_width=100)
    call global_logger%configure(level = NONE_LEVEL)
    ! for logging in verlet3 module.  Start false, then modify, if needed
    call verlet3_init(.false.)

! Process command-line args
    i = 0
    DO
        i = i + 1
        IF (i .gt. command_argument_count()) exit

        CALL get_command_argument(i, arg)
        ! Delta - defaults to 0.1 (see above)
        IF (arg == "-v") THEN
            call global_logger%configure(level = ALL_LEVEL)
            ! for logging in verlet3 module.  DO I NEED THIS?
!            call verlet3_init(.true.)
        ELSE IF (arg == "-d") THEN
            ! delta
            i = i + 1
            CALL get_command_argument(i, arg)
            write(log_msg, '(A, A)') "arg -d: ", arg
            ! CALL global_logger%log_warning(log_msg)
            read (arg, '(f33.32)') delta_cli
            IF (delta_cli <= 0.0) THEN
                print *, "Delta too small!"
                STOP __LINE__ - 1
            END IF
        ELSE IF (arg == "-s") THEN
            ! Num steps - defaults to 0 (see above): Should be in data file
            i = i + 1
            CALL get_command_argument(i, arg)
            read (arg, '(I5)') nk_cli
        ELSE IF (arg == "-t") THEN
            ! tau
            i = i + 1
            CALL get_command_argument(i, arg)
            read (arg, '(f1.2)') tau_cli
        ELSE IF (arg == "-f") THEN
            ! filename
            i = i + 1
            CALL get_command_argument(i, arg)
            file_name = trim(arg)
            INQUIRE (FILE=file_name, EXIST=OK)
            if (.NOT. OK) THEN
                write(log_msg, '(A, A)') "ERROR!!  File does not exist: ", file_name
                ! CALL global_logger%log_warning(log_msg)
                STOP __LINE__ - 1
            END IF
        ELSE IF (arg == "-X") THEN
            write(log_msg, '(A)') "No output in XYZ format"
            CALL global_logger%log_warning(log_msg)
            XYZ = .FALSE.
        ELSE
            call global_logger%configure(indent=.true., max_width=110)
            call global_logger%configure(time_stamp=.false.)
            call global_logger%configure(level = ALL_LEVEL)
            write(log_msg, '(A, A)') "Unknown Arg used: ", arg
            CALL global_logger%log_error(log_msg)
            write(log_msg, '(A)') "Usage: verlet3 [ -h ] [ -f data_file ]  [ -d delta ] [ -s steps ] [-X]"
            CALL global_logger%log_error(log_msg)
            write(log_msg, '(A)') "DEFAULTS:"
            CALL global_logger%log_error(log_msg)
            write(log_msg, '(A20, A)') "data_file: ", file_name
            CALL global_logger%log_error(log_msg)
            write(log_msg, '(A20, F0.9)') "delta: ", delta
            CALL global_logger%log_error(log_msg)
            write(log_msg, '(A20, F0.9)') "tau: ", tau
            CALL global_logger%log_error(log_msg)
            write(log_msg, '(A20, I0)') "steps: ", nk
            CALL global_logger%log_error(log_msg)
            write(log_msg, '(A20, L1)') "Ouput XYZ data: ", XYZ
            CALL global_logger%log_error(log_msg)
            !write(log_msg, '(A, I0)') "main +", __LINE__ + 1
            !CALL global_logger%log_warning(log_msg)
            STOP __LINE__ - 1
        END IF
    END DO

! Read atoms, etc from data file
    open (UNIT=11, FILE=file_name, STATUS="old", ACTION="read")
    read (unit=11, FMT=*) nk, tau
    read (unit=11, FMT=*) delta
    read (unit=11, FMT=*) num_atoms
    !write(log_msg, '(A)') nk, tau, sigma, epsilon, num_atoms
    write(log_msg, '(A, I0)') "Number of atoms: ", num_atoms
    ! CALL global_logger%log_warning(log_msg)

    ! Let the command-line `steps` override nk, if it's set
    if (nk_cli > 0) then
        nk = nk_cli
    end if
    ! same for delta
    if (delta_cli > 0) then
        delta = delta_cli
    end if
    ! nd tau
    if (tau_cli > 0) then
        tau = tau_cli
    end if

    write(log_msg, '(A, F5.1)') "Will use delta: ", delta
    ! CALL global_logger%log_warning(log_msg)
    write(log_msg, '(A, A)') "Will use file: ", file_name
    ! CALL global_logger%log_warning(log_msg)
    write(log_msg, '(A, I0)') "Will use num steps: ", nk
    ! CALL global_logger%log_warning(log_msg)
    write(log_msg, '(A, F5.1)') "Will use tau: ", tau
    ! CALL global_logger%log_warning(log_msg)
    write(log_msg, '(A, L)') "Create XYZ file: ", XYZ
    ! CALL global_logger%log_warning(log_msg)

    ! Allocate arrays for position, velocity, force & mass
    ! Position
    allocate (x(num_atoms, 3))
    ! Velocity
    allocate (v(num_atoms, 3))
    ! Force
    allocate (f(num_atoms, 3))
    ! Force
    allocate (fnext(num_atoms, 3))
    ! mass
    allocate (mass(num_atoms))

    ! read in info for particles
    do i = 1, num_atoms, 1
        read (unit=11, FMT=*) mass(i), ax, ay, az, vx, vy, vz
        x(i, :) = (/ax, ay, az/)
        v(i, :) = (/vx, vy, vz/)
        write(log_msg, '(A, I0, A)') 'Particle ', i, ': '
        ! CALL global_logger%log_information(log_msg)
        write(log_msg, '(A, F5.1)') '  Mass: ', mass(i)
        ! CALL global_logger%log_information(log_msg)
        write(log_msg, '(A, 3F5.1)') '  Starting Position: ', x(i, :)
        ! CALL global_logger%log_information(log_msg)
        write(log_msg, '(A, 3F5.1)') '  Initial Velocity: ', v(i, :)
        ! CALL global_logger%log_information(log_msg)
    end do
    close (unit=11)

    write(log_msg, '(A)') "initial output for XYZ data file -"
    ! CALL global_logger%log_warning(log_msg)
    ! https://en.wikipedia.org/wiki/XYZ_file_format
    if (XYZ) then
        print *, num_atoms
        print *, "Initial Positions"
        print *, "atom1", x(1, :)
        print *, "atom2", x(2, :)
        print *, "atom3", x(3, :)
    end if

! Steps here found Chap 7, on p. 61 of ``Chemistry at the Fronteir...'' by Rampino
  ! Step 0: Calculate initial force on all particles
    ! Initial Force for particles
    call compute_force(x, delta, f)
    write(log_msg, '(A)') " ***** Initial Positions / Forces: "
    ! CALL global_logger%log_warning(log_msg)
    do i = 1, size(x, 1), 1
        write(log_msg, '(A, I0)') "  Atom: ", i
        CALL global_logger%log_information(log_msg)
        write(log_msg, '(A)') "          X    Y    Z"
        CALL global_logger%log_information(log_msg)
        write(log_msg, '(A, 3F15.9)') "     x: ", x(i, :)
        CALL global_logger%log_information(log_msg)
        write(log_msg, '(A, 3F15.9)') "     f: ", f(i, :)
        CALL global_logger%log_information(log_msg)
    end do

  ! Iterate!  For nk # of steps
    do k = 1, nk, 1
        do atom_num = 1, 3, 1
            write(log_msg, '(A, I0, A, 3F5.1)') " ***** x before: ", atom_num, " - ", x(atom_num, :) 
            ! CALL global_logger%log_warning(log_msg)
        end do
        do atom_num = 1, size(x, 1), 1
            do dimn = 1, size(x, 2), 1
                write(log_msg, '(A)') "NUMS: "
                ! CALL global_logger%log_warning(log_msg)
                write(log_msg, '(A, F15.10)') "      x: ", x(atom_num, dimn)
                ! CALL global_logger%log_warning(log_msg)
                write(log_msg, '(A, F15.10)') "      t: ", tau
                ! CALL global_logger%log_warning(log_msg)
                write(log_msg, '(A, F15.10)') "      v: ", v(atom_num, dimn)
                ! CALL global_logger%log_warning(log_msg)
                write(log_msg, '(A, F15.10)') "      f: ", f(atom_num, dimn)
                ! CALL global_logger%log_warning(log_msg)
                write(log_msg, '(A, F15.10)') "      m: ", mass(atom_num)
                ! CALL global_logger%log_warning(log_msg)
                write(log_msg, '(A, F15.10)') "     l1: ", tau * v(atom_num, dimn)
                ! CALL global_logger%log_warning(log_msg)
                write(log_msg, '(A, F15.10)') "     l2: ", x(atom_num, dimn) + tau * v(atom_num, dimn)
                ! CALL global_logger%log_warning(log_msg)
                write(log_msg, '(A, F15.10)') "   l3-1: ", f(atom_num, dimn)
                ! CALL global_logger%log_warning(log_msg)
                write(log_msg, '(A, F15.10)') "   l3-2: ", 2 * mass(atom_num)
                ! CALL global_logger%log_warning(log_msg)
                write(log_msg, '(A, F15.10)') "     l3: ", f(atom_num, dimn) / (2 * mass(atom_num)) * tau**2
                ! CALL global_logger%log_warning(log_msg)

  ! Step 1: Calculate x_{k+1}^{particle}
                x(atom_num, dimn) = x(atom_num, dimn) + tau * v(atom_num, dimn) + &
                                     tau**2 * (f(atom_num, dimn) / (2 * mass(atom_num)))
                write(log_msg, '(A, I0, A, I0, A, F5.1)') "  new x: ", atom_num, ", ", dimn, ", ", x(atom_num, dimn)
                ! CALL global_logger%log_warning(log_msg)
            end do
        end do
        do atom_num = 1, size(x, 1)
            write(log_msg, '(A, I0, 3F5.1)') " ***** x AFTER: ", atom_num, x(atom_num, :) 
            ! CALL global_logger%log_warning(log_msg)
        end do
  ! Step 2: Calculate new force for each dimension - f_{k+1}^{particle, dimension}
        call compute_force(x, delta, fnext)
        do atom_num = 1, size(x, 1)
            write(log_msg, '(A, I0, 3F15.9)') "fnext: ", atom_num, fnext(atom_num, :) 
            ! CALL global_logger%log_warning(log_msg)
        end do
  ! Step 3: Calculate velocity for each dimension - v_{k+1}^{particle, dimension}
        ! A
        do atom_num = 1, size(x, 1)
            do dimn = 1, size(x, 2)
                v(atom_num, dimn) = v(atom_num, dimn) + &
                                    (tau / (2 * mass(atom_num))) * &
                                    (f(atom_num, dimn) + fnext(atom_num, dimn))
            end do
        end do

        ! Print out the coordinates in XYZ format, if requested (with -x)
        if (XYZ) then
            print *, num_atoms
            print *, "step:", k
            do atom_num = 1, size(x, 1)
                write(*, '(A, I0, A5, 3F20.10)') "atom", atom_num, "", x(atom_num, :)
            end do
        end if

  ! Step 4: Assign the value of f_{k+1}^{particle} to f_{k}^{particle}
        do atom_num = 1, 3, 1
          do dimn = 1, 3, 1
            f(atom_num, dimn) = fnext(atom_num, dimn)
            write(log_msg, '(I5, I5, F20.15)') atom_num, dimn, f(atom_num, dimn) 
            ! CALL global_logger%log_warning(log_msg)
          end do
        end do

    end do

! Cleanup
    deallocate (x)
    deallocate (v)
    deallocate (f)
    deallocate (fnext)
    deallocate (mass)

end program main
