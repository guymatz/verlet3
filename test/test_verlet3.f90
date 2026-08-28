module test_verlet3
  use stdlib_logger
  use verlet3
  use kinds, ONLY: wp => dp
  use testdrive, only : error_type, unittest_type, new_unittest, check
  implicit none
  private

  public :: collect_verlet3
  type :: test_fixture_type
    character(len=100) :: log_msg
    real(KIND=wp), DIMENSION(3, 3) :: x, v, f, fnext, points, forces
    real(KIND=wp), DIMENSION(3, 3) :: expected_forces
    real(KIND=wp), DIMENSION(3) :: mass
    real(KIND=wp), DIMENSION(3) :: ser, expected_get_ser
    real(KIND=wp), DIMENSION(3, 3) :: ser_delta, expected_get_ser_delta
    real(KIND=wp) :: distance_with_delta
    real(KIND=wp) :: tau, delta
    real(KIND=wp) :: tol ! tolerance for assertEqual tests
    real(KIND=wp) :: real_result
    real(KIND=wp), DIMENSION(3) :: real_result_3d
    real(KIND=wp), DIMENSION(3, 3) :: real_result_3x3
    real(KIND=wp) :: expected_eudist

    ! test distance between atoms 1 & 2 in dimension 1, etc
    real(KIND=wp) :: expected_eudist_with_delta_12_1
    real(KIND=wp) :: expected_eudist_with_delta_13_2
    real(KIND=wp) :: expected_eudist_with_delta_23_3
  end type test_fixture_type

contains

  subroutine collect_verlet3(testsuite)
    !> Collection of tests

    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
                  new_unittest("eudist_p1_p2", test_eudist_p1_p2), &
                  new_unittest("eudist_p2_p1", test_eudist_p2_p1), &
                  new_unittest("eudist_with_delta_12_1", test_eudist_with_delta_12_1), &
                  new_unittest("eudist_with_delta_13_2", test_eudist_with_delta_13_2), &
                  new_unittest("eudist_with_delta_23_3", test_eudist_with_delta_23_3), &
                  new_unittest("get_ser", test_get_ser), &
                  new_unittest("compute_force", test_compute_force) &
    ]

  end subroutine collect_verlet3

  subroutine setup(fx)
    type(test_fixture_type), intent(out) :: fx

  !  for logging
    call global_logger%configure(indent=.true., max_width=100)
    call global_logger%configure(level = NONE_LEVEL)

  ! Initialization
    fx%tol = 0.0001_wp
    fx%tau = 0.2_wp                                     ! tau
    fx%delta = 0.5_wp                                          ! delta
    !delta = 1                                          ! delta
    ! atom A
    fx%mass(1) = 1.0080_wp
    fx%x(1, :) = (/-10.0, 0.0, 0.0/)          ! x, y, z
    fx%v(1, :) = (/1.0, 0.0, 0.0/)        ! vx, vy, vz
    ! atom B
    fx%mass(2) = 1.0080_wp                             ! m, x, y, z, vx, vy, vz
    fx%x(2, :) = (/10.0, 0.0, 0.0/)           ! m, x, y, z, vx, vy, vz
    fx%v(2, :) = (/0.0, 0.0, 0.0/)        ! m, x, y, z, vx, vy, vz
    ! atom C
    fx%mass(3) = 1.0080_wp                           ! m, x, y, z, vx, vy, vz
    fx%x(3, :) = (/11.40065, 0.0, 0.0/)           ! m, x, y, z, vx, vy, vz
    !x(3, :) = (/11.0, 0.0, 0.0/)           ! m, x, y, z, vx, vy, vz
    fx%v(3, :) = (/0.0, 0.0, 0.0/)        ! m, x, y, z, vx, vy, vz

    fx%expected_eudist = 20.0_wp

    fx%expected_eudist_with_delta_12_1 = 19.5_wp
    fx%expected_eudist_with_delta_13_2 = 21.40649_wp
    fx%expected_eudist_with_delta_23_3 = 1.487219_wp

    fx%expected_get_ser = (/20.0_wp, 21.40065_wp, 1.40065_wp/)
    !fx%expected_get_ser_delta
    fx%expected_forces(1, :) = (/7.5123e-007, -6.7676E-009, -6.7677E-009/)
    fx%expected_forces(2, :) = (/5.026129, 0.07112365, 0.07112365/)
    fx%expected_forces(3, :) = (/1.505611, 0.07112365, 0.07112365/)

  end subroutine setup


  subroutine test_eudist_p1_p2(error)
    type(error_type), allocatable, intent(out) :: error
    type(test_fixture_type) :: fx

    call setup(fx)

    fx%real_result = eudist(fx%x(1, :), fx%x(2, :))

    write(fx%log_msg, '(A, F5.1, A, F5.1)') "test_eudist_p1_p2: ", fx%real_result, " =? ", fx%expected_eudist
    CALL global_logger%log_warning(fx%log_msg)

    call check(error, fx%expected_eudist, fx%real_result, thr=fx%tol)
  end subroutine test_eudist_p1_p2

  subroutine test_eudist_p2_p1(error)
    type(error_type), allocatable, intent(out) :: error
    type(test_fixture_type) :: fx

    call setup(fx)

    fx%real_result = eudist(fx%x(2, :), fx%x(1, :))
    call check(error, fx%expected_eudist, fx%real_result, thr=fx%tol)
  end subroutine test_eudist_p2_p1

  subroutine test_eudist_with_delta_12_1(error)
    type(error_type), allocatable, intent(out) :: error
    type(test_fixture_type) :: fx

    call setup(fx)

    fx%real_result = eudist_with_delta(fx%x(1, :), fx%x(2, :), fx%delta, 1)

    write(fx%log_msg, '(A, F5.1, A, F5.1)') "test_eudist_with_delta_12_1: ", fx%real_result, " =? ", fx%expected_eudist_with_delta_12_1
    CALL global_logger%log_warning(fx%log_msg)

    call check(error, fx%expected_eudist_with_delta_12_1, fx%real_result, thr=fx%tol)
  end subroutine test_eudist_with_delta_12_1

  subroutine test_eudist_with_delta_13_2(error)
    type(error_type), allocatable, intent(out) :: error
    type(test_fixture_type) :: fx

    call setup(fx)

    fx%real_result = eudist_with_delta(fx%x(1, :), fx%x(3, :), fx%delta, 2)
    call check(error, fx%expected_eudist_with_delta_13_2, fx%real_result, thr=fx%tol)
  end subroutine test_eudist_with_delta_13_2

  subroutine test_eudist_with_delta_23_3(error)
    type(error_type), allocatable, intent(out) :: error
    type(test_fixture_type) :: fx

    call setup(fx)

    fx%real_result = eudist_with_delta(fx%x(2, :), fx%x(3, :), fx%delta, 3)
    call check(error, fx%expected_eudist_with_delta_23_3, fx%real_result, thr=fx%tol)
  end subroutine test_eudist_with_delta_23_3

  subroutine test_get_ser(error)
    type(error_type), allocatable, intent(out) :: error
    type(test_fixture_type) :: fx
    integer :: i

    call setup(fx)

    call get_ser(fx%x(1, :), fx%x(2, :), fx%x(3, :), fx%ser)

    write(fx%log_msg, '(A, 3F5.1, A, 3F5.1)') "test_get_ser: ", fx%ser, " =? ", fx%expected_get_ser
    CALL global_logger%log_warning(fx%log_msg)

    do i = 1, size(fx%ser)
        call check(error, fx%expected_get_ser(i), fx%ser(i), thr=fx%tol)
    end do
  end subroutine test_get_ser

  subroutine test_compute_force(error)
    type(error_type), allocatable, intent(out) :: error
    type(test_fixture_type) :: fx
    integer :: i, j

    call setup(fx)

    call compute_force(fx%x, fx%delta, fx%forces)

    do i = 1, size(fx%forces, 1)
        write(fx%log_msg, '(A, I0, A, 3F5.1, A, 3F5.1)') "test_compute_force: ", i, " - ", fx%forces(i, :), " =? ", fx%expected_forces(i, :)
        CALL global_logger%log_warning(fx%log_msg)
        do j = 1, size(fx%forces(i, :), 1)
            call check(error, fx%expected_forces(i, j), fx%forces(i, j), thr=fx%tol)
        end do
    end do
  end subroutine test_compute_force

!test_compute_force
!        @assertEqual((/7.5123e-007, -6.7676E-009, -6.7677E-009/), forces(1, :), tolerance=tol)
!        @assertEqual((/123.386101, -6.26126, -6.26126/), forces(2, :), tolerance=tol)
!        @assertEqual((/-11.178783, -6.26126, -6.26126/), forces(3, :), tolerance=tol)
!setup
end module test_verlet3

program tester
  use, intrinsic :: iso_fortran_env, only : error_unit
  use testdrive, only : run_testsuite
  use test_verlet3
  implicit none
  integer :: stat

  stat = 0
  call run_testsuite(collect_verlet3, error_unit, stat)

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop
  end if

end program tester
