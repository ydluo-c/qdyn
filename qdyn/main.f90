!Quasi-Dynamic
! for 3D problem, kernel instability occurs while L/dw > 100
! MvdE: this is a terrible location for a warning like that...
! TODO: create a warning entry into log if L/dw > 100

program main

  use problem_class
  use input
  use initialize
  use solver
  use my_mpi, only: init_mpi
  use derivs_all
  use unittests, only: init_tests, kernel_export
  use logger, only : init_log
  use constants, only : DEBUG, RESTART, VERBOSE, FID_LOG

  type(problem_type), pointer :: pb
  integer :: i
  character(len=32) :: arg
  allocate(pb)

  odepb => pb

  ! Loop over potential commandline arguments
  do i = 1, iargc()
    ! Get next argument
    call getarg(i, arg)

    ! Toggle debug
    if (arg == "debug") DEBUG = .true.

    ! Toggle restart
    if (arg == "restart") RESTART = .true.

    ! Toggle verbosity
    if (arg == "verbose") VERBOSE = .true.

    ! Unit tests
    if (arg == "test") then
      DEBUG = .true.
      RESTART = .false.
      VERBOSE = .true.
      call init_mpi()
      call init_log()
      call init_tests(pb)
      stop
    endif

    ! Export kernel
    if (arg == "kernel") then
      call kernel_export(pb)
      stop
    endif

  enddo

  ! If debugging: set verbose to false
  if (VERBOSE .and. DEBUG) VERBOSE = .false.
  ! If (still) verbose: change to stdout
  if (VERBOSE) FID_LOG = 6
  call init_mpi()
  call init_log()
  call read_main(pb)
  call init_all(pb)
  call solve(pb)

end program main
