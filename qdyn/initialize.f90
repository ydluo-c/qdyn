module initialize

  implicit none
  private

  public :: init_all

contains

!=============================================================
subroutine init_all(pb)

  use problem_class
  use mesh, only : init_mesh, mesh_get_size, read_mesh_nodes
  use constants, only : PI, SOLVER_TYPE
  use my_mpi, only: is_MPI_master
  use logger, only : log_msg
  use fault_stress, only : init_kernel
  use friction, only : set_theta_star, compute_velocity_RSF
  use friction_cns, only : compute_velocity, dphi_dt
  use solver, only : init_rk45
  use diffusion_solver, only: init_tp
  use ode_rk45_2, only : init_rk45_2
  use output, only: initialize_output, log_write_header
!!$  use omp_lib

  type(problem_type), intent(inout) :: pb

  integer :: n
  character(255) :: msg
!  integer :: TID, NTHREADS

  ! Allocate scalar (pointer) quantities
  allocate(pb%ivmax)
  allocate(pb%pot, pb%pot_rate)
  pb%ivmax = 0
  pb%pot = 0d0
  pb%pot_rate = 0d0

  call init_mesh(pb%mesh)

 ! number of equations
 ! Initial number of neqs is defined in problem_class.f90
  pb%neqs = pb%neqs + pb%features%localisation
  if (pb%features%stress_coupling == 1 .and. pb%mesh%dim == 2) then
    pb%neqs = pb%neqs + 1
  endif
  
  ! If permeability change is asked for
  if (pb%features%var_k == 1) then 
    pb%neqs = pb%neqs + 1
  endif

 ! dt_max & perturbation
 ! if periodic loading, set time step smaller than a fraction of loading period
  if (pb%Aper /= 0.d0 .and. pb%Tper > 0.d0) then
    if (pb%dt_max > 0.d0) then
      pb%dt_max = min(pb%dt_max,0.2d0*pb%Tper)
    else
      pb%dt_max = 0.2d0*pb%Tper
    endif
  endif
  if (pb%Tper > 0.d0) then
    pb%Omper = 2.d0*PI/pb%Tper
  else
    pb%Omper = 0.d0
  endif

 ! impedance
  if (pb%beta > 0d0) then
    pb%zimpedance = 0.5d0*pb%smu/pb%beta
  else
    pb%zimpedance = 0.d0
  endif

  ! Log impedance
  write(msg, "(a, e15.3)") "Impedance = ", pb%zimpedance
  call log_msg(msg)
  !---------------------- impedance ------------------

  !---------------------- ref_value ------------------
  n = mesh_get_size(pb%mesh)
  ! SEISMIC: initialise pb%tau in input.f90 to be compatible with CNS model
  allocate ( pb%dtau_dt(n), pb%theta_star(n) )
  pb%dtau_dt = 0d0
  call set_theta_star(pb)

  ! SEISMIC: initialise thermal pressurisation model (diffusion_solver.f90)
  allocate(pb%P(pb%mesh%nn))
  pb%P = 0d0
  
! If fluid diffusion 
  if (pb%features%fluid_diff == 1) then
    pb%P = pb%fluid_diff%P_a
    pb%fluid_diff%P_temp = pb%P
  endif
  
  ! If variable permeability
  if (pb%features%var_k == 1) then
    pb%var_k%kstar = pb%var_k%kmin + (pb%fluid_diff%permeability-pb%var_k%kmin) * exp(abs((pb%sigma-pb%P))/abs(pb%var_k%Snk))
  endif

  
  
  
! If thermal pressurisation  
  if (pb%features%tp == 1) then
    allocate(pb%T(pb%mesh%nn))
    pb%T = 0d0
    call init_tp(pb)
    call log_msg("Spectral mesh initiated")
  endif

  ! SEISMIC: the CNS model has the initial shear stress defined in the
  ! input file, so we can skip the initial computation of friction
  if (pb%i_rns_law /= 3) then
    pb%v = compute_velocity_RSF(pb%tau, pb%sigma-pb%P, pb%theta, pb)
  endif
  if (pb%i_rns_law == 3) then
    pb%v = compute_velocity(pb%tau, pb%sigma-pb%P, pb%theta, pb%theta2, pb)
    if (pb%features%tp == 1) then
      call dphi_dt( pb%v, pb%tau, pb%sigma-pb%P, pb%theta, pb%theta2, &
                    pb%dtheta_dt, pb%dtheta2_dt, pb)
    endif
  endif

  call init_kernel( pb%lam, pb%smu, pb%mesh, pb%kernel, pb%D, pb%H, &
                    pb%features%stress_coupling, pb%finite, pb%test%test_mode)

  call initialize_output(pb)

  ! SEISMIC: initialise Runge-Kutta ODE solver, if selected
  if (SOLVER_TYPE == 2) then
    call init_rk45(pb)
  endif
  if (SOLVER_TYPE == 3) then
    call init_rk45_2(pb)
  endif

  call log_msg("Initialization completed")
  call log_write_header(pb)

  ! Info about threads
!!$OMP PARALLEL PRIVATE(NTHREADS, TID)
!!$  TID = OMP_GET_THREAD_NUM()
!!$  write(6,*) 'Thread index = ', TID
!!$OMP BARRIER
!!$  if (TID == 0) then
!!$    NTHREADS = OMP_GET_NUM_THREADS()
!!$    write(6,*) 'Total number of threads = ', NTHREADS
!!$  end if
!!$OMP END PARALLEL

end subroutine init_all

end module initialize
