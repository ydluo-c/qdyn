!qdyn input all parameters

module input

  implicit none
  private

  public :: read_main

contains
!=====================================================================
! Read input file qdyn.in
! If MPI, there is one qdyn*.in file per processor
!
subroutine read_main(pb)

  use problem_class
  use mesh, only : read_mesh_parameters, read_mesh_nodes, mesh_get_size
  use output, only : ot_read_stations
  use logger, only : log_msg
  use my_mpi, only : my_mpi_tag, is_MPI_parallel
  use constants, only : FAULT_TYPE, SOLVER_TYPE, FID_IN, RESTART

  type(problem_type), intent(inout) :: pb

  double precision, dimension(:), allocatable :: read_buf
  integer :: i, j, n, N_cols, N_creep

  call log_msg("Start reading input...")

  open(unit=FID_IN, file='qdyn'//trim(my_mpi_tag())//'.in', action='read')

  call read_mesh_parameters(FID_IN, pb%mesh)
  call log_msg("   Mesh input complete")

  if (pb%mesh%dim==1) read(FID_IN, *) pb%finite
  read(FID_IN, *) pb%itheta_law
  read(FID_IN, *) pb%i_rns_law

  ! Read restart time and step
  allocate(pb%time, pb%it)
  read(FID_IN, *) pb%time, pb%it

  ! Read number of faults
  read(FID_IN, *) pb%nfault
  ! Allocate potency (rate): 1 for each fault
  allocate(pb%pot_fault(pb%nfault), pb%pot_rate_fault(pb%nfault))
  pb%pot_fault = 0d0
  pb%pot_rate_fault = 0d0

  ! CRP: allocate ivmax_fault and vmax_fault
  allocate(pb%ivmax_fault(pb%nfault), pb%vmax_fault(pb%nfault))
  pb%ivmax_fault = 0
  pb%vmax_fault = 0d0

  ! SEISMIC: various simulation features can be turned on (1) or off (0)
  if (pb%i_rns_law == 3) then
    read(FID_IN, *) pb%cns_params%N_creep
  endif
  read(FID_IN, *) pb%features%stress_coupling, pb%features%tp, pb%features%fluid_diff,pb%features%var_k , pb%features%localisation
  read(FID_IN, *) pb%ntout_log, pb%ot%ntout, pb%ox%ntout, pb%ot%ic, pb%ox%nxout, pb%ox%nwout, &
             pb%ox%nxout_dyn, pb%ox%nwout_dyn, pb%ox%i_ox_seq, pb%ox%i_ox_dyn
  read(FID_IN, *) pb%beta, pb%smu, pb%lam, pb%D, pb%H, pb%ot%v_th
  read(FID_IN, *) pb%Tper, pb%Aper
  read(FID_IN, *) pb%dt_try, pb%dt_max,pb%tmax, pb%acc
  read(FID_IN, *) pb%NSTOP
  read(FID_IN, *) pb%DYN_FLAG, pb%DYN_SKIP
  read(FID_IN, *) pb%DYN_M, pb%DYN_th_on, pb%DYN_th_off
  read(FID_IN, *) FAULT_TYPE, SOLVER_TYPE
  
  ! In case of fluid diffusion, read the number of source
  if (pb%features%fluid_diff == 1) then
      read(FID_IN, *) pb%fluid_diff%nb_source
      ! If there are some source, allocate and read the associated parameters
      if (pb%fluid_diff%nb_source>=1) then
          ! Allocate the arrays
          allocate( pb%fluid_diff%t_injection_beg(pb%fluid_diff%nb_source),  &
                    pb%fluid_diff%t_injection_end(pb%fluid_diff%nb_source),  &
                    pb%fluid_diff%index_injection(pb%fluid_diff%nb_source),  &
                    pb%fluid_diff%Q(pb%fluid_diff%nb_source))
          
          ! Read the parameters
          read(FID_IN, *) pb%fluid_diff%Q
          read(FID_IN, *) pb%fluid_diff%index_injection
          read(FID_IN, *) pb%fluid_diff%t_injection_beg
          read(FID_IN, *) pb%fluid_diff%t_injection_end
          
      endif
  endif
  
  call log_msg("  Flags input complete")

  n = mesh_get_size(pb%mesh) ! number of nodes in this processor
  allocate ( pb%tau(n), pb%sigma(n), pb%v(n), pb%theta(n),  &
             pb%v_star(n), pb%ot%iot(n), pb%ot%iasp(n), &
             pb%dc(n), pb%coh(n), pb%v_pl(n) )

  !JPA if MPI, read only the nodes of this processor
  ! <SEISMIC>
  ! Read input parameters for the CNS model. These parameters are (in order):
  !   sigma:            applied normal stress (effective if no TP)
  !   tau:              initial shear stress
  !   theta:            initial porosity (phi in CNS formulation)
  !   v_pl:             load-point velocity
  !   a_tilde:          coefficient of logarithmic rate dependence
  !   mu_tilde_star:    reference friction coefficient at y_gr_star
  !   y_gr_star:        reference granular fow strain rate
  !   H:                dilatancy geometric factor
  !   phi_c:            critical state porosity
  !   phi0:             lower cut-off porosity
  !
  !   (for each out of N creep mechanisms)
  !   A:                creep law kinetic constant
  !   n:                stress exponent
  !   m:                porosity exponent
  !
  !   L:                total thickness of the fault zone
  !   iot:              flag for other time series output
  !   iasp:             flag for identification purposes
  !
  ! Note that these parameters are material (gouge) properties, and are
  ! generally not spatically uniform, and hence are allocatable
  ! See friction.f90 for a description of and references to the CNS model
  ! See user manual for detailed definitions of the above parameters
  ! (this is a permanent TODO)
  ! </SEISMIC>

  ! If the CNS model is selected
  if (pb%i_rns_law == 3) then
    N_creep = pb%cns_params%N_creep
    ! Define the number of columns to read: 13 plus 3 for each creep mechanism
    N_cols = 13 + 3*N_creep
    ! Allocate read buffer
    allocate( read_buf(N_cols) )
    ! Allocate CNS parameters
    allocate( pb%cns_params%a_tilde(n), pb%cns_params%mu_tilde_star(n), &
              pb%cns_params%y_gr_star(n), pb%cns_params%H(n), &
              pb%cns_params%phi_c(n), pb%cns_params%phi0(n), &
              pb%cns_params%L(n) )
    ! Allocate CNS parameters for N creep mechanisms
    allocate( pb%cns_params%A(N_creep*n), pb%cns_params%n(N_creep*n), &
              pb%cns_params%m(N_creep*n) )
    do i=1,n

      ! Read data into buffer
      read(FID_IN, *) read_buf(:)

      ! General fault parameters
      pb%sigma(i) = read_buf(1)
      pb%tau(i) = read_buf(2)
      pb%theta(i) = read_buf(3)
      pb%v_pl(i) = read_buf(4)

      ! Granular flow parameters
      pb%cns_params%a_tilde(i) = read_buf(5)
      pb%cns_params%mu_tilde_star(i) = read_buf(6)
      pb%cns_params%y_gr_star(i) = read_buf(7)
      pb%cns_params%H(i) = read_buf(8)

      ! Porosity parameters
      pb%cns_params%phi_c(i) = read_buf(9)
      pb%cns_params%phi0(i) = read_buf(10)

      ! Creep mechanism parameters (3 per mechanism)
      do j=1,N_creep
        ! Kinetic constant
        pb%cns_params%A((i-1)*N_creep+j) = read_buf(11+(j-1)*3)
        ! Stress exponent
        pb%cns_params%n((i-1)*N_creep+j) = read_buf(12+(j-1)*3)
        ! Porosity exponent
        pb%cns_params%m((i-1)*N_creep+j) = read_buf(13+(j-1)*3)
      enddo

      ! Fault zone thickness
      pb%cns_params%L(i) = read_buf(N_cols-2)

      ! Output parameters
      pb%ot%iot(i) = int(read_buf(N_cols-1))
      pb%ot%iasp(i) = int(read_buf(N_cols))

      ! Set dummy values to unused parameters
      pb%dc(i) = 1.0
      pb%V(i) = 0.0

    end do
    ! Deallocate the read buffer
    deallocate(read_buf)
  ! Else, the RSF model is selected
  else
    allocate ( pb%a(n), pb%inv_a(n), pb%b(n), pb%v1(n), &
               pb%v2(n), pb%mu_star(n), pb%inv_visc(n) )
    do i=1,n
      read(15,*)pb%sigma(i), pb%tau(i), pb%theta(i),  &
                pb%a(i), pb%b(i), pb%dc(i), pb%v1(i), &
                pb%v2(i), pb%mu_star(i), pb%v_star(i), &
                pb%ot%iot(i), pb%ot%iasp(i), pb%coh(i), pb%v_pl(i), &
                pb%inv_visc(i)
    end do

    ! Precompute 1/a to avoid many divisions later on
    pb%inv_a = 1.0 / pb%a
  endif
    
  ! <SEISMIC>
  ! Read input parameters for the localisation model (CNS only).
  ! These parameters and corresponding units are (in order):
  !
  if (pb%features%localisation == 1) then
    ! Raise an error if the CNS model is not selected
    if (pb%i_rns_law /= 3) then
      call log_msg("Localisation of shear strain is compatible only with the CNS model (i_rns_law = 3)")
      stop
    endif

    N_creep = pb%cns_params%N_creep
    ! Define the number of columns to read: 2 plus 3 for each creep mechanism
    N_cols = 2 + 3*N_creep
    ! Allocate read buffer
    allocate( read_buf(N_cols) )
    allocate( pb%cns_params%lambda(n), pb%theta2(n) )
    allocate( pb%cns_params%A_bulk(N_creep*n), pb%cns_params%n_bulk(N_creep*n),&
              pb%cns_params%m_bulk(N_creep*n) )

    do i=1,n

      ! Read data into buffer
      read(FID_IN, *) read_buf(:)

      ! Degree of localisation (0 <= lambda <= 1)
      pb%cns_params%lambda(i) = read_buf(1)
      ! Porosity in bulk gouge
      pb%theta2(i) = read_buf(2)

      ! Creep mechanism parameters (3 per mechanism)
      do j=1,N_creep
        ! Kinetic constant
        pb%cns_params%A_bulk((i-1)*N_creep+j) = read_buf(3+(j-1)*3)
        ! Stress exponent
        pb%cns_params%n_bulk((i-1)*N_creep+j) = read_buf(4+(j-1)*3)
        ! Porosity exponent
        pb%cns_params%m_bulk((i-1)*N_creep+j) = read_buf(5+(j-1)*3)
      enddo

    end do
    ! Deallocate the read buffer
    deallocate(read_buf)
  ! Else, set lambda and theta2 to dummy values
  else
    allocate(pb%cns_params%lambda(n), pb%theta2(n))
    do i=1,n
      pb%cns_params%lambda(i) = 1d0
      pb%theta2(i) = 0d0
    end do
  endif

  ! End reading localisation model parameters
  ! </SEISMIC>

  ! <SEISMIC>
  ! Read input parameters for the thermal pressurisation (TP) model.
  ! These parameters and corresponding units are (in order):
  !
  !   rhoc:   density times specific heat capacity of host rock [J/K/m^3]
  !   beta:   bulk compressibility [1/Pa]
  !   eta:    dynamic viscocity [Pa s]
  !   w:      half-width of shear distribution profile [m]
  !   k_t:    thermal conductivity [J/s/K/m]
  !   k_p:    intrinsic hydraulic permeability [m^2]
  !   l:      nett thermal expansion coefficient [1/K]
  !   P_a:    ambient fluid pressure [Pa]
  !   T_a:    ambient temperature [K]
  !
  if (pb%features%tp == 1) then
    allocate (  pb%tp%rhoc(n), pb%tp%beta(n), pb%tp%eta(n), pb%tp%w(n), &
                pb%tp%k_t(n), pb%tp%k_p(n), pb%tp%l(n), &
                pb%tp%P_a(n), pb%tp%T_a(n), pb%tp%dilat_factor(n) )
    do i=1,n
      read(FID_IN, *) pb%tp%rhoc(i), pb%tp%beta(i), pb%tp%eta(i), pb%tp%w(i), &
                      pb%tp%k_t(i), pb%tp%k_p(i), pb%tp%l(i), &
                      pb%tp%P_a(i), pb%tp%T_a(i), pb%tp%dilat_factor(i)
    end do
  endif
  
  ! End reading TP model parameters
  ! </SEISMIC>

  !#########################################################################################
  ! Read input parameters for the fluid diffusion model. 
  if (pb%features%fluid_diff == 1) then

    allocate (  pb%fluid_diff%rhof(n), pb%fluid_diff%beta(n), pb%fluid_diff%eta(n), &
                pb%fluid_diff%phi(n), pb%fluid_diff%permeability(n),pb%fluid_diff%P_a(n), &
                pb%fluid_diff%P_temp(n), pb%fluid_diff%P_dot_temp(n))
                
    allocate (pb%fluid_diff%permeability_x(n+1)) ! Permeability at the interface 
    do i=1,n
      read(FID_IN, *) pb%fluid_diff%rhof(i), pb%fluid_diff%beta(i), pb%fluid_diff%eta(i), &
                      pb%fluid_diff%phi(i),  pb%fluid_diff%permeability(i), pb%fluid_diff%P_a(i) 
    end do
  endif
  
  !#########################################################################################
  ! Read input parameters for the variable permeability model. 
  
  if (pb%features%var_k == 1) then

    allocate (  pb%var_k%kmin(n), pb%var_k%kmax(n), pb%var_k%L1(n), &
                pb%var_k%T1(n),pb%var_k%kstar(n),pb%var_k%Snk(n),pb%var_k%dkstar_dt(n))
                
    do i=1,n
      read(FID_IN, *) pb%var_k%kmin(i), pb%var_k%kmax(i), pb%var_k%L1(i), &
                       pb%var_k%T1(i),pb%var_k%Snk(i)
    end do
  endif
  
  !#########################################################################################
 
  ! CRP: Instead, call read_mesh_nodes for all mesh types so the fault label can
  ! be written in the outputs for all fault dimensionalities

  call read_mesh_nodes(FID_IN, pb%mesh)

  ! Overwrite slip if restart with time and slip of last simulation
  allocate(pb%slip(n))
  pb%slip = 0d0
  if(RESTART) then
    pb%slip = pb%mesh%restart_slip
  endif

  close(FID_IN)
  call log_msg("Input complete")

end subroutine read_main

end module input
