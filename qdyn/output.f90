module output

! OUTPUT: This module manages outputs

  use logger, only : log_msg, log_debug
  use constants, only : DEBUG

  implicit none
  private

  character(255) :: msg

  public :: ot_read_stations, initialize_output, write_output, log_write_header

contains

!===============================================================================
! Initialise the output modules

subroutine initialize_output(pb)

  use problem_class
  use my_mpi, only: is_MPI_master, is_MPI_parallel
  use constants, only : PI
  type (problem_type) :: pb
  integer :: nbase, nobj

  ! The spirit of the output module is as follows:
  ! First, a container of pointers is created that is shared by all the output
  ! submodules (time series, snapshots, etc). Then, for each output type, the
  ! indices pointing to the appropriate quantities are stored in a list. When
  ! a particular output is requested, the corresponding submodule calls for the
  ! appropiate quantity stored in the list.

  ! Number of objects in containers
  ! Base number
  nbase = pb%nobj
  ! Total number
  nobj = nbase
  ! If thermal pressurisation is requested: add 2 more objects to total
  if (pb%features%tp == 1) then
    nobj = nobj + 2
  endif
  
  
  ! If fluid diffusion is requested: add 1 more objects to total
  if (pb%features%fluid_diff == 1) then
    nobj = nobj + 1
  endif

  ! Overwrite number of objects to output
  pb%nobj = nobj
  
    
  ! Allocate containers
  allocate(pb%objects_glob(nobj))
  allocate(pb%objects_loc(nobj))

  ! If parallel: initialise (allocate) global quantities
  if (is_MPI_parallel()) then
    call init_pb_global(pb)
  ! If serial: point global quantities to local ones
  else
    ! Mesh parameters
    pb%mesh%xglob => pb%mesh%x
    pb%mesh%yglob => pb%mesh%y
    pb%mesh%zglob => pb%mesh%z
    pb%mesh%fault_label_glob => pb%mesh%fault_label
    ! Mechanical quantities
    pb%v_glob => pb%v
    pb%theta_glob => pb%theta
    pb%tau_glob => pb%tau
    pb%dtau_dt_glob => pb%dtau_dt
    pb%slip_glob => pb%slip
    pb%sigma_glob => pb%sigma
    ! If thermal pressurisation is requested: P and T
    if (pb%features%tp == 1) then
      pb%P_glob => pb%P
      pb%T_glob => pb%T
    endif
     ! If fluid diffusion is requested: P
    if (pb%features%fluid_diff == 1) then
      pb%P_glob => pb%P
    endif
    ! Max rupture stats
    pb%tau_max_glob => pb%tau_max
    pb%t_rup_glob => pb%t_rup
    pb%v_max_glob => pb%v_max
    pb%t_vmax_glob => pb%t_vmax
  endif

  ! Assign global quantities (for output)

  ! [integer scalar] step, index of max slip rate
  pb%objects_glob(1)%i => pb%it
  pb%objects_glob(2)%i => pb%ivmax
  ! [double scalar] time, potency, potency rate
  pb%objects_glob(1)%s => pb%time
  pb%objects_glob(2)%s => pb%pot
  pb%objects_glob(3)%s => pb%pot_rate
  ! [double vector] spatial coordinates
  pb%objects_glob(1)%v => pb%mesh%xglob
  pb%objects_glob(2)%v => pb%mesh%yglob
  pb%objects_glob(3)%v => pb%mesh%zglob
  ! [double vector] mechanical quantities
  pb%objects_glob(4)%v => pb%v_glob
  pb%objects_glob(5)%v => pb%theta_glob
  pb%objects_glob(6)%v => pb%tau_glob
  pb%objects_glob(7)%v => pb%dtau_dt_glob
  pb%objects_glob(8)%v => pb%slip_glob
  pb%objects_glob(9)%v => pb%sigma_glob
  pb%objects_glob(10)%v => pb%pot_fault
  pb%objects_glob(11)%v => pb%pot_rate_fault
  ! CRP: [double vector] quantity related to vmax of each fault
  pb%objects_glob(12)%v => pb%vmax_fault
  ! [double vector] fault label
  pb%objects_glob(1)%vi => pb%mesh%fault_label_glob
  ! CRP: [double vector] quantity related to ivmax of each fault
  pb%objects_glob(2)%vi => pb%ivmax_fault

  ! If thermal pressurisation is requested, add P and T
  if (pb%features%tp == 1) then
    ! [double vector] pressure, temperature
    pb%objects_glob(nbase+1)%v => pb%P_glob
    pb%objects_glob(nbase+2)%v => pb%T_glob
  endif
   ! If fluid diffusion is requested, add P 
  if (pb%features%fluid_diff == 1) then
    ! [double vector] pressure, temperature
    pb%objects_glob(nbase+1)%v => pb%P_glob
  endif

  ! Assign local quantities (which need to be synchronised)
  ! Scalars do not need to be synchronised, so can be skipped

  ! [double vector] spatial coordinates
  pb%objects_loc(1)%v => pb%mesh%x
  pb%objects_loc(2)%v => pb%mesh%y
  pb%objects_loc(3)%v => pb%mesh%z
  ! [double vector] mechanical quantities
  pb%objects_loc(4)%v => pb%v
  pb%objects_loc(5)%v => pb%theta
  pb%objects_loc(6)%v => pb%tau
  pb%objects_loc(7)%v => pb%dtau_dt
  pb%objects_loc(8)%v => pb%slip
  pb%objects_loc(9)%v => pb%sigma

  ! If thermal pressurisation is requested, add P and T
  if (pb%features%tp == 1) then
    ! [double vector] pressure, temperature
    pb%objects_loc(nbase+1)%v => pb%P
    pb%objects_loc(nbase+2)%v => pb%T
  endif
   ! If fluid diffusion is requested, add P 
  if (pb%features%fluid_diff == 1) then
    ! [double vector] pressure, temperature
    pb%objects_loc(nbase+1)%v => pb%P
  endif

  ! Init ot, ox, screen
  call ot_init(pb)
  call ox_init(pb)

end subroutine initialize_output

!===============================================================================
! Trigger output modules

subroutine write_output(pb)

  use problem_class
  type (problem_type) :: pb
  logical :: last_call

  last_call = (pb%it == pb%itstop)

  if (DEBUG) then
    write(msg, *) "log_write"
    call log_debug(msg, pb%it)
  endif

  ! Call a print to screen if requested
  if (mod(pb%it, pb%ntout_log) == 0 .or. last_call) then
    call log_write(pb)
  endif

  if (DEBUG) then
    write(msg, *) "ot_write"
    call log_debug(msg, pb%it)
  endif

  ! Write time series if requested
  if (mod(pb%it, pb%ot%ntout) == 0 .or. last_call) then
    call ot_write(pb)
  endif

  if (DEBUG) then
    write(msg, *) "ox_write"
    call log_debug(msg, pb%it)
  endif

  ! ox_write has its own internal checks for which output to write
  call ox_write(pb)

  if (DEBUG) then
    write(msg, *) "end write_output"
    call log_debug(msg, pb%it)
  endif

end subroutine write_output

!=====================================================================
! Write log header
subroutine log_write_header(pb)

  use problem_class
  use constants, only : PI
  use my_mpi, only : is_MPI_master

  type (problem_type), intent(inout) :: pb
  double precision :: K
  character(255) :: msg

  if (.not. is_MPI_master()) return

  ! SEISMIC: skip calculating critical stiffness for CNS model
  ! Is not very useful
  if (pb%i_rns_law == 3) return

  ! TODO: what's this?
  if (pb%ot%ic < 1) return

  call log_msg("Values at selected point of the fault:")

  K = pb%mesh%Lfault
  if (pb%mesh%dim == 1) then
    if (.not. pb%kernel%k2f%finite) K = pb%mesh%W
  endif

  K = PI*pb%smu/K
  if (pb%mesh%dim < 2) then
    write(msg, *) "K/Kc = ", K/(pb%sigma(pb%ot%ic)*(pb%b(pb%ot%ic)-pb%a(pb%ot%ic))/pb%dc(pb%ot%ic))
    call log_msg(msg)
    write(msg, *) "K/Kb = ", K/(pb%sigma(pb%ot%ic)*pb%b(pb%ot%ic)/pb%dc(pb%ot%ic))
    call log_msg(msg)
  end if

  call log_msg("     it,   dt (sec), time (yrs), v_max (m/s), sigma_max (MPa)")

end subroutine log_write_header


!=====================================================================
! Write to log
subroutine log_write(pb)

  use constants, only : YEAR
  use problem_class
  use my_mpi, only : is_MPI_parallel, is_mpi_master, max_allproc

  type (problem_type), intent(in) :: pb
  double precision :: sigma_max, sigma_max_glob
  character(255) :: msg

  ! Synchronise max sigma across nodes
  sigma_max = maxval(pb%sigma)
  if (is_MPI_parallel()) then
    call max_allproc(sigma_max,sigma_max_glob)
  else
    sigma_max_glob = sigma_max
  endif

  ! Log message
  write(msg, "(i7,4(x,e11.3))")  &
    pb%it, pb%dt_did, pb%time/YEAR, pb%vmaxglob, sigma_max_glob/1.0D6
  call log_msg(msg)

end subroutine log_write

!=====================================================================
! CRP: Is this subrutine used somewhere?
subroutine ot_read_stations(ot)

  use problem_class, only : ot_type
  use constants, only: FID_STATIONS

  type (ot_type), intent(inout) :: ot

  integer :: nsta,ista

  open(unit=FID_STATIONS, file='stations.dat', action='read')
  read(FID_STATIONS, *) nsta
  allocate(ot%xsta(nsta),ot%ysta(nsta),ot%zsta(nsta))
  do ista=1,nsta
    read(FID_STATIONS, *) ot%xsta(ista), ot%ysta(ista), ot%zsta(ista)
  enddo
  close(FID_STATIONS)

end subroutine ot_read_stations

!=====================================================================
! Write ot file header
subroutine ot_init(pb)

  use problem_class
  use constants, only:  BIN_OUTPUT, FID_IASP, FID_OT, FID_VMAX, FID_FAULT, &
                        FILE_IASP, FILE_OT, FILE_VMAX, FILE_FAULT, RESTART
  use my_mpi, only: is_MPI_parallel, is_MPI_master, gather_allvi_root
  use mesh, only: mesh_get_size, nnLocal_perproc, nnoffset_glob_perproc

  type (problem_type), intent(inout) :: pb
  integer :: i, id, iasp_count, iot_count, n, niasp, niot, nnGlobal, Ncols
  integer, dimension(pb%mesh%nnglob) :: iasp_buf, iot_buf
  integer, allocatable, dimension(:) :: iasp_list, iot_list

  character(len=100) :: tmp
  character(len=100), allocatable :: iot_name

  ! Number of mesh elements
  n = mesh_get_size(pb%mesh)
  nnGlobal = pb%mesh%nnglob

  ! Allocate slip rate history for IASP output
  allocate(pb%ot%v_pre(nnGlobal), pb%ot%v_pre2(nnGlobal))
  pb%ot%v_pre = 0.d0
  pb%ot%v_pre2 = 0.d0

  ! Number of ot output quantities
  pb%ot%not = 12
  ! Number of ot_vmax output quantities
  pb%ot%not_vmax = 10
  ! If thermal pressurisation is requested, add 2 more
  if (pb%features%tp == 1) then
    pb%ot%not = pb%ot%not + 2
    pb%ot%not_vmax = pb%ot%not_vmax + 2
  endif
    ! If fluid diffusion is requested, add 1 more
  if (pb%features%fluid_diff == 1) then
    pb%ot%not = pb%ot%not + 1
    pb%ot%not_vmax = pb%ot%not_vmax + 1
  endif
  ! Allocate space in array of pointers
  allocate(pb%ot%fmt(pb%ot%not))
  allocate(pb%ot%fmt_vmax(pb%ot%not_vmax))

  ! Default output format
  pb%ot%fmt = "(e15.7)"
  ! Simulation step is an integer
  pb%ot%fmt(1) = "(i15)"
  ! Time needs higher precision
  pb%ot%fmt(2) = "(e24.14)"
  ! Fault label is an integer
  pb%ot%fmt(pb%ot%not) = "(i15)"

  ! Default vmax output format
  pb%ot%fmt_vmax = "(e15.7)"
  ! Simulation step is an integer
  pb%ot%fmt_vmax(1) = "(i15)"
  ! Time needs higher precision
  pb%ot%fmt_vmax(2) = "(e24.14)"
  ! vmax location is an integer
  pb%ot%fmt_vmax(3) = "(i15)"
  ! Fault label is an integer
  pb%ot%fmt_vmax(pb%ot%not_vmax) = "(i15)"

  ! If parallel:
  if (is_MPI_parallel()) then
    ! Combine local OT indices
    call gather_allvi_root( pb%ot%iot, n, iot_buf, nnLocal_perproc, &
                            nnoffset_glob_perproc, nnGlobal)
    ! Combine local IASP indices
    call gather_allvi_root( pb%ot%iasp, n, iasp_buf, nnLocal_perproc, &
                            nnoffset_glob_perproc, nnGlobal)
  ! If serial:
  else
    ! Point global to local indices
    iot_buf = pb%ot%iot
    iasp_buf = pb%ot%iasp
  endif

  if (is_MPI_master()) then

    ! Get the number of OT output indices
    niot = sum(iot_buf)
    ! Allocate space for a list of iot indices
    allocate(iot_list(niot))
    iot_count = 0
    ! Loop over all potential OT locations
    do i=1,nnGlobal
      ! If location is OT location
      if (iot_buf(i) == 1) then
        iot_count = iot_count + 1
        iot_list(iot_count) = i
      endif
    enddo

    ! Overwrite IOT in pb%ot
    pb%ot%iot = iot_list

    ! Get the number of IASP output indices
    niasp = sum(iasp_buf)
    ! Allocate space for a list of iasp indices
    allocate(iasp_list(niasp))
    iasp_count = 0
    ! Loop over all potential OT locations
    do i=1,nnGlobal
      ! If location is OT location
      if (iasp_buf(i) == 1) then
        iasp_count = iasp_count + 1
        iasp_list(iasp_count) = i
      endif
    enddo

    ! Overwrite IASP in pb%ot
    pb%ot%iasp = iasp_list

  endif
  
  ! Open files, write headers (if not restarting the simulation with time of last simulation)
  if (is_MPI_master()) then

    ! Time series output: loop over all OT locations
    do i=1,niot
      id = FID_OT + pb%ot%iot(i) - 1
      ! Write headers
      write(tmp, "(a, i0)") FILE_OT, pb%ot%iot(i) - 1
      iot_name = trim(tmp)
      if (.not. RESTART) then
        open(id, file=iot_name, status="replace")
        write(id, "(a)") "# macroscopic values:"
        write(id, "(a)") "# 1=step, 2=t, 3=pot, 4=pot_rate"
        write(id, "(a)") "# values at selected point:"
        write(id, "(a)") "# 5=V, 6=theta, 7=tau, 8=dtau_dt, 9=slip, 10=sigma"
        if (pb%features%tp == 1) then
          write(id, "(a)") "# 12=P, 13=T"
        endif
         if (pb%features%fluid_diff == 1) then
          write(id, "(a)") "# 12=P"
        endif
        write(id, "(a)") "# last=fault_label"
        close(id)
      endif 
    enddo

    ! Vmax output
    if (.not. RESTART) then
      open(FID_VMAX, file=FILE_VMAX, status="replace")
      write(FID_VMAX, "(a)") "# values at max(V) location:"
      write(FID_VMAX, "(a)") "# 1=step, 2=t, 3=ivmax, 4=v, 5=theta, 6=tau, 7=dtau_dt, 8=slip, 9=sigma"
      if (pb%features%tp == 1) then
        write(FID_VMAX, "(a)") "# 11=P, 12=T"
      endif
      if (pb%features%fluid_diff == 1) then
        write(FID_VMAX, "(a)") "# 11=P"
      endif
      write(FID_VMAX, "(a)") "# last=fault_label"
      close(FID_VMAX)
    endif

    ! IASP output
    if (.not. RESTART) then
      open(FID_IASP, file=FILE_IASP, status="replace")
      write(FID_IASP, "(a)") "# Seismicity record:"
      write(FID_IASP, "(a)") "# 1=i, 2=t, 3=v, 4=fault_label"
      close(FID_IASP)
    endif

    ! Fault output (potency, potency rate)
    if (.not. RESTART) then
      open(FID_FAULT, file=FILE_FAULT, status="replace")
      write(FID_FAULT, "(a)") "# fault values:"
      write(FID_FAULT, "(a)", advance="no") "# 1=step, 2=t "

      ! number of column for output
      Ncols = 5

      do i=1, pb%nfault
        write(FID_FAULT, "(a, i0, a, i0)", advance="no") " ", Ncols, "=pot_", i
        Ncols = Ncols+1
        write(FID_FAULT, "(a, i0, a, i0)", advance="no") " ", Ncols, "=pot_rate_", i
        Ncols = Ncols+1
        !CRP: write headers for ivmax_fault and vmax_fault
        write(FID_FAULT, "(a, i0, a, i0)", advance="no") " ", Ncols, "=vmax_", i
        Ncols = Ncols+1
        write(FID_FAULT, "(a, i0, a, i0)", advance="no") " ", Ncols, "=ivmax_", i
        Ncols = Ncols+1
      enddo

      close(FID_FAULT)

    endif

  endif

end subroutine ot_init


!=====================================================================
! write ox file header
subroutine ox_init(pb)

  use problem_class
  use constants, only: FID_OX, FILE_OX, FID_OX_LAST, FILE_OX_LAST, RESTART
  use my_mpi, only : is_MPI_parallel, is_mpi_master

  type (problem_type), intent(inout) :: pb

  integer :: count_w, count_x, n

  ! Number of mesh elements
  n = pb%mesh%nn
  ! Number of ox elements
  count_w = ceiling(pb%mesh%nw / float(pb%ox%nwout))
  count_x = ceiling(pb%mesh%nx / float(pb%ox%nxout))
  pb%ox%count = count_w * count_x
  ! Initial potency
  pb%ox%pot_pre = 0.d0
  ! Initial snapshot count
  pb%ox%seq_count = 0

  ! ---------------------------------------------------------------------------
  ! Prepare data structure and headers for ox, ox_dyn, and dyn output

  ! Determine the number of objects
  pb%ox%nox = 12
  pb%ox%nrup = 2
  if (pb%features%tp == 1) then
    pb%ox%nox = pb%ox%nox + 2
  endif
    if (pb%features%fluid_diff == 1) then
    pb%ox%nox = pb%ox%nox + 1
  endif
  ! Allocate space for dynamic output that needs to be carried over
  allocate(pb%ox%objects_rup(pb%ox%nrup))
  ! Allocate space for the output format
  allocate(pb%ox%fmt(pb%ox%nox))
  ! Default output format
  pb%ox%fmt = "(e20.12)"
  ! Simulation step is an integer
  pb%ox%fmt(1) = "(i15)"
  ! Time needs higher precision
  pb%ox%fmt(2) = "(e24.14)"
  ! Theta may require even higher precision
  pb%ox%fmt(7) = "(e20.12)"
  ! tau and sigma
  pb%ox%fmt(8) = "(e20.12)"
  pb%ox%fmt(11) = "(e20.12)"
  ! Fault label is an integer
  pb%ox%fmt(pb%ox%nox) = "(i15)" 

  ! Case 1: serial execution. All quantities are local
  if (.not. is_MPI_parallel()) then
    ! Assign local tau/v to rupture object
    pb%ox%objects_rup(1)%v => pb%tau
    pb%ox%objects_rup(2)%v => pb%v

  ! Case 2: parallel execution with output master. All quantities are global
  elseif (is_MPI_parallel() .and. is_MPI_master()) then
    ! Assign global tau/v to rupture object
    pb%ox%objects_rup(1)%v => pb%tau_glob
    pb%ox%objects_rup(2)%v => pb%v_glob

  endif

  ! Create ox file (if not split over multiple files)
  if (pb%ox%i_ox_seq == 0) then
    if (.not. RESTART) then
      open(FID_OX, file=FILE_OX, status="replace")
      close(FID_OX)
    endif
  endif

  ! Create ox file for last snapshot 
  if (pb%ox%i_ox_seq == 0) then
    open(FID_OX_LAST, file=FILE_OX_LAST, status="replace")
    close(FID_OX_LAST)
  endif

  ! Define headers
  pb%ox%header = '# step t x y z v theta tau tau_dot slip sigma fault_label'
  if (pb%features%tp == 1) then
    pb%ox%header = '# step t x y z v theta tau tau_dot slip sigma P T fault_label'
  endif
    if (pb%features%fluid_diff == 1) then
    pb%ox%header = '# step t x y z v theta tau tau_dot slip sigma P fault_label'
  endif

  ! If ox_dyn is requested
  if (pb%ox%i_ox_dyn == 1) then
    ! Allocate rupture quantities
    allocate ( pb%t_rup(n), pb%tau_max(n), pb%v_max(n), pb%t_vmax(n) )
    pb%t_rup = 0d0
    pb%tau_max = 0d0
    pb%v_max = 0d0
    pb%t_vmax = 0d0

    ! If parallel, the global quantities are already allocated during
    ! init_pb_global. If serial, point these global quantities to local ones.
    if (.not. is_MPI_parallel()) then
      pb%tau_max_glob => pb%tau_max
      pb%t_rup_glob => pb%t_rup
      pb%v_max_glob => pb%v_max
      pb%t_vmax_glob => pb%t_vmax
    endif
  endif

  ! If dynamic output is requested (ox_dyn or QSB)
  if ((pb%ox%i_ox_dyn == 1) .or. (pb%DYN_FLAG == 1)) then
    ! Number of dynamic events (start at 0)
    pb%ox%dyn_count = 0
    ! Dynamic event status (0 = no, 1 = yes)
    pb%ox%dyn_stat = 0
  endif

  ! End preparing data structures/headers
  ! ---------------------------------------------------------------------------

end subroutine ox_init

!=====================================================================
! Export timeseries
subroutine ot_write(pb)

  use problem_class
  use constants, only:  BIN_OUTPUT, FID_IASP, FID_OT, FID_VMAX, FID_FAULT, &
                        FILE_IASP, FILE_OT, FILE_VMAX, FILE_FAULT
  use my_mpi, only: is_MPI_master, is_MPI_parallel

  type (problem_type), intent(inout) :: pb
  integer :: i, id, iot, istart, ivmax, k, n, niasp, niot

  character(len=100) :: tmp
  character(len=100), allocatable :: iot_name


  ! If parallel: do sync
  if (is_MPI_parallel()) then
    if (DEBUG) then
      write(msg, *) "pb_global"
      call log_debug(msg, pb%it)
    endif
    call pb_global(pb)
  endif

  if (DEBUG) then
    write(msg, *) "get_ivmax"
    call log_debug(msg, pb%it)
  endif

  ! Get the maximum slip rate
  call get_ivmax(pb)
  ivmax = pb%ivmax

  if (DEBUG) then
    write(msg, *) "calc_potency"
    call log_debug(msg, pb%it)
  endif

  ! TODO: merge potency calculations
  !! Calculate potency (rate)
  call calc_potency(pb)

  ! Get the maximum slip rate per fault and its corresponding index
  call get_ivmax_fault(pb)
  ! ivmax_fault = pb%ivmax_fault_glob
  ! vmax_fault = pb%vmax_fault_glob


  ! If this is master proc: write output
  if (is_MPI_master()) then

    ! Number of OT locations
    niot = size(pb%ot%iot)

    if (DEBUG) then
      write(msg, *) "write ot"
      call log_debug(msg, pb%it)
    endif

    ! Loop over OT locations
    do iot=1,niot
      id = FID_OT + pb%ot%iot(iot) - 1
      ! Format file name
      write(tmp, "(a, i0)") FILE_OT, pb%ot%iot(iot) - 1
      iot_name = trim(tmp)
      ! Open OT file
      open(id, file=iot_name, status="old", access="append")
      ! Write data (one line)
      call write_ot_lines(id, pb%ot%fmt, pb%objects_glob, pb%ot%iot(iot), pb)
      ! Close file
      close(id)
    enddo

    ! Skip the first 3 elements of the pointer container (x/y/z)
    istart = 4
    ! Size of the container
    k = pb%ot%not_vmax

    if (DEBUG) then
      write(msg, *) "write vmax"
      call log_debug(msg, pb%it)
    endif

    ! Write vmax data
    open(FID_VMAX, file=FILE_VMAX, status="old", access="append")
    ! Write step
    write(FID_VMAX, pb%ot%fmt_vmax(1), advance="no") pb%objects_glob(1)%i
    ! Write time
    write(FID_VMAX, pb%ot%fmt_vmax(2), advance="no") pb%objects_glob(1)%s
    ! Write index of vmax location
    write(FID_VMAX, pb%ot%fmt_vmax(3), advance="no") pb%objects_glob(2)%i
    ! Write vector data (at vmax location), except last one
    do i=istart,k-1
      ! Write quantity, do not advance to next line
      write(FID_VMAX, pb%ot%fmt_vmax(i), advance="no") pb%objects_glob(i)%v(ivmax)
    enddo
    ! Write last ot output quantity, advance to next line
    write(FID_VMAX, pb%ot%fmt_vmax(k)) pb%objects_glob(1)%vi(ivmax)
    ! Close file
    close(FID_VMAX)

    if (DEBUG) then
      write(msg, *) "write IASP"
      call log_debug(msg, pb%it)
    endif

    ! IASP output slip velocity maxima at selected nodes
    niasp = size(pb%ot%iasp)
    ! Open file
    open(FID_IASP, file=FILE_IASP, status="old", access="append")
    ! Loop over all IASP locations
    do i=1,niasp
      ! Mesh index of IASP location
      n = pb%ot%iasp(i)
      ! Write critera:
      !   - v1 > v_th
      !   - v2 < v1
      !   - v1 > v0
      if ((pb%ot%v_pre(n) >= pb%ot%v_th) .and. &
          (pb%v_glob(n) < pb%ot%v_pre(n)) .and. &
          (pb%ot%v_pre(n) >= pb%ot%v_pre2(n))) then
          write(FID_IASP, "(i10, 2e24.16)") n, pb%time, pb%ot%v_pre(n), pb%mesh%fault_label(n)
      endif
    enddo
    close(FID_IASP)
    ! Update v0, v1
    pb%ot%v_pre2 = pb%ot%v_pre
    pb%ot%v_pre = pb%v_glob

    if (DEBUG) then
      write(msg, *) "write fault data"
      call log_debug(msg, pb%it)
    endif

    ! Write fault data (potency, potency rate and slip_dt)
    open(FID_FAULT, file=FILE_FAULT, status="old", position="append")
    call write_fault_lines(FID_FAULT, pb%ot%fmt, pb%objects_glob, pb)
    close(FID_FAULT)
  endif

end subroutine ot_write

!=====================================================================
! Export snapshots
! Do not confuse with subrutine write_ox (writing snapshots)

subroutine ox_write(pb)

  use problem_class
  use constants, only:  FID_MW, FID_OX, FID_OX_LAST, FID_OX_DYN, FID_QSB_POST, FID_QSB_PRE, &
                        FILE_OX, FILE_OX_LAST, FILE_OX_DYN_MAX, FILE_OX_DYN_POST, FILE_OX_DYN_PRE
  use my_mpi, only: is_MPI_parallel, is_mpi_master

  type (problem_type), intent(inout) :: pb

  integer :: iw, ix, n, nwout_dyn, nxout_dyn, unit
  logical ::  call_gather, close_unit, dynamic, falling_edge, last_call, &
              MPI_master, skip, rising_edge, write_ox, write_ox_last, write_ox_dyn, &
              write_ox_seq, write_QSB
  double precision, dimension(pb%mesh%nnglob) :: tau, v
  character(len=100) :: tmp
  character(len=100), allocatable :: file_name

  nwout_dyn = pb%ox%nwout_dyn
  nxout_dyn = pb%ox%nxout_dyn

  !---------------------------------------------------------------------------
  ! Perform series of checks to see which operations need to be executed

  ! Is this the last call (last simulation step)?
  last_call = (pb%it == pb%itstop)

  ! Is this proc MPI master? (default .true. for serial)
  MPI_master = is_MPI_master()

  ! Check if we're crossing dynamic thresholds from below (rising edge)
  ! or from above (falling edge)
  rising_edge = ((pb%ox%dyn_stat == 0) .and. &
                (pb%v_glob(pb%ivmax) >= pb%DYN_th_on))
  falling_edge =  ((pb%ox%dyn_stat == 1) .and. &
                  (pb%v_glob(pb%ivmax) <= pb%DYN_th_off))
  ! If we're crossing a dynamic threshold, update global dynamic state
  if (rising_edge .and. MPI_master) then
    ! Set global dynamic state to 1
    pb%ox%dyn_stat = 1
  endif
  ! We're dynamic as long as dyn_stat == 1
  dynamic = (pb%ox%dyn_stat == 1)
  skip = (pb%ox%dyn_count <= pb%dyn_skip)

  ! Should this proc write ox, ox_dyn, or QSB output?
  write_ox = (mod(pb%it, pb%ox%ntout) == 0 .or. last_call)
  !write_ox_last = last_call
  ! CRP: trying to debug empty write_ox_last
  write_ox_last = (mod(pb%it, pb%ox%ntout) == 0 .or. last_call)
  write_ox_dyn =  ((pb%ox%i_ox_dyn == 1) .and. dynamic) .and. .not. skip
  write_ox_seq = write_ox .and. (pb%ox%i_ox_seq == 1)
  write_QSB = ((pb%DYN_FLAG == 1) .and. dynamic) .and. .not. skip

  ! Call an MPI gather when:
  !  1. Output is requested (either regular ox or other flavour)
  !  2. AND parallel execution
  call_gather = (write_ox .or. write_ox_last .or. write_ox_dyn .or. write_QSB) .and. &
                is_MPI_parallel()

  ! Update write output only if this proc is MPI master
  write_ox = write_ox .and. MPI_master
  write_ox_last = write_ox_last .and. MPI_master
  write_ox_dyn = write_ox_dyn .and. MPI_master
  write_ox_seq = write_ox_seq .and. MPI_master
  write_QSB = write_QSB .and. MPI_master

  ! Does the output unit need to be closed?
  close_unit = last_call

  if (DEBUG) then
    write(msg, *) "pb_global"
    call log_debug(msg, pb%it)
  endif

  ! End checks
  !---------------------------------------------------------------------------
  ! Perform an MPI global gather

  if (call_gather) then
    call pb_global(pb)
  endif

  if (DEBUG) then
    write(msg, *) "get_ivmax"
    call log_debug(msg, pb%it)
  endif

  ! Get the maximum slip rate
  call get_ivmax(pb)

  ! End MPI global gather
  !---------------------------------------------------------------------------
  ! Write regular ox output

  if (write_ox) then

    pb%ox%seq_count = pb%ox%seq_count + 1

    ! Check if data needs to be written to separate file
    if (write_ox_seq) then
      write(tmp, "(a, a, i0)") FILE_OX, "_", pb%ox%seq_count
      file_name = trim(tmp)
      open(FID_OX, file=file_name, status="replace")
    else
      open(FID_OX, file=FILE_OX, status="old", position="append")
    endif

    if (DEBUG) then
      write(msg, *) "write_ox_lines (ox)"
      call log_debug(msg, pb%it)
    endif

    ! Write ox data
    call write_ox_lines(FID_OX, pb%ox%fmt, pb%objects_glob, &
                        pb%ox%nxout, pb%ox%nwout, pb)
    close(FID_OX)

  endif

  ! End regular ox output
  !---------------------------------------------------------------------------
  ! Write ox_last output
  
  ! Snapshot of the last time-step of the simulation. It is useful to have a snapshot with the 
  ! full mesh resolution in case of needing to re-start the simulation not from
  ! the beggining. This avoids the need to generate an output of all the snapshots
  ! with the full resolution
  
  if (write_ox_last) then
    open(FID_OX_LAST, file=FILE_OX_LAST, status="replace")

    if (DEBUG) then
      write(msg, *) "write_ox_lines (last)"
      call log_debug(msg, pb%it)
    endif
    
    ! Write ox data (in full resolution)
    call write_ox_lines(FID_OX_LAST, pb%ox%fmt, pb%objects_glob, 1, 1, pb)
    close(FID_OX_LAST)
  
  endif

  !---------------------------------------------------------------------------
  ! Write ox_dyn output

  if (write_ox_dyn) then

    ! If we cross the dynamic threshold from below
    if (rising_edge) then

      ! Collect output quantities
      do iw=1, pb%mesh%nw, nwout_dyn
        do ix=1, pb%mesh%nx, nxout_dyn
          n = (iw - 1) * pb%mesh%nx + ix
          pb%tau_max_glob(n) = pb%tau_glob(n)
          pb%v_max_glob(n) = pb%v_glob(n)
        enddo
      enddo
      pb%t_rup_glob = pb%time
      pb%t_vmax_glob = pb%time

      ! Define unit for rising edge ox output
      unit = FID_OX_DYN + 1
      ! Write ox data
      write(tmp, "(a, i0)") FILE_OX_DYN_PRE, pb%ox%dyn_count
      file_name = trim(tmp)
      open(unit, file=file_name, status="replace")

      if (DEBUG) then
        write(msg, *) "write_ox_lines (dyn; rising)"
        call log_debug(msg, pb%it)
      endif

      call write_ox_lines(unit, pb%ox%fmt, pb%objects_glob, &
                          pb%ox%nxout_dyn, pb%ox%nwout_dyn, pb)
      ! Close output unit
      close(unit)

    ! If we cross the dynamic threshold from above
    elseif (falling_edge) then

      ! Skip the first dyn_skip events
      if (pb%ox%dyn_count > pb%dyn_skip) then
        ! First: falling edge ox output
        ! Define unit for falling edge ox output
        unit = FID_OX_DYN + 2
        write(tmp, "(a, i0)") FILE_OX_DYN_POST, pb%ox%dyn_count
        file_name = trim(tmp)
        ! Write ox data
        open(unit, file=file_name, status="replace")

        if (DEBUG) then
          write(msg, *) "write_ox_lines (dyn; falling)"
          call log_debug(msg, pb%it)
        endif

        call write_ox_lines(unit, pb%ox%fmt, pb%objects_glob, &
                            pb%ox%nxout_dyn, pb%ox%nwout_dyn, pb)
        ! Close output unit
        close(unit)

        ! Second: max rupture stats
        ! Define unit for max stats output
        unit = FID_OX_DYN + 3
        write(tmp, "(a, i0)") FILE_OX_DYN_MAX, pb%ox%dyn_count
        file_name = trim(tmp)
        open(unit, file=file_name, status="replace")
        write(unit,'(a)') '# x y z t_rup tau_max t_vmax vmax'
        ! WARNING: nx should be nxglob, but this has not been implemented yet
        ! In general, QDYN will break when FFT_TYPE == 2
        do iw=1, pb%mesh%nw, nwout_dyn
          do ix=1, pb%mesh%nx, nxout_dyn
            n = (iw - 1) * pb%mesh%nx + ix
            write(unit, '(3e15.7,4e28.20)') &
              pb%mesh%xglob(n), pb%mesh%yglob(n), pb%mesh%zglob(n), &
              pb%t_rup_glob(n), pb%tau_max_glob(n), pb%t_vmax_glob(n), &
              pb%v_max_glob(n), pb%mesh%fault_label_glob(n)
          enddo
        enddo
        ! Close output unit
        close(unit)
      endif

    ! No threshold crossing, but still dynamic
    else

      ! WARNING: %v is a placeholder for any vector, so do not replace
      ! %v with %tau (i.e. the below is not a bug...)
      tau = pb%ox%objects_rup(1)%v
      v = pb%ox%objects_rup(2)%v
      ! Update max stress/slip rate for each ox element
      do iw=1, pb%mesh%nw, nwout_dyn
        do ix=1, pb%mesh%nx, nxout_dyn
          n = (iw - 1) * pb%mesh%nx + ix
          ! Mark location of rupture front (max tau at given location)
          if (tau(n) > pb%tau_max_glob(n)) then
            pb%tau_max_glob(n) = tau(n)
            pb%t_rup_glob(n) = pb%time
          endif
          ! Mark location of max slip rate (at given location)
          if (v(n) > pb%v_max_glob(n)) then
            pb%v_max_glob(n) = v(n)
            pb%t_vmax_glob(n) = pb%time
          endif
        enddo
      enddo
    endif

  endif

  ! End ox_dyn output
  !---------------------------------------------------------------------------
  ! Write QSB output
  ! Note that for QSB output, the data for each element needs to be written to
  ! the output file, so that nxout = nwout = 1

  if (write_QSB) then

    if (rising_edge) then

      unit = FID_QSB_PRE
      open(unit, file='DYN_PRE.txt', status='REPLACE')
      call write_ox_lines(unit, pb%ox%fmt, pb%objects_glob, 1, 1, pb)
      pb%ox%pot_pre = pb%pot
      close(unit)

    else if (falling_edge) then

      unit = FID_QSB_PRE
      open(unit, file='DYN_POST.txt', status='REPLACE')
      call write_ox_lines(unit, pb%ox%fmt, pb%objects_glob, 1, 1, pb)
      close(unit)

    endif

  endif

  ! End QSB output
  !---------------------------------------------------------------------------
  ! Check dynamic state

  ! TODO: this needs to be implemented somewhere (not in QSB!!!)
  ! ! Check for seismic moment
  ! if ((pb%pot - pb%ox%pot_pre) * pb%smu >= pb%DYN_M) then
  !   ! TODO: replace dyn_count with QSB count !!
  !   pb%ox%dyn_count = pb%ox%dyn_count + 1
  !   if (pb%ox%dyn_count > pb%DYN_SKIP) then
  !     pb%itstop = pb%it
  !     write(FID_MW, '(3i10,e24.14)') pb%it, pb%ivmax, pb%ox%dyn_count, pb%time
  !   endif
  ! endif

  if (falling_edge .and. MPI_master) then
    ! Set global dynamic state back to 0
    pb%ox%dyn_stat = 0
    ! Increment event counter
    pb%ox%dyn_count = pb%ox%dyn_count + 1
  endif

  ! End check dynamic state
  !---------------------------------------------------------------------------

  ! All done

  if (DEBUG) then
    write(msg, *) "end write_ox"
    call log_debug(msg, pb%it)
  endif

end subroutine ox_write

!=====================================================================
! Write ot data to file
subroutine write_ot_lines(unit, fmt, objects, iot, pb)

  use problem_class
  type (problem_type), intent(inout) :: pb
  integer :: i, iot, istart, k, unit
  character(len=16), dimension(pb%nobj) :: fmt
  type(optr), dimension(pb%nobj) :: objects

  ! Skip the first 3 elements of the pointer container (x/y/z)
  istart = 4
  ! Size of the container
  k = pb%ot%not

  ! Write step
  write(unit, fmt(1), advance="no") objects(1)%i
  ! Write time
  write(unit, fmt(2), advance="no") objects(1)%s
  ! Write pot, pot_rate
  write(unit, fmt(3), advance="no") objects(2)%s
  write(unit, fmt(4), advance="no") objects(3)%s
  ! Loop over all ot output quantities
  do i=istart,k-3
    ! Write ot output quantity, do not advance to next line
    write(unit, fmt(i+1), advance="no") objects(i)%v(iot)
  enddo
  ! Write fault ID, advance to next line
  write(unit, fmt(k)) objects(1)%vi(iot)

end subroutine write_ot_lines

!=====================================================================
! Write fault data (potency, potency rate and delta slip) to file
subroutine write_fault_lines(unit, fmt, objects, pb)

  use problem_class
  type (problem_type), intent(inout) :: pb
  integer :: i, unit
  character(len=16), dimension(pb%nobj) :: fmt
  type(optr), dimension(pb%nobj) :: objects


  ! Write step
  write(unit, fmt(1), advance="no") objects(1)%i
  ! Write time
  write(unit, fmt(2), advance="no") objects(1)%s

  do i=1, pb%nfault
    ! CRP: Write pot_fault, pot_rate_fault, ivmax_fault and vmax_fault (do not advance to next line)
    write(unit, fmt(3), advance="no") objects(10)%v(i) !pot_fault
    write(unit, fmt(4), advance="no") objects(11)%v(i) !pot_rate_fault
    write(unit, "(e15.7)", advance="no") objects(12)%v(i) !vmax_fault
    ! If we're at the last fault, write vmax_fault and advance to next line
    if (i == pb%nfault) then
      write(unit, "(i15)") objects(2)%vi(i) !ivmax_fault
    ! Else do not advance
    else
      write(unit, "(i15)", advance="no") objects(2)%vi(i)
    endif

  enddo

end subroutine write_fault_lines

!=====================================================================
! Write ox data to file
subroutine write_ox_lines(unit, fmt, objects, nxout, nwout, pb)

  use problem_class
  type (problem_type), intent(inout) :: pb
  integer :: iox, iwout, ixout, k, n, nw, nwout, nx, nxout, unit
  character(len=16), dimension(pb%nobj) :: fmt
  type(optr), dimension(pb%nobj) :: objects

  k = pb%ox%nox

  ! Number of nodes to loop over (either local or global)
  nx = pb%mesh%nxglob
  nw = pb%mesh%nwglob

  ! Write header
  write(unit,'(a)') trim(pb%ox%header)
  ! Loop over all ox elements
  do iwout=1, nw, nwout
    do ixout=1, nx, nxout
      n = (iwout - 1) * nx + ixout
      ! Write step
      write(unit, fmt(1), advance="no") objects(1)%i
      ! Write time
      write(unit, fmt(2), advance="no") objects(1)%s
      ! Loop over all ox output quantities
      do iox=1,9
        ! Write ox output quantity, do not advance to next line
        write(unit, fmt(iox+2), advance="no") objects(iox)%v(n)
      enddo
      ! Add P/T if needed 
      ! There is a bug here? should be iox+3 and iox+4? (13 and 14 knowing that nbase is 12)
      ! iox should be incremented by one for each column feature
      ! fmt has length = 14, but the last element is the fault number 
      if (pb%features%tp == 1) then
        write(unit, fmt(iox+2+1), advance="no") objects(iox+1)%v(n)
        write(unit, fmt(iox+2+2), advance="no") objects(iox+2)%v(n)
      endif
      
        ! Add P/T if needed
      if (pb%features%fluid_diff == 1) then
        !print*,'iox',iox
!         print*,'pb%ox%nox',pb%ox%nox
       !  print*,'iox',iox
!         print*,'pb%ot%not',pb%ot%not
!         print*,'pb%ox%fmt',pb%ox%fmt
        write(unit, fmt(iox+2), advance="no") objects(iox+3)%v(n)
        ! print*,'n',n
!         print*,'objects(iox+1)%v',objects(iox+1)%v(n)

      endif
      ! Write fault number, advance to next line
      write(unit, fmt(k)) objects(1)%vi(n)
    enddo
  enddo
  
 !  print*,'nxout',nxout
!   print*,'nwout',nwout
!   print*,'objects(iox+1)%v',objects(11)%v
  
end subroutine write_ox_lines

!=====================================================================

! Compute the potency (rate) on the fault, defined as:
! 0d:  [pot/pot_rate] = [slip/v] * L
! 1d:  [pot/pot_rate] = sum([slip/v] * dx)
! 2d:  [pot/pot_rate] = sum([slip/v] * dx * dw)

subroutine calc_potency(pb)

  use problem_class
  use my_mpi, only: is_MPI_parallel, sum_allreduce

  type(problem_type), intent(inout) :: pb
  double precision, dimension(pb%mesh%nn) :: area
  double precision, dimension(pb%nfault) :: pot, pot_rate
  double precision, dimension(1) :: buf
  integer :: iw, ix, n, lbl

  pot = 0d0
  pot_rate = 0d0
  area = 0d0

  ! Step 1: define area vector
  do iw=1, pb%mesh%nw
    do ix=1, pb%mesh%nx
      n = (iw - 1) * pb%mesh%nx + ix
      area(n) = pb%mesh%dx(ix) * pb%mesh%dw(iw)
    end do
  end do

  ! Step 2: multiply [slip/v] * area and sum
  ! This is done for each fault separately
  do n=1, pb%mesh%nn
    lbl = pb%mesh%fault_label(n)
    pot(lbl) = pot(lbl) + pb%slip(n) * area(n)
    pot_rate(lbl) = pot_rate(lbl) + pb%v(n) * area(n)
  end do

  ! Step 3: synchronise all faults across all nodes
  do lbl=1, pb%nfault
    if (is_MPI_parallel()) then
      buf(1) = pot(lbl)
      call sum_allreduce(buf(1), 1)
      pot(lbl) = buf(1)

      buf(1) = pot_rate(lbl)
      call sum_allreduce(buf(1), 1)
      pot_rate(lbl) = buf(1)

    endif
  enddo

  pb%pot = sum(pot)
  pb%pot_fault = pot
  pb%pot_rate = sum(pot_rate)
  pb%pot_rate_fault = pot_rate

  ! Done

end subroutine calc_potency

!=====================================================================
! Get the location of the maximum slip rate

subroutine get_ivmax(pb)

  use problem_class

  type(problem_type), intent(inout) :: pb

  pb%ivmax = maxloc(pb%v_glob, 1)

end subroutine get_ivmax

!=====================================================================

! CRP: Get the location of the maximum slip rate in each fault

subroutine get_ivmax_fault(pb)

  use problem_class
  ! use my_mpi, only: is_MPI_parallel, sum_allreduce

  type(problem_type), intent(inout) :: pb

  integer, dimension(pb%nfault) :: ivmax_fault
  double precision, dimension(pb%nfault) :: vmax_fault
  ! double precision, dimension(1) :: buf
  integer :: n, lbl, mesh_size, i_max_v
  double precision :: max_v

  ! Initialize ivmax_fault and vmax_fault to 0 for all faults
  ivmax_fault = 0
  vmax_fault = 0d0 
  
  ! mesh size
  mesh_size=size(pb%v_glob)

  do lbl =1, pb%nfault
    max_v = 0d0
    i_max_v = 0

    do n=1, mesh_size
      if (pb%mesh%fault_label_glob(n) == lbl) then
        if (pb%v_glob(n) > max_v) then
          max_v = pb%v_glob(n)
          i_max_v = n 
        end if
      end if
    end do

    vmax_fault(lbl) = max_v
    ivmax_fault(lbl) = i_max_v

  end do


  if (DEBUG) then
    write(msg, *) "ivmax_fault=", ivmax_fault
    call log_debug(msg, pb%it)
    write(msg, *) "vmax_fault=", vmax_fault
    call log_debug(msg, pb%it)
    write(msg, *) "size pb%ivmax_fault=", size(pb%ivmax_fault)
    call log_debug(msg, pb%it)
    write(msg, *) "size pb%vmax_fault=", size(pb%vmax_fault)
    call log_debug(msg, pb%it)
  endif

  ! Store the final vmax_fault array in pb%vmax_fault
  pb%vmax_fault = vmax_fault

  if (DEBUG) then
    write(msg, *) " after storing pb%vmax_fault=", pb%vmax_fault
    call log_debug(msg, pb%it)
  endif


    ! Store the final ivmax_fault array in pb%ivmax_fault
  pb%ivmax_fault = ivmax_fault

  if (DEBUG) then
    write(msg, *) " after storing pb%ivmax_fault=", pb%ivmax_fault
    call log_debug(msg, pb%it)
  endif

end subroutine get_ivmax_fault


!=====================================================================
! Initiate MPI global gather

subroutine init_pb_global(pb)

  use problem_class

  type(problem_type), intent(inout) :: pb
  integer :: n

  n = pb%mesh%nnglob

  ! Check if global variables are already allocated (should not be!)
  if (.not. pb%allocated_glob) then

    ! Base quantities
    allocate( pb%v_glob(n), pb%theta_glob(n), pb%tau_glob(n), &
              pb%dtau_dt_glob(n), pb%slip_glob(n), pb%sigma_glob(n))

    pb%v_glob = 0d0
    pb%theta_glob = 0d0
    pb%tau_glob = 0d0
    pb%dtau_dt_glob = 0d0
    pb%slip_glob = 0d0
    pb%sigma_glob = 0d0

    ! If thermal pressurisation is requested, allocate P and T
    if (pb%features%tp == 1) then
      allocate(pb%P_glob(n), pb%T_glob(n))
      pb%P_glob = 0d0
      pb%T_glob = 0d0
    endif
    
    ! If fluid diffusion is requested, allocate P
    if (pb%features%fluid_diff == 1) then
      allocate(pb%P_glob(n))
      pb%P_glob = 0d0
    endif

    ! Allocate rupture max stats
    allocate( pb%tau_max_glob(n), pb%t_rup_glob(n), &
              pb%v_max_glob(n), pb%t_vmax_glob(n))

    pb%tau_max_glob = 0d0
    pb%t_rup_glob = 0d0
    pb%v_max_glob = 0d0
    pb%t_vmax_glob = 0d0

    ! Update flag
    pb%allocated_glob = .true.

  ! If global quantities are already allocated, something is wrong!
  else
    call log_msg("Error in output.f90::init_pb_global")
    call log_msg("Global quantities are allocated while they shouldn't be")
    stop "Terminating..."
  endif

end subroutine init_pb_global

!=====================================================================
! Collect global fault nodes to master processor for outputs

subroutine pb_global(pb)

  use problem_class
  use mesh, only: nnLocal_perproc, nnoffset_glob_perproc
  use my_mpi, only: gather_allvdouble_root

  type(problem_type), intent(inout) :: pb
  integer :: i, k, nnLocal, nnGlobal

  k = pb%nobj
  nnLocal = pb%mesh%nn
  nnGlobal = pb%mesh%nnglob

  ! Loop over proc-specific quantities
  ! Skip x, y, z
  do i=4,k

    ! Skip fault potency (rate) and vmax_fault
    ! TODO: replace this hardcoded stuff with something
    ! more elegant...
    if (i == 10 .or. i == 11 .or. i == 12) cycle

    ! Initialise to zero
    pb%objects_glob(i)%v = 0d0

    ! Call MPI gather
    call gather_allvdouble_root(  pb%objects_loc(i)%v, nnLocal, &
                                  pb%objects_glob(i)%v, nnLocal_perproc, &
                                  nnoffset_glob_perproc, nnGlobal)
  enddo

end subroutine pb_global
!=====================================================================


end module output
