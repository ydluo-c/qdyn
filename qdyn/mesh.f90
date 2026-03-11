module mesh

  use logger, only : log_msg

  implicit none
  private

  type mesh_type
    integer :: dim = 0  ! dim = 1, 2 ,3 ~xD
    integer :: nx, nw, nn, nnglob, nwglob, nxglob ! along-strike, along-dip, total grid number
    double precision :: Lfault, W, Z_CORNER ! fault length, width, lower-left corner z (follow Okada's convention)
    double precision, allocatable :: DIP_W(:) ! along-dip grid size and dip (adjustable), nw count
    ! Local mesh coordinates
    double precision, pointer :: x(:) => null(), y(:) => null(), z(:) => null()
    double precision, pointer :: restart_slip(:) => null()
    integer, pointer :: fault_label(:) => null()
    ! Local dip, and along-strike and along-dip size of grid cells (nx*nw count)
    ! TODO: grid sizes should be constant in order for FFT to work!
    double precision, pointer :: dip(:) => null(), dx(:) => null(), dw(:) => null()
    ! Global mesh coordinates (nx*nwglobal count)
    double precision, pointer :: xglob(:) => null(), yglob(:) => null(), zglob(:) => null()
    double precision, pointer :: dipglob(:) => null(), dwglob(:) => null()
    double precision, pointer :: restart_slip_glob(:) => null()
    integer, pointer :: fault_label_glob(:) => null()
  end type mesh_type

  ! SEISMIC: discretisation of spectral domain for diffusion solver
  type spectral_mesh_type
    double precision, dimension(:), allocatable :: lw, F_inv ! dimensionless wavenumbers, inv Fourier kernel
    double precision :: Dlogl=0.8, lw_max=10.0 ! logarithmic grid spacing, max dimensionless wavenumber
    integer :: Nl=20 ! number of mesh elements
  end type spectral_mesh_type

  integer, allocatable, save :: nnLocal_perproc(:),nnoffset_glob_perproc(:)

  public :: mesh_type, read_mesh_parameters, read_mesh_nodes, init_mesh, mesh_get_size
  public :: nnLocal_perproc, nnoffset_glob_perproc
  public :: spectral_mesh_type

contains

!=============================================================
subroutine read_mesh_parameters(iin,m)

  type(mesh_type), intent(inout) :: m
  integer, intent(in) :: iin

  integer :: i

  ! problem dimension (1D, 2D or 3D), mesh type
  read(iin,*) m%dim

  !spring-block or 2d problem
  select case (m%dim)

  case(0,1)
    read(iin,*) m%nn
    read(iin,*) m%Lfault, m%W

  case(2) !3d problem
    read(iin,*) m%nx,m%nw
    m%nn = m%nx * m%nw
    read(iin,*) m%Lfault, m%W , m%Z_CORNER   ! JPA m%W is not used in this case, remove it
    allocate(m%dw(m%nw), m%DIP_W(m%nw))
    do i=1,m%nw
      read(iin,*) m%dw(i), m%DIP_W(i)
    end do

  case default
    call log_msg("mesh dimension should be 0, 1 or 2")
    stop

  end select

  if (m%dim==0) m%nn = 1

end subroutine read_mesh_parameters

!=============================================================
subroutine read_mesh_nodes(iin,m)

  type(mesh_type), intent(inout) :: m
  integer, intent(in) :: iin

  integer :: i

  allocate(m%x(m%nn), m%y(m%nn), m%z(m%nn), m%dip(m%nn), m%fault_label(m%nn), m%restart_slip(m%nn))
  do i=1,m%nn
    read(iin,*) m%x(i),m%y(i),m%z(i),m%dip(i), m%fault_label(i), m%restart_slip(i)
  enddo

end subroutine read_mesh_nodes

!=============================================================

integer function mesh_get_size(m) result(n)
  type(mesh_type), intent(in) :: m
  n = m%nn
end function mesh_get_size

!=============================================================

subroutine init_mesh(m)

  type(mesh_type), intent(inout) :: m

  call log_msg("Initializing mesh...")

  m%nnglob = m%nn

  select case (m%dim)
    case(0); call init_mesh_0D(m)
    case(1); call init_mesh_1D(m)
    case(2); call init_mesh_2D(m)
  end select

end subroutine init_mesh

!--------------------------------------------------------
! Spring-block System
subroutine init_mesh_0D(m)

  type(mesh_type), intent(inout) :: m

  call log_msg("Spring-block System")
  allocate(m%x(m%nn), m%y(m%nn), m%z(m%nn), m%fault_label(m%nn), m%restart_slip(m%nn))
  allocate(m%dx(1), m%dw(1))
  m%nx = 1
  m%nw = 1
  m%nxglob = m%nx
  m%nwglob = m%nw
  m%dx = m%Lfault
  m%dw = 1d0
  m%x = 0d0
  m%y = 0d0
  m%z = 0d0
  m%fault_label = 1d0
  m%restart_slip= 0d0

end subroutine init_mesh_0D

!--------------------------------------------------------
! 1D fault, uniform grid
subroutine init_mesh_1D(m)

  type(mesh_type), intent(inout) :: m

  integer :: i

  call log_msg("1D fault, uniform grid")
  allocate(m%dx(m%nn), m%dw(1))
  m%nx = m%nn
  m%nw = 1
  m%dx = m%Lfault/m%nn
  m%dw = 1d0
  allocate(m%x(m%nn), m%y(m%nn), m%z(m%nn), m%fault_label(m%nn), m%restart_slip(m%nn))
  do i=1,m%nn
    m%x(i) = (i-m%nn*0.5d0-0.5d0)*m%dx(1)
    ! Assuming nn is even (usually a power of 2),
    ! the center of the two middle elements (i=nn/2 and nn/2+1)
    ! are located at x=-dx/2 and x=dx/2, respectively
    m%y(i) = 0d0
    m%z(i) = 0d0
    m%fault_label(i) = 1d0
    m%restart_slip(i) = 0d0
  enddo

  m%xglob => m%x
  m%yglob => m%y
  m%zglob => m%z
  m%nxglob = m%nx
  m%nwglob = m%nw
  m%fault_label_glob => m%fault_label
  m%restart_slip_glob => m%restart_slip

end subroutine init_mesh_1D

!--------------------------------------------------------
  ! 2D fault, uniform grid along-strike
  ! Assumptions:
  !   + the fault trace is parallel to x
  !   + the lower "left" corner of the fault is at (0,0,Z_CORNER)
  !   + z is negative downwards (-depth)
  ! Storage scheme: faster index runs along-strike (x)
subroutine init_mesh_2D(m)

  use constants, only : PI
  use my_mpi, only: is_MPI_parallel, my_mpi_NPROCS, my_mpi_tag, &
                    gather_alli, gather_allvdouble, gather_allvi
  use my_mpi, only: is_MPI_master

  type(mesh_type), intent(inout) :: m
  integer :: iproc, nwLocal, nnLocal, nwGlobal, nnGlobal, NPROCS
  integer, allocatable :: nwLocal_perproc(:), nwoffset_glob_perproc(:)

  call log_msg("2D fault, uniform grid along-strike")

  ! Number of elements along-strike is the same for all procs (serial or parallel)
  m%nxglob = m%nx

  ! MvdE: is m%dx used elsewhere in the code? For the FFT?
  allocate(m%dx(m%nx))
  m%dx = m%Lfault/m%nx

  ! Serial mode: point global parameters to local parameters
  if (.not.is_MPI_parallel()) then
    m%xglob => m%x
    m%yglob => m%y
    m%zglob => m%z
    m%dipglob => m%dip
    m%dwglob => m%dw
    m%nwglob = m%nw
    m%fault_label_glob => m%fault_label
    m%restart_slip_glob => m%restart_slip

  ! Parallel mode: gather mesh node locations for computing the kernel
  else

    NPROCS = my_mpi_NPROCS()

    ! Allocations
    allocate(nwLocal_perproc(0:NPROCS-1))
    allocate(nwoffset_glob_perproc(0:NPROCS-1))
    allocate(nnLocal_perproc(0:NPROCS-1))
    allocate(nnoffset_glob_perproc(0:NPROCS-1))
    
    ! Number of elements along-dip (this proc)
    nwLocal = m%nw
    ! Collecte nw from other nodes
    call gather_alli(nwLocal, nwLocal_perproc)
    
    ! Get offsets for all nodes (in terms of nw)
    do iproc = 0, NPROCS - 1
      nwoffset_glob_perproc(iproc) = sum(nwLocal_perproc(0:iproc)) - nwLocal_perproc(iproc)
    enddo

    ! Global number of elements along-dip
    nwGlobal = sum(nwLocal_perproc)
    m%nwglob = nwGlobal

    ! Local number of nodes
    nnLocal = m%nx * nwLocal
    ! Number of nodes for local to all procs
    nnLocal_perproc = nwLocal_perproc * m%nx
    ! Offsets for all procs
    nnoffset_glob_perproc = nwoffset_glob_perproc * m%nx
    ! Global number of nodes
    nnGlobal = nwGlobal * m%nx
    m%nnglob = nnGlobal

   ! Allocate global mesh for computing the kernel and for outputs
    allocate(m%xglob(nnGlobal),m%yglob(nnGlobal),&
             m%zglob(nnGlobal),m%dwglob(nwGlobal),m%dipglob(nnGlobal), m%fault_label_glob(nnGlobal), m%restart_slip_glob(nnGlobal))

    ! Gather mesh nodes
    call gather_allvdouble(m%x,   nnLocal, m%xglob,   nnLocal_perproc, nnoffset_glob_perproc, nnGlobal)
    call gather_allvdouble(m%y,   nnLocal, m%yglob,   nnLocal_perproc, nnoffset_glob_perproc, nnGlobal)
    call gather_allvdouble(m%z,   nnLocal, m%zglob,   nnLocal_perproc, nnoffset_glob_perproc, nnGlobal)
    call gather_allvdouble(m%dip, nnLocal, m%dipglob, nnLocal_perproc, nnoffset_glob_perproc, nnGlobal)
    call gather_allvdouble(m%dw,  nwLocal, m%dwglob,  nwLocal_perproc, nwoffset_glob_perproc, nwGlobal)
    call gather_allvi(m%fault_label,  nnLocal, m%fault_label_glob,  nnLocal_perproc, nnoffset_glob_perproc, nnGlobal)
    call gather_allvdouble(m%restart_slip,  nnLocal, m%restart_slip_glob,  nnLocal_perproc, nnoffset_glob_perproc, nnGlobal)

endif

end subroutine init_mesh_2D

end module mesh
