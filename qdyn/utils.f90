! Collection of auxiliary functions

module utils

  use problem_class, only : problem_type
  implicit none
  public
 !   double precision, allocatable, save :: Vec
 !   public :: Vec

contains

!===============================================================================
! Helper routine to pack variables for solver
! Storage conventions:
!
! theta = yt(1::pb%neqs)
! v or tau = yt(2::pb%neqs) (depending on friction law)
! + feature specific variables
!
! dtheta/dt = dydt(1::pb%neqs)
! dv/dt or dtau/dt = dydt(2::pb%neqs)
! + feature specific variables
!===============================================================================
subroutine pack(yt, theta, main_var, sigma, theta2, slip, pb, ideriv)

  type(problem_type), intent(inout) :: pb
  double precision, dimension(pb%neqs*pb%mesh%nn), intent(out) :: yt
  double precision, dimension(pb%mesh%nn), intent(in) :: theta, main_var, sigma
  double precision, dimension(pb%mesh%nn), intent(in) :: theta2, slip
  integer :: nmax, ind
  logical, optional :: ideriv


  ! pb%neqs is defined in problem_class.f90
  nmax = pb%neqs*pb%mesh%nn ! The number of equation for the time solver
  
  yt(1:nmax:pb%neqs) = theta
  yt(2:nmax:pb%neqs) = main_var
  yt(3:nmax:pb%neqs) = slip

  ! initialise the index to 3 because of the first three main variables
  ind = 3
  ! For each features that need ODE solver
  if (pb%features%stress_coupling == 1) then
    ind = ind + 1
    yt(ind:nmax:pb%neqs) = sigma
  endif
  if (pb%features%localisation == 1) then
    ind = ind + 1
    yt(ind:nmax:pb%neqs) = theta2
  endif
  
   ! Add for the permeability change. Be careful because pack and unpack are used for both yt and dydt.
  if (pb%features%var_k == 1) then
    ind = ind + 1
    ! If the optional argument is given
    if (present(ideriv)) then
        if (ideriv) then
             yt(ind:nmax:pb%neqs) = pb%var_k%dkstar_dt
        else 
             yt(ind:nmax:pb%neqs) = pb%var_k%kstar
        endif
    ! Default value if there is no ideriv
    else 
       yt(ind:nmax:pb%neqs) = pb%var_k%kstar
    endif
  endif
  
end subroutine pack

!===============================================================================
! Helper routine to unpack variables from solver
!===============================================================================
subroutine unpack(yt, theta, main_var, sigma, theta2, slip, pb, ideriv)

  type(problem_type), intent(inout) :: pb
  double precision, dimension(pb%neqs*pb%mesh%nn), intent(in) :: yt
  double precision, dimension(pb%mesh%nn) :: theta, main_var, sigma
  double precision, dimension(pb%mesh%nn) :: theta2, slip
  integer :: nmax, ind
  logical, optional :: ideriv

  ! pb%neqs is defined in problem_class.f90
  nmax = pb%neqs*pb%mesh%nn
  
  theta = yt(1:nmax:pb%neqs)
  main_var = yt(2:nmax:pb%neqs)
  slip = yt(3:nmax:pb%neqs)

  ! initialise the index to 3 because of the first three main variable
  ind = 3 
  ! For each features that need ODE solver
  if (pb%features%stress_coupling == 1) then
    ind = ind + 1
    sigma = yt(ind:nmax:pb%neqs)
  else
    sigma = pb%sigma
  endif
  if (pb%features%localisation == 1) then
    ind = ind + 1 
    theta2 = yt(ind:nmax:pb%neqs)
  else
    theta2 = 0d0
  endif
  
  
  ! Add for the permeability change. Be careful because pack and unpack are used for both yt and dydt.
  if (pb%features%var_k == 1) then
    ind = ind + 1
    if (present(ideriv)) then
        if (ideriv) then
            pb%var_k%dkstar_dt = yt(ind:nmax:pb%neqs)
        else 
            pb%var_k%kstar = yt(ind:nmax:pb%neqs)
        endif
    ! Default value if there is no ideriv
    else 
        pb%var_k%kstar = yt(ind:nmax:pb%neqs)
    endif
  endif
  

end subroutine unpack

!
!=====================================================================
!
subroutine save_array(x,y,z,V,iproc,typ,nw,nx)

  double precision, dimension(nw*nx), intent(in) :: V
  double precision, dimension(nw*nx), intent(in) :: x,y,z
  character(len=256) :: fileproc
  character(len=16)   :: typ !PG, 'loc' or 'glo'

  integer :: i,iproc,nx,nw

  write(fileproc,'(a,i6.6,a)') 'snap_v_',iproc,typ

  open(101,file=fileproc(1:len_trim(fileproc)),status='replace',form='formatted',action='write')

  do i=1,nw*nx
      write(101,'(4(D15.7))') x(i),y(i),z(i),V(i)
  enddo

  close(101)

end subroutine save_array

!------------------------------------------------


subroutine save_vectorV(x,y,z,V,iproc,typ,nw,nx)

  double precision, dimension(nw,nx), intent(in) :: V
  double precision, dimension(nw,nx), intent(in) :: x,y,z
  character(len=256) :: fileproc
  character(len=16)   :: typ !PG, 'loc' or 'glo'

  integer :: i,j,iproc,nx,nw

  write(fileproc,'(a,i6.6,a)') 'snap_v_',iproc,typ

  open(101,file=fileproc(1:len_trim(fileproc)),status='replace',form='formatted',action='write')

  do i=1,nw
    do j=1,nx
      write(101,'(4(D15.7))') x(i,j),y(i,j),z(i,j),V(i,j)
    enddo
  enddo

  close(101)

end subroutine save_vectorV

! -----------------------------------------------
subroutine save_vector(V,iproc,typ,nw,nx)

  double precision, dimension(nw,nx), intent(in) :: V
  character(len=256) :: fileproc
  character(len=16)   :: typ !PG, 'loc' or 'glo'

  integer :: i,j,iproc,nx,nw

  write(fileproc,'(a,i6.6,a)') 'snap_v_',iproc,typ

  open(101,file=fileproc(1:len_trim(fileproc)),status='replace',form='formatted',action='write')

  do i=1,nw
    do j=1,nx
      write(101,'(D15.7)') V(i,j)
    enddo
  enddo

  close(101)

end subroutine save_vector

! -----------------------------------------------
subroutine save_vector3(V,iproc,typ,nwloc,nwglob,nx)

  double precision, dimension(nwloc,nwglob,nx), intent(in) :: V
  character(len=256) :: fileproc
  character(len=16)   :: typ !PG, 'loc' or 'glo'

  integer :: i,j,k,iproc,nx,nwloc,nwglob

  write(fileproc,'(a,i6.6,a)') 'snap_v_',iproc,typ

  open(101,file=fileproc(1:len_trim(fileproc)),status='replace',form='formatted',action='write')

 do k=1,nwloc
  do i=1,nwglob
    do j=1,nx
      write(101,'(D15.7)') V(k,i,j)
    enddo
  enddo
 enddo

  close(101)

end subroutine save_vector3

end module utils
