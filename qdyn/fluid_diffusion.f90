module fluid_diffusion
use problem_class, only : problem_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! FINITE VOLUME METHOD 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   
!   Compute Pressure P for new time using implicit scheme
!   TODO: change ds that is wrong for the FVM (it is the size of element, not the distance between the two points of FVM)
!   global input: permeablity, P
implicit none
private
public :: compute_P 


contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   
!   Compute Pressure P for new time using implicit scheme
!   It is using conjugate gradient method with a sparse matrix representation
!
subroutine compute_P(dt,pb,ier)
type(problem_type), intent(inout) :: pb
double precision      :: dt
integer               :: source_id,iteration_count
integer               :: ier
double precision, dimension(pb%mesh%nn) :: pk,rk,rk1,xk      ! Tmp array for implicit solver
integer                                     :: k                              ! Integer for for-loop  
! 
! Variable related to fluid diffusion
double precision                               :: alpha_k, beta_k,max_imp, errmax,ds,tol_solver      ! Tmp array for implicit solver
double precision, dimension(pb%mesh%nn)  :: b1,diag, Pscal                          ! Tmp array for implicit solver

! Define tol_solver
tol_solver = 1d-8


! Compute gridsize 
ds = pb%mesh%dx(1)

! Initialisation
ier = 0 ! At the beginning there is no problem
!
!######### Compute the permeability based on permeability_star
if (pb%features%var_k == 1) then
    pb%fluid_diff%permeability = (pb%var_k%kstar-pb%var_k%kmin)*exp(-abs(pb%sigma-pb%fluid_diff%P_temp)/abs(pb%var_k%Snk)) + pb%var_k%kmin
endif


!######### Compute the new permeability_x
! Make the harmonic average
pb%fluid_diff%permeability_x(2:pb%mesh%nn) = 2.*pb%fluid_diff%permeability(1:pb%mesh%nn-1) &
                                               *pb%fluid_diff%permeability(2:pb%mesh%nn)/ &
    (pb%fluid_diff%permeability(1:pb%mesh%nn-1)+pb%fluid_diff%permeability(2:pb%mesh%nn))
pb%fluid_diff%permeability_x(1) = 0.0
pb%fluid_diff%permeability_x(pb%mesh%nn+1) = 0.0


!
! Preconditioner
call compute_diag(dt,diag,pb)

!######### Compute the vector b
do k=1,pb%mesh%nn
! Compute vector b1
b1(k) = pb%P(k)/diag(k)
enddo


    ! For each source 
    do source_id=1,pb%fluid_diff%nb_source
    
  
        ! If time of injection
        if ((pb%time+dt>=pb%fluid_diff%t_injection_beg(source_id) ).and.(pb%time+dt< pb%fluid_diff%t_injection_end(source_id))) then
            ! print*,'dt',dt
!             print*,'pb%fluid_diff%Q(source_id)',pb%fluid_diff%Q(source_id)
!             print*,'pb%fluid_diff%phi',pb%fluid_diff%phi(pb%fluid_diff%index_injection(source_id))
!             print*,'pb%fluid_diff%beta',pb%fluid_diff%beta(pb%fluid_diff%index_injection(source_id))
!             print*,'ds',ds
!             print*,'diag(pb%fluid_diff%index_injection(source_id))',diag(pb%fluid_diff%index_injection(source_id))
            b1(pb%fluid_diff%index_injection(source_id)) = & 
                                   b1(pb%fluid_diff%index_injection(source_id)) &
                                 + (dt*pb%fluid_diff%Q(source_id) &
                                 /(pb%fluid_diff%phi(pb%fluid_diff%index_injection(source_id)) &
                                 *pb%fluid_diff%beta(pb%fluid_diff%index_injection(source_id))) &
                                 /ds)/&
                                 diag(pb%fluid_diff%index_injection(source_id))
        end if
    enddo
  
! if (pb%time .gt. 10.0) then
! print*,'dt',dt
! print*,'pb%fluid_diff%Q(pb%fluid_diff%index_injection)',pb%fluid_diff%Q(pb%fluid_diff%index_injection)
! print*,'pb%fluid_diff%phi(pb%fluid_diff%index_injection',pb%fluid_diff%phi(pb%fluid_diff%index_injection)
! print*,'pb%fluid_diff%beta(pb%fluid_diff%index_injection(source_id)',pb%fluid_diff%beta(pb%fluid_diff%index_injection)
! print*,'ds',ds
! print*,'diag(pb%fluid_diff%index_injection',diag(pb%fluid_diff%index_injection)
! 
! print*,'maxval b1',maxval(b1)
! print*,'pb%fluid_diff%Q',pb%fluid_diff%Q
! stop
! endif
!if (pb%time .gt. 10.0) stop


! 
!######### Conjugate gradient
! Initialisation
xk = pb%P
rk = b1-compute_MatrixVector(xk,dt,diag,pb)
!print*,'rk ini', dsqrt(sum(rk**2))/pb%mesh%nn
pk = rk
iteration_count = 0
!
! Tentative of scaling for the error
Pscal = abs(pb%fluid_diff%P_temp) +dt*abs(pb%fluid_diff%P_dot_temp)+ tol_solver
Pscal = 1.0
errmax=maxval(abs(rk(:)/Pscal(:)))/tol_solver !Evaluate accuracy.
!
! 
do while (errmax >= 1.0)

! print*,'rk',dsqrt(sum(rk**2))/pb%mesh%nn

! Count the number of iteration
iteration_count = iteration_count + 1

alpha_k = sum(rk**2)/sum(pk*compute_MatrixVector(pk,dt,diag,pb))

xk = xk + alpha_k*pk
rk1 = rk - alpha_k*compute_MatrixVector(pk,dt,diag,pb)

! Update pk
beta_k = sum(rk1**2)/sum(rk**2)
pk = rk1 + beta_k*pk


rk = rk1

if(iteration_count>200)then
  ier = 1
    !print*,' did no converge (compute_P)','dt = ',dt
  return
end if


errmax=maxval(abs(rk(:)/Pscal(:)))/tol_solver !Evaluate accuracy.

end do



!######### Calculate new P_dot at new time
if (dt<1e-9) then
pb%fluid_diff%P_dot_temp = 0.
else
pb%fluid_diff%P_dot_temp = (xk-pb%P)/dt
endif

! Update P
pb%fluid_diff%P_temp = xk

if (maxval(abs(pb%fluid_diff%P_temp))>huge(1.)) then
 ier = 1
 end if
 
 
 
 if (dt>50) then
  !   print*,'dt',dt
!     print*,'diag',diag
!     print*,'b1',b1
!     print*,'pb%P',pb%P
!     print*,'pb%P_temp',pb%fluid_diff%P_temp
!  



endif
 
 
 
 
 
! 
end subroutine
! 
! 
! 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   
!   Compute Pressure P for new time using implicit scheme
!   Function that return the result of the matrix vector multiplication with A 
!
function compute_MatrixVector(vector,dt,diag,pb)
type(problem_type), intent(inout) :: pb
double precision, dimension(pb%mesh%nn)    ::  compute_MatrixVector,vector
double precision                           :: dt,A1,A2,A3
integer      :: crossing_id, fault_id1, fault_id2,node_id1,node_id2
double precision                                :: permeability_x_crossing_left, permeability_x_crossing_right,ds
double precision, dimension(pb%mesh%nn)  :: diag   ! Tmp array for implicit solver
integer                               :: k

! Calculate gridsize
ds = pb%mesh%dx(1) ! Constant gridsize for 1D fault in 2D medium 
! print*,'ds1',ds

!-----------------------------------------------------------------------------------------
! Without coupling between faults
!-----------------------------------------------------------------------------------------

! Fist and last componants
A2 = 1.+(pb%fluid_diff%permeability_x(2))*dt &
                /(pb%fluid_diff%phi(2)*pb%fluid_diff%eta(2)*pb%fluid_diff%beta(2)*ds**2)
A3 = -(pb%fluid_diff%permeability_x(2))*dt   &
           /(pb%fluid_diff%phi(2)*pb%fluid_diff%eta(2)*pb%fluid_diff%beta(2)*ds**2)
             
   
compute_MatrixVector(1) = (A2*vector(1)+A3*vector(2))/diag(1)
           
A1 = -(pb%fluid_diff%permeability_x(pb%mesh%nn))*dt   &
                               /(pb%fluid_diff%phi(pb%mesh%nn)*pb%fluid_diff%eta(pb%mesh%nn) &
                               *pb%fluid_diff%beta(pb%mesh%nn)*ds**2)
A2 = 1.+(pb%fluid_diff%permeability_x(pb%mesh%nn))*dt &
                                  /(pb%fluid_diff%phi(pb%mesh%nn)*pb%fluid_diff%eta(pb%mesh%nn) &
                                 *pb%fluid_diff%beta(pb%mesh%nn)*ds**2)

compute_MatrixVector(pb%mesh%nn) = (A1*vector(pb%mesh%nn-1)+A2*vector(pb%mesh%nn))/diag(pb%mesh%nn)

! For each element on the diagonal
do k=2,pb%mesh%nn-1
    A1 = -(pb%fluid_diff%permeability_x(k))*dt/(pb%fluid_diff%phi(k)*pb%fluid_diff%eta(k)*pb%fluid_diff%beta(k)*ds**2)
    
    
    A2 = 1.+(pb%fluid_diff%permeability_x(k)+pb%fluid_diff%permeability_x(k+1))   &
                    *dt/(pb%fluid_diff%phi(k)*pb%fluid_diff%eta(k)*pb%fluid_diff%beta(k)*ds**2)


    A3 = -(pb%fluid_diff%permeability_x(k+1))*dt/(pb%fluid_diff%phi(k+1)*pb%fluid_diff%eta(k+1)*pb%fluid_diff%beta(k+1)*ds**2)
    
    
    compute_MatrixVector(k) = (A1*vector(k-1)+A2*vector(k)+A3*vector(k+1))/diag(k)
    
end do 



end function
! 
! 
! 
! 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   
!   Compute Pressure P for new time using implicit scheme
!   Function that return the result of the matrix vector multiplication with A 
!
subroutine compute_diag(dt,diag,pb)
type(problem_type), intent(inout) :: pb
double precision                           :: dt,ds
double precision, dimension(pb%mesh%nn)  :: diag   ! Tmp array for implicit solver
integer                               :: k


!
! Calculate ds
ds = pb%mesh%dx(1)
! print*,'ds2 =',ds

!
! compute_MatrixVector(:) = 0._dp
diag = 0.

! Fist and last componants
diag(1) = 1.+(pb%fluid_diff%permeability_x(2))*dt &
                /(pb%fluid_diff%phi(2)*pb%fluid_diff%eta(2)*(pb%fluid_diff%beta(2))*ds**2)
           

diag(pb%mesh%nn) = 1.+(pb%fluid_diff%permeability_x(pb%mesh%nn))*dt &
                                  /(pb%fluid_diff%phi(pb%mesh%nn)*pb%fluid_diff%eta(pb%mesh%nn)*pb%fluid_diff%beta(pb%mesh%nn)*ds**2)


! For each element on the diagonal
do k=2,pb%mesh%nn-1  
    diag(k) = 1.+(pb%fluid_diff%permeability_x(k)+pb%fluid_diff%permeability_x(k+1))   &
                    *dt/(pb%fluid_diff%phi(k)*pb%fluid_diff%eta(k)*pb%fluid_diff%beta(k)*ds**2)   
end do 




end subroutine

end module