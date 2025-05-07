subroutine congrad(chi,gradientm,xin55,d0)
!-----------------------------------------------------------------------------------------
! Conjugate gradient minimisation algorithm (spectral space)
! suitable for a quadratic cost-function J (e.g. incremental 4D-Var)
!
! Call of:
!
! - cost_function : Diagnostic of J at each iteration (optional)
! - burgers_tl    : Dynamical model for the time propagation of increments
! - burgers_ad    : Adjoint of dynamical model for retroprogation of gradients
! - hopt_tl       : Conversion of spectral increments into observation space
! - hopt_ad       : Adjoint of observation operator = conversion of gradient from
!                   observation space to model space (including spectral transform)
!
! Inputs:
! - chi           : control vector in transformed space (chi = Ldx 
!                   where L is the square-root of B^{-1})
! - gradientm     : initial gradient of the cost-function
! - xin55         : initial conditions for each timeslot over the assimilation window
! - d0            : innovation vector = difference between simulated and 
!                   real observations for each timeslot
!
! Outputs:
! - chi           : control vector at the end of the minimisation (analysis state)
! - gradientm     : final gradient of the cost-function 
!
!  The four last routines are necessary to compute the optimal step in the
!  descent direction and the gradient of the cost-function with respect 
!  to the initial conditions at the beginning of the assimilation window
!
!
!                                                           J.-F. Mahfouf (05/2025)
!-------------------------------------------------------------------------------------------
 use const
 use mod_vars
 use spectral_vars
 
 implicit none
 
 complex, dimension(-mm:mm), intent(inout)     :: chi       ! control vector in spectral space
 complex, dimension(-mm:mm), intent(inout)     :: gradientm ! gradient of cost-function
 complex, dimension(-mm:mm,nslots), intent(in) :: xin55     ! initial conditions for each time slot
 complex, dimension(-mm:mm,nslots), intent(in) :: d0        ! innovation vector
 
 integer    :: islot  ! current time slot
 
 complex, dimension(-mm:mm) :: xin,  xout     ! prognostic variable (input/output)
 complex, dimension(-mm:mm) :: xin5, xout5    ! prognostic variable (input/output) - initial trajectory
 complex, dimension(-mm:mm) :: xadd
 complex, dimension(-mm:mm) :: gradient       ! gradient of cost-function is spectral space
 complex, dimension(-mm:mm) :: dk, dkm, zdkt  ! for descent directions
 complex, dimension(-mm:mm) :: zvar           ! temporary variable
 
 real, dimension(nobs,nslots) :: yo_save
 
 real, dimension(nobs) :: yo5
 
 real :: zsum, xj, xjb, xjo, alphak, betak, zgr1, zgr2
 
 integer :: niter_max, m, k
!
 dkm(:) = (0.,0.)
 betak = 0.0
 niter_max = 20
!
! Initial conditions at the beginning of the assimilation window
! 
 xin5(:) = xin55(:,1)
 
 call cost_function(chi,d0,xin5,xj,xjo,xjb)
!
! Start minimisation algorithm - conjugate gradient descent method
!  
 do k=1,niter_max  ! start iteration for conjugate gradient 
 
! Define a descent direction

   dk(:) = gradientm(:) + betak*dkm(:)  
  
! Define the step in the descent direction
   
   call chavarin(dk,zvar)
   
   do islot=1,nslots  
     call burgers_tl(xin55(:,islot),zvar,xout5,xout,npdt)
     call hopt_tl(xout5,xout,yo5,yo_save(:,islot))
     zvar(:) = xout(:)             
     yo_save(:,islot) = yo_save(:,islot)/sigmao**2   
   enddo
   
   xin(:) = (0.0,0.0)
   xadd(:) = (0.0,0.0)
     
   do islot=nslots,1,-1     
     call hopt_ad(xin55(:,islot),xin,yo5,yo_save(:,islot))
     xin(:) = xin(:) + xadd(:) 
     call burgers_ad(xin55(:,islot),xout,xout5,xin,npdt)
     xadd(:) = xout(:)  
   enddo
   
   call chavarin_ad(zdkt,xout)

   zsum = 0.0
   zgr1 = 0.0   
   
   do m=-mm,mm
     zsum = zsum + conjg(dk(m))*(dk(m) + zdkt(m))
     zgr1 = zgr1 + conjg(gradientm(m))*gradientm(m) 
   enddo  
   alphak = zgr1/zsum   
   
! Define a new value of the state vector 
 
   chi(:) = chi(:) + alphak*dk(:)   
 
! New gradient and new beta factor for descent direction

   gradient(:) = gradientm(:) - alphak*(dk(:) + zdkt(:))  
  
   zgr2 = 0.0
   do m=-mm,mm
     zgr2 = zgr2 + conjg(gradient(m))*gradient(m) 
   enddo  
   betak = zgr2/zgr1 
!
! Swapp gradient and direction 
!   
   gradientm(:) = gradient(:)
   dkm(:) = dk(:)
  
   call cost_function(chi,d0,xin5,xj,xjo,xjb)  

   write (173,*) k,xj,zgr2
 
 enddo 
 
 return
end subroutine congrad
