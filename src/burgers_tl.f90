subroutine burgers_tl(xin5,xin,xout5,xout,npdt1)
!------------------------------------------------------------------------
! 1D Burgers equation solved in spectral space - tangent-linear version
!
! Temporal schemes : Euler forward for advection 
!                    Euler backward for diffusion
!
! Inputs : 
! - xin5  : spectral coefficients at initial time (trajectory)
! - xin   : spectral coefficients at initial time (perturbation)
! - npdt1 : number of time steps for temporal integration 
! 
! Outputs :
! - xout5 : spectral coefficients after npdt1 time steps (trajectory)
 !- xout  : spectral coefficients after npdt1 time steps (perturbation)
!
! Transform method for non-linear terms
!
! External routines:
! 
! - fft_d : direct spectral transform
! - fft_i : inverse spectral transform
!
!                                                J.-F. Mahfouf (05/2025)
!-------------------------------------------------------------------------
 use const
 use mod_vars
 
 implicit none
 
 complex, dimension (-mm:mm), intent(in)  :: xin5  ! trajectory
 complex, dimension (-mm:mm), intent(in)  :: xin
 complex, dimension (-mm:mm), intent(out) :: xout5 ! trajectory
 complex, dimension (-mm:mm), intent(out) :: xout
 integer, intent(in) :: npdt1
 
 integer :: m, nstep
 complex :: ztend_u5, ztend_u
 
!
! Initial conditions in spectral space 
!
 u5_m(:) = xin5(:)
 u_m(:)  = xin(:)

 do nstep = 0,npdt1
! 
! Back to physical space for non-linear terms
! 
   call fft_i(u5_m(:),u5)
   call fft_i(u_m(:),u)
 
   u25(:) = 0.5*u5(:)*u5(:) 
   u2(:)  = u(:)*u5(:) 
 
! Direct FFT for current time step 
   
   call fft_d(u25,u25_m)  ! for non linear terms 
   call fft_d(u2,u2_m)    ! for non linear terms 
        
   do m=-mm,mm
     ztend_u5 = -j*float(m)*u25_m(m)/a
     ztend_u  = -j*float(m)*u2_m(m)/a
     u5_m(m) = (u5_m(m) + dt*ztend_u5)/(1.0 + dt*(float(m)/a)**2*kdiff)
     u_m(m)  = (u_m(m)  + dt*ztend_u)/(1.0  + dt*(float(m)/a)**2*kdiff)
   enddo  
   
 enddo
!
! Store prognostic variable after npdt1 time steps  
! 
 xout5(:) = u5_m(:) 
 xout(:)  = u_m(:) 
 
 return

end subroutine burgers_tl
