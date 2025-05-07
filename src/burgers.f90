subroutine burgers(xin,xout,npdt1)
!----------------------------------------------------------
! 1D Burgers equation solved in spectral space
!
! Temporal schemes : Euler forward for advection 
!                    Euler backward for diffusion
!
! xin  : spectral coefficients at initial time
! xout : spectral coefficients after npdt1 time steps
!
! Transform method for non-linear terms
!
! External routines:
! 
! - fft_d : direct spectral transform
! - fft_i : inverse spectral transform
!
!                             J.-F. Mahfouf (05/2025)
!----------------------------------------------------------
 use const
 use mod_vars
 
 implicit none
 
 complex, dimension (-mm:mm), intent(in)  :: xin
 complex, dimension (-mm:mm), intent(out) :: xout
 integer, intent(in) :: npdt1
 
 integer :: m, nstep
 complex :: ztend_u
!
! Initial conditions in spectral space 
!
 u_m(:) = xin(:)

 do nstep = 0,npdt1
! 
! Back to physical space for non-linear terms
! 
   call fft_i(u_m(:),u)
 
   u2(:) = 0.5*u(:)*u(:) 
 
! Direct FFT for current time step 
   
   call fft_d(u2,u2_m)  ! for non linear terms 
        
   do m=-mm,mm
     ztend_u = -j*float(m)*u2_m(m)/a
     u_m(m) = (u_m(m) + dt*ztend_u)/(1.0 + dt*(float(m)/a)**2*kdiff)
   enddo  
   
 enddo
!
! Store prognostic variable after npdt1 time steps  
! 
 xout(:) = u_m(:) 
 
 return

end subroutine burgers
