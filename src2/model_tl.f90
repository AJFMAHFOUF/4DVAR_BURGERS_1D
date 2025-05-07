subroutine burgers_tl(xin5,xin,xout5,xout,dt,npdt)

 use const
 use mod_vars
 
 implicit none
 
 complex, dimension (-mm:mm), intent(in)  :: xin5  ! trajectory
 complex, dimension (-mm:mm), intent(in)  :: xin
 complex, dimension (-mm:mm), intent(out) :: xout5 ! trajectory
 complex, dimension (-mm:mm), intent(out) :: xout
 real, intent(in) :: dt
 integer, intent(in) :: npdt
 
 integer :: m, nstep, i
 complex :: ztend_u5, ztend_u
 complex, dimension(-mm:mm) :: filter5, filter
!
! Initial conditions in spectral space 
!
 u5_m(:,1) = xin5(:)
 u_m(:,1)  = xin(:)

 do nstep = 0,npdt
! 
! Back to physical space for non-linear terms
! 
   if (nstep == 0) then
     call fft_i(u5_m(:,1),u5)
     call fft_i(u_m(:,1),u)
   else
     call fft_i(u5_m(:,2),u5)
     call fft_i(u_m(:,2),u)
   endif  
 
   u25(:) = 0.5*u5(:)*u5(:) 
   u2(:)  = u(:)*u5(:)
 
! Direct FFT for current time step 
   
   call fft_d(u25,u25_m)  ! for non linear terms 
   call fft_d(u2,u2_m)
  
   if (nstep > 0) then
   
! Solve Burgers equation in spectral space with implicit diffusion   (leapfrog)
   
     do m=-mm,mm
       !ztend_u = -j*float(m)*20.0*u_m(m,1)/a ! advection with constant velocity
       ztend_u5 = -j*float(m)*u25_m(m)/a
       ztend_u  = -j*float(m)*u2_m(m)/a
       ztend_u5 = -j*float(m)*20.0*u5_m(m,2)/a
       ztend_u  = -j*float(m)*20.0*u_m(m,2)/a
       u5_m(m,3) = (u5_m(m,1) + 2.0*dt*ztend_u5)/(1.0 + 2.0*dt*(float(m)/a)**2*kdiff)
       u_m(m,3)  = (u_m(m,1)  + 2.0*dt*ztend_u)/(1.0  + 2.0*dt*(float(m)/a)**2*kdiff)
     enddo
 
! Apply Robert Asselin Williams filter to remove 2*dt noise
     
     filter5(:) = u5_m(:,1) - 2.0*u5_m(:,2) + u5_m(:,3)
     filter(:)  = u_m(:,1)  - 2.0*u_m(:,2)  + u_m(:,3)
     u5_m(:,2)  = u5_m(:,2) + nu*wk*filter5(:)
     u_m(:,2)   = u_m(:,2)  + nu*wk*filter(:)
     u5_m(:,3)  = u5_m(:,3) - nu*(1.0-wk)*filter5(:)  
     u_m(:,3)   = u_m(:,3)  - nu*(1.0-wk)*filter(:)  
  
! Swap time steps    
   
     u5_m(:,1) = u5_m(:,2)
     u_m(:,1)  = u_m(:,2)
     u5_m(:,2) = u5_m(:,3)
     u_m(:,2)  = u_m(:,3)
        
   else
   
     do m=-mm,mm
       !ztend_u = -j*float(m)*20.0*u_m(m,1)/a  ! advection with constant velocity
       ztend_u5  = -j*float(m)*u25_m(m)/a
       ztend_u   = -j*float(m)*u2_m(m)/a
       ztend_u5 = -j*float(m)*20.0*u5_m(m,1)/a
       ztend_u  = -j*float(m)*20.0*u_m(m,1)/a
       u5_m(m,2) = (u5_m(m,1) + dt*ztend_u5)/(1.0 + dt*(float(m)/a)**2*kdiff)
       u_m(m,2)  = (u_m(m,1)  + dt*ztend_u)/(1.0  + dt*(float(m)/a)**2*kdiff)
     enddo  
       
   endif     

 enddo
!
! Store prognostic variable after npdt time steps  
! 
 xout5(:) = u5_m(:,2) 
 xout(:) = u_m(:,2) 
 
 return

end subroutine burgers_tl
