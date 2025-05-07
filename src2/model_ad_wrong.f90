subroutine burgers_ad(xin5,xin,xout5,xout,dt,npdt)

 use const
 use mod_vars
 
 implicit none
 
 complex, dimension (-mm:mm), intent(in)  :: xin5  ! trajectory
 complex, dimension (-mm:mm), intent(out) :: xin
 complex, dimension (-mm:mm), intent(out) :: xout5 ! trajectory
 complex, dimension (-mm:mm), intent(in)  :: xout
 real, intent(in) :: dt
 integer, intent(in) :: npdt
 
 integer :: m, nstep
 complex :: ztend_u5, ztend_u
 complex, dimension(-mm:mm) :: filter5, filter
 
 real, dimension(nlon,0:npdt) :: u55 ! for storage of trajectory
 
!----------------------------
!  Trajectory computations
!----------------------------
 
!
! Initial conditions in spectral space 
!
 u5_m(:,1) = xin5(:)
! 
! Back to physical space for non-linear terms
! 
 call fft_i(u5_m(:,1),u5)

 do nstep = 0,npdt
 
   u25(:) = 0.5*u5(:)*u5(:) 
   
! Store trajectory

   u55(:,nstep) = u5(:) 
 
! Direct FFT for current time step 
   
   call fft_d(u25,u25_m)  ! for non linear terms 
  
   if (nstep > 0) then
   
! Solve Burgers equation in spectral space with implicit diffusion   (leapfrog)
   
     do m=-mm,mm
       ztend_u5 = -j*float(m)*u25_m(m)/a
       u5_m(m,3) = (u5_m(m,1) + 2.0*dt*ztend_u5)/(1.0 + 2.0*dt*(float(m)/a)**2*kdiff)
     enddo
 
! Apply Robert Asselin Williams filter to remove 2*dt noise
     
     filter5(:) = u5_m(:,1) - 2.0*u5_m(:,2) + u5_m(:,3)
     u5_m(:,2)  = u5_m(:,2) + nu*wk*filter5(:)
     u5_m(:,3)  = u5_m(:,3) - nu*(1.0-wk)*filter5(:)  
  
! Swap time steps    
   
     u5_m(:,1) = u5_m(:,2)
     u5_m(:,2) = u5_m(:,3)
        
   else
   
     do m=-mm,mm
       ztend_u5 = -j*float(m)*u25_m(m)/a
       u5_m(m,2) = (u5_m(m,1) + dt*ztend_u5)/(1.0 + dt*(float(m)/a)**2*kdiff)
     enddo  
       
   endif     

! Inverse FFT for next time step   
   
   call fft_i(u5_m(:,2),u5)
   
 enddo
!
! Store prognostic variable after npdt time steps  
! 
 xout5(:) = u5_m(:,2) 
 
!-----------------------------
!  Adjoint computations
!-----------------------------
! 
! Initialisations 

 u_m(:,:) = (0.,0.)
 ztend_u = (0.,0.)  
 filter(:) = (0.,0.)
 u(:) = 0.0
 u2(:) = 0.0
 u2_m(:) = (0.,0.)
 xin(:) = (0.,0.)
!
! Store prognostic variable after npdt time steps  
! 
 u_m(:,2) = u_m(:,2) + xout(:)
! xout(:) = (0.,0.) - not used after (therefore not modified by the subroutine)
! 
 do nstep = npdt,0,-1
 
! Adjoint of inverse FFT for next time step    

  !call fft_d(u,u_m(:,2))
  !u_m(:,2) = u_m(:,2)*float(nlon)
  
  
  if (nstep > 0) then
  
! Swap time steps

    u_m(:,3) = u_m(:,3) + u_m(:,2)
    u_m(:,2) = (0.,0.)
    u_m(:,2) = u_m(:,2) + u_m(:,1)
    u_m(:,1) = (0.,0.) 

! Apply Robert Asselin Williams filter to remove 2*dt noise
     
     filter(:) = filter(:) - nu*(1.0-wk)*u_m(:,3)
     filter(:) = filter(:) + nu*wk*u_m(:,2)
     u_m(:,1) = u_m(:,1) + filter(:)
     u_m(:,2) = u_m(:,2) - 2.0*filter(:)
     u_m(:,3) = u_m(:,3) + filter(:)
     filter(:) = (0.,0.)
     
! Solve Burgers equation in spectral space with implicit diffusion   (leapfrog)
   
     do m=-mm,mm
       u_m(m,1) = u_m(m,1) + u_m(m,3)/(1.0 + 2.0*dt*(float(m)/a)**2*kdiff)
       ztend_u = ztend_u + 2.0*dt*u_m(m,3)/(1.0 + 2.0*dt*(float(m)/a)**2*kdiff)
       u_m(m,3) = (0.,0.)
       u2_m(m) = u2_m(m) + ztend_u*conjg(-j*float(m)/a) ! conjugate transpose
       ztend_u = (0.,0.)
     enddo  
     
  else
      
    do m=-mm,mm
       u_m(m,1) = u_m(m,1) + u_m(m,2)/(1.0 + dt*(float(m)/a)**2*kdiff)
       ztend_u = ztend_u + u_m(m,2)*dt/(1.0 + dt*(float(m)/a)**2*kdiff)
       u_m(m,2) = (0.,0.)
       u2_m(m) = u2_m(m) + ztend_u*conjg(-j*float(m)/a) ! conjugate transpose
       ztend_u =  (0.,0.)
    enddo  
  
  endif
  
! Adjoint of direct FFT for current time step 
   
   call fft_i(u2_m,u2)  
   u2(:) = u2(:)/float(nlon)
   u2_m(:) = (0.,0.)
   
   u(:) = u(:) + u2(:)*u55(:,nstep)
   u2(:) = 0.0
   
   if (nstep > 0) then
     call fft_d(u,u_m(:,2)) 
     u_m(:,2) = u_m(:,2)*float(nlon)
     u(:) = (0.,0.)
   endif  
 
 enddo
! 
! Back to physical space for non-linear terms (adjoint)
! 
 call fft_d(u,u_m(:,1)) 
 u_m(:,1) = u_m(:,1)*float(nlon)
 u(:) = 0.0
!
! Initial conditions in spectral space 
!
 xin(:) = xin(:) + u_m(:,1) 
 u_m(:,1) = (0.,0)
 
 return

end subroutine burgers_ad
