subroutine burgers_ad(xin5,xin,xout5,xout,npdt1)
!---------------------------------------------------------------------------
! 1D Burgers equation solved in spectral space -  adjoint version
!	
! Temporal schemes : Leapfrog for advection + RAW filter 
!                    Euler backward (implicit) for diffusion
!
! Inputs : 
! - xin5  : spectral coefficients at initial time (trajectory)
! - xout  : gradient w.r.t. spectral coefficients at final time (adjoint)
! - npdt1 : number of time steps for temporal integration 
! 
! Outputs :
! - xout5 : spectral coefficients after npdt1 time steps (trajectory)
 !- xin   : gradient w.r.t. spectral coefficients at initial time (adjoint)
!
! Transform method for non-linear terms
!
! External routines:
! 
! - fft_d : direct spectral transform
! - fft_i : inverse spectral transform
!
!                                                J.-F. Mahfouf (05/2025)
!-----------------------------------------------------------------------------
 use const
 use mod_vars
 
 implicit none
 
 complex, dimension (-mm:mm), intent(in)  :: xin5  ! trajectory
 complex, dimension (-mm:mm), intent(out) :: xin
 complex, dimension (-mm:mm), intent(out) :: xout5 ! trajectory
 complex, dimension (-mm:mm), intent(in)  :: xout
 integer, intent(in) :: npdt1
 
 integer :: m, nstep, i
 complex :: ztend_u, ztend_u5
 complex, dimension (-mm:mm) :: filter, filter5
 
 real, dimension(nlon,0:npdt) :: u55 ! for storage of trajectory
 
 complex, dimension (-mm:mm) :: zvar_m ! temporary variable for adjoint of direct FFT
 real, dimension (nlon)      :: zvar   ! temporary variable for adjoint of inverse FFT

!-------------------------------
!  Trajectory computations 
!-------------------------------

!
! Initial conditions in spectral space 
!
 u5_m(:,1) = xin5(:)
  
 do nstep = 0,npdt1
!
! Back to physical space for non linear terms 
! 
   if (nstep == 0) then
     call fft_i(u5_m(:,1),u5)
   else
     call fft_i(u5_m(:,2),u5)
   endif
   
   u25(:) = 0.5*u5(:)**2
   u55(:,nstep) = u5(:) ! storage of trajectory
!
!  Direct FFT for non linear terms
!
   call fft_d(u25,u25_m)    
 
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
       ztend_u5 = -j*float(m)*20.0*u5_m(m,1)/a
       u5_m(m,2) = (u5_m(m,1) + dt*ztend_u5)/(1.0 + dt*(float(m)/a)**2*kdiff)
     enddo  
       
   endif     

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
 xin(:) = (0.,0.)
 ztend_u = (0.,0.)
 filter(:) = (0.,0.) 
 u(:) = 0.0
 u2(:) = 0.0
 u2_m(:) = (0.,0.)
 zvar(:) = 0.0
 zvar_m(:) = (0.,0.)
!
! Store prognostic variable after npdt time steps  
! 
 u_m(:,2) = u_m(:,2) + xout(:)
! xout(:) = (0.,0.) - not used after (therefore not modified by the subroutine)
! 
 do nstep = npdt1,0,-1
 
     if (nstep > 0) then
 
! Swap variable time step

       u_m(:,3) = u_m(:,3) + u_m(:,2)
       u_m(:,2) = (0.,0.) 
       u_m(:,2) = u_m(:,2) + u_m(:,1)
       u_m(:,1) = (0.,0.) 
       
! Apply Robert Asselin Williams filter to remove 2*dt noise
       
       filter(:) = filter(:) - nu*(1.0-wk)*u_m(:,3)
       filter(:) = filter(:) + nu*wk*u_m(:,2)
       u_m(:,1)  = u_m(:,1) + filter(:)
       u_m(:,2)  = u_m(:,2) - 2.0*filter(:)
       u_m(:,3)  = u_m(:,3) + filter(:)
       filter(:) = (0.,0.)
!
! Solve advection equation in spectral space (leapfrog scheme)
!
       do m=-mm,mm
         u_m(m,1) = u_m(m,1) + u_m(m,3)/(1.0 + 2.0*dt*(float(m)/a)**2*kdiff) 
         ztend_u  = ztend_u + 2.0*dt*u_m(m,3)/(1.0 + 2.0*dt*(float(m)/a)**2*kdiff)
         u_m(m,3) = (0.,0.)
         u2_m(m)  = u2_m(m) + ztend_u*conjg(-j*float(m)/a)
         ztend_u  = (0.,0.)
         u_m(m,2) = u_m(m,2) + conjg(-j*float(m)*20.0)*ztend_u/a
         ztend_u  = (0.,0.)
       enddo
       
    else    
     
       do m=-mm,mm
         u_m(m,1) = u_m(m,1) + u_m(m,2)/(1.0 + dt*(float(m)/a)**2*kdiff) 
         ztend_u  = ztend_u + dt*u_m(m,2)/(1.0 + dt*(float(m)/a)**2*kdiff)
         u_m(m,2) = (0.,0.)
         u2_m(m)  = u2_m(m) + ztend_u*conjg(-j*float(m)/a)
         ztend_u  = (0.,0.)
         u_m(m,1) = u_m(m,1) + conjg(-j*float(m)*20.0)*ztend_u/a
         ztend_u  = (0.,0.)
       enddo
         
    endif 
    
    call fft_i(u2_m,zvar) 
    u2(:) = u2(:) + zvar(:)/float(nlon)
    u2_m(:) = (0.,0.)
    
    u(:) = u(:) + u2(:)*u55(:,nstep)
    u2(:) = 0.0
    
    if (nstep == 0) then
      call fft_d(u,zvar_m)
      u_m(:,1) = u_m(:,1) + zvar_m(:)*float(nlon)
      u(:) = 0.
    else
      call fft_d(u,zvar_m)
      u_m(:,2) = u_m(:,2) + zvar_m(:)*float(nlon)
      u(:) = 0.
    endif
    
 enddo  
!   
! Initial conditions in spectral space 
!
 xin(:) = xin(:) + u_m(:,1) 
 u_m(:,1) = (0.,0)
  
 return

end subroutine burgers_ad
