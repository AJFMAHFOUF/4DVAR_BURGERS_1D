subroutine burgers_ad(xin5,xin,xout5,xout,npdt1)
!------------------------------------------------------------------------
! 1D Burgers equation solved in spectral space -  adjoint version
!	
! Temporal schemes : Euler forward for advection 
!                    Euler backward for diffusion
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
!-------------------------------------------------------------------------
 use const
 use mod_vars
 
 implicit none
 
 complex, dimension (-mm:mm), intent(in)   :: xin5  ! trajectory
 complex, dimension (-mm:mm), intent(out)  :: xin
 complex, dimension (-mm:mm), intent(out)  :: xout5 ! trajectory
 complex, dimension (-mm:mm), intent(in)   :: xout
 integer, intent(in) :: npdt1
 
 integer :: m, nstep
 complex :: ztend_u5, ztend_u
 
 complex, dimension(-mm:mm) :: zvar_m
 real, dimension(nlon)      :: zvar
 
 real, dimension(nlon,0:npdt) :: u55
 
!
! Initial conditions in spectral space 
!
 u5_m(:) = xin5(:)

 do nstep = 0,npdt1
! 
! Back to physical space for non-linear terms
! 
   call fft_i(u5_m(:),u5)
 
   u25(:) = 0.5*u5(:)*u5(:) 
   u55(:,nstep) = u5(:) 
 
! Direct FFT for current time step 
   
   call fft_d(u25,u25_m)  ! for non linear terms 
        
   do m=-mm,mm
     ztend_u5 = -j*float(m)*u25_m(m)/a
     u5_m(m) = (u5_m(m) + dt*ztend_u5)/(1.0 + dt*(float(m)/a)**2*kdiff)
   enddo  
   
 enddo
!
! Store prognostic variable after npdt1 time steps  
! 
 xout5(:) = u5_m(:) 
 
! Adjoint computations

 u_m(:) = (0.,0.)
 ztend_u = (0.,0.)
 u2_m(:) = (0.,0.) 
 u2(:) = 0.0
 u(:) = 0.0
 xin(:) = (0.,0.)
! 
 
 u_m(:) = u_m(:) + xout(:)
 !xout(:) = (0.,0.)
 
 do nstep = npdt1,0,-1
    
    do m=-mm,mm
      ztend_u = ztend_u + dt*u_m(m)/(1.0  + dt*(float(m)/a)**2*kdiff)
      u_m(m) = u_m(m)/(1.0  + dt*(float(m)/a)**2*kdiff)
      u2_m(m) = u2_m(m) + ztend_u*j*float(m)/a
      ztend_u = (0.,0.)
    enddo  
 
    call fft_i(u2_m,zvar)
    u2(:) = u2(:) + zvar(:)/float(nlon)
    u2_m(:) = (0.,0.)
      
    u(:)  = u(:) + u2(:)*u55(:,nstep)
    u2(:) = 0.0   
 
    call fft_d(u,zvar_m)
    u_m(:) = u_m(:) + zvar_m(:)*float(nlon)
    u(:) = 0.0
  
 enddo
 
 xin(:) = xin(:) + u_m(:)
 u_m(:) = (0.,0.)
 
 return

end subroutine burgers_ad
