subroutine simpleburgers_ad(xin5,xin,xout5,xout,dt,npdt)

 use const
 use mod_vars
 
 implicit none
 
 complex, dimension (-mm:mm), intent(in)  :: xin5  ! trajectory
 complex, dimension (-mm:mm), intent(out) :: xin
 complex, dimension (-mm:mm), intent(out) :: xout5 ! trajectory
 complex, dimension (-mm:mm), intent(in)  :: xout
 real, intent(in) :: dt
 integer, intent(in) :: npdt
 
 integer :: m, nstep, i
 complex :: ztend_u, ztend_u5
 
 real, dimension(nlon,0:npdt) :: u55 ! for storage of trajectory
 
 complex, dimension (-mm:mm) :: zvar_m ! temporary variable for adjoint of direct FFT
 real, dimension(nlon)       :: zvar   ! temporary variable for adjoint of inverse FF
!
! Initial conditions in spectral space 
!
 u5_m(:,1) = xin5(:)

 do nstep = 0,npdt
! 
! Back to physical space for non-linear terms
! 
   call fft_i(u5_m(:,1),u5)
 
   u25(:) = 0.5*u5(:)*u5(:) 
   u55(:,nstep) = u5(:) 
 
! Direct FFT for current time step 
   
   call fft_d(u25,u25_m)  ! for non linear terms 
        
   do m=-mm,mm
     ztend_u5 = -j*float(m)*u25_m(m)/a
     u5_m(m,2) = u5_m(m,1) + dt*ztend_u5
   enddo  
       
   u5_m(:,1) = u5_m(:,2)
   
 enddo
!
! Store prognostic variable after npdt time steps  
! 
! xout5(:) = u5_m(:,2) 

  xout5(:) = u25_m(:) 
 
! Adjoint computations

 u_m(:,:) = (0.,0.)
 xin(:) = (0.,0.)
 ztend_u = (0.,0.)
 u(:) = 0.0
 u2(:) = 0.0
 u2_m(:) = (0.,0.)

! u_m(:,2) = u_m(:,2) + xout(:)

 u2_m(:) = u2_m(:) + xout(:)
 
 do nstep=npdt,0,-1
 
   u_m(:,2) = u_m(:,2) + u_m(:,1) 
   u_m(:,1) = (0.,0.) 
 
   !do m=-mm,mm
   !  ztend_u  = ztend_u + dt*u_m(m,2)
   !  u_m(m,1) = u_m(m,1) + u_m(m,2)
   !  u_m(m,2) = (0.,0.)
   !  !u2_m(m) = u2_m(m) + ztend_u*conjg(-j*float(m)/a) 
   !  ztend_u = (0.,0.)    
   !enddo  
   
   !call fft_i(u2_m,zvar)  ! for non linear terms 
   !u2(:) = u2(:) + zvar(:)/float(nlon)
   !zvar(:) = 0.0
   
   !u(:) = u(:) + u2(:)*u55(:,nstep)
   !u2(:) = 0.0
   
   !call fft_d(u,zvar_m)
   !u_m(:,1) = u_m(:,1) + zvar_m(:)*float(nlon)
   !zvar_m(:) = (0.,0.)
 
 enddo 
  
 xin(: ) = xin(:) + u_m(:,1) 
 u_m(:,1) = (0.,0.) 
 
 return

end subroutine simpleburgers_ad
