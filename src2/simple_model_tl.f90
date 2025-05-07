subroutine simpleburgers_tl(xin5,xin,xout5,xout,dt,npdt)

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
 
!
! Initial conditions in spectral space 
!
 u5_m(:,1) = xin5(:)
 u_m(:,1) = xin(:)

 do nstep = 0,npdt
! 
! Back to physical space for non-linear terms
! 
   call fft_i(u5_m(:,1),u5)
   call fft_i(u_m(:,1),u)
 
   u25(:) = 0.5*u5(:)*u5(:) 
   u2(:) = u(:)*u5(:) 
 
! Direct FFT for current time step 
   
   call fft_d(u25,u25_m)  ! for non linear terms 
   call fft_d(u2,u2_m)  ! for non linear terms 
        
   do m=-mm,mm
     !ztend_u = -j*float(m)*20.0*u_m(m,1)/a  ! advection with constant velocity
     ztend_u5 = 0.! -j*float(m)*u25_m(m)/a
     ztend_u  = 0.! -j*float(m)*u2_m(m)/a
     u5_m(m,2) = u5_m(m,1) + dt*ztend_u5
     u_m(m,2)  = u_m(m,1)  + dt*ztend_u
   enddo  
       
   u5_m(:,1) = u5_m(:,2)
   u_m(:,1)  = u_m(:,2)
   
 enddo
!
! Store prognostic variable after npdt time steps  
! 
! xout5(:) = u5_m(:,2) 
! xout(:)  = u_m(:,2) 

  xout5(:) = u25_m(:)
  xout(:)  = u2_m(:) 
 
 return

end subroutine simpleburgers_tl
