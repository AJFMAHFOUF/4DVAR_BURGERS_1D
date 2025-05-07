subroutine simpleburgers_tl2(xin5,xin,xout5,xout,dt,npdt)

 use const
 use mod_vars
 
 implicit none
 
 complex, dimension (-mm:mm), intent(in)  :: xin5  ! trajectory
 complex, dimension (-mm:mm), intent(in)  :: xin
 complex, dimension (-mm:mm), intent(out) :: xout5 ! trajectory
 complex, dimension (-mm:mm), intent(out) :: xout
 real, intent(in) :: dt
 integer, intent(in) :: npdt
 
 integer :: m, nstep, i, j1
 complex :: ztend_u5, ztend_u
 
!
! Initial conditions in spectral space 
!
 u5_m(:,1) = xin5(:)
 u_m(:,1)  = xin(:)

 do nstep = 0,npdt
! 
! Back to physical space for non-linear terms
! 
!   call fft_i(u5_m(:,1),u5)
!   call fft_i(u_m(:,1),u)
   
   do j1=1,nlon
     u5(j1) = 0.0
     u(j1)  = 0.0
     do m=-mm,mm
       u5(j1) = u5(j1) + real(u5_m(m,1)*cexp(2.0*j*pi*float(m)*(float(j1) - 1.0)/float(nlon)))
       u(j1)  = u(j1)  + real(u_m(m,1)*cexp(2.0*j*pi*float(m)*(float(j1) - 1.0)/float(nlon)))
     enddo
   enddo
 
   u25(:) = 0.5*u5(:)*u5(:) 
   u2(:)  = u(:)*u5(:) 
 
! Direct FFT for current time step 
   
!   call fft_d(u25,u25_m)  ! for non linear terms 
!   call fft_d(u2,u2_m)    ! for non linear terms 
   
   do m=-mm,mm
     u25_m(m) = 0.0
     u2_m(m)  = 0.0
     do j1=1,nlon
       u25_m(m) =  u25_m(m) + u25(j1)*cexp(-2.0*j*float(m)*pi*(float(j1) - 1.0)/float(nlon)) 
       u2_m(m)  =  u2_m(m)  + u2(j1)*cexp(-2.0*j*float(m)*pi*(float(j1) - 1.0)/float(nlon)) 
     enddo
     u25_m(m) = u25_m(m)/float(nlon)
     u2_m(m)  = u2_m(m)/float(nlon)
   enddo
        
   do m=-mm,mm
     ztend_u5 = -j*float(m)*u25_m(m)/a
     ztend_u  = -j*float(m)*u2_m(m)/a
     u5_m(m,1) = (u5_m(m,1) + dt*ztend_u5)/(1.0 + dt*(float(m)/a)**2*kdiff)
     u_m(m,1)  = (u_m(m,1)  + dt*ztend_u)/(1.0  + dt*(float(m)/a)**2*kdiff)
   enddo  
   
 enddo
!
! Store prognostic variable after npdt time steps  
! 
 xout5(:) = u5_m(:,1) 
 xout(:)  = u_m(:,1) 
 
 return

end subroutine simpleburgers_tl2
