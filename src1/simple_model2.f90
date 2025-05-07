subroutine simpleburgers2(xin,xout,dt,npdt)

 use const
 use mod_vars
 
 implicit none
 
 complex, dimension (-mm:mm), intent(in)  :: xin
 complex, dimension (-mm:mm), intent(out) :: xout
 real, intent(in) :: dt
 integer, intent(in) :: npdt
 
 integer :: m, nstep, j1
 complex :: ztend_u
 complex, dimension(-mm:mm) :: filter
!
! Initial conditions in spectral space 
!
 u_m(:,1) = xin(:)

 do nstep = 0,npdt
! 
! Back to physical space for non-linear terms
! 
  !call fft_i(u_m(:,1),u)
   
   do j1=1,nlon
     u(j1) = 0.0
     do m=-mm,mm
       u(j1) = u(j1) + real(u_m(m,1)*cexp(2.0*j*pi*float(m)*(float(j1) - 1.0)/float(nlon)))
     enddo
   enddo   
 
   u2(:) = 0.5*u(:)*u(:) 
 
! Direct FFT for current time step 
   
  !call fft_d(u2,u2_m)  ! for non linear terms 
   
   do m=-mm,mm
     u2_m(m) = 0.0
     do j1=1,nlon
       u2_m(m) =  u2_m(m) + u2(j1)*cexp(-2.0*j*float(m)*pi*(float(j1) - 1.0)/float(nlon)) 
     enddo
     u2_m(m) = u2_m(m)/float(nlon)
   enddo
        
   do m=-mm,mm
     ztend_u = -j*float(m)*u2_m(m)/a
     u_m(m,1) = (u_m(m,1) + dt*ztend_u)/(1.0 + dt*(float(m)/a)**2*kdiff)
   enddo  
   
 enddo
!
! Store prognostic variable after npdt time steps  
! 
 xout(:) = u_m(:,1) 
 
 return

end subroutine simpleburgers2
