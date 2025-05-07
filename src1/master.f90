program master

 use const
 use mod_vars
 use spectral_vars
 use fft99_mod
 
 implicit none
 
 integer :: i, ii, m
 
 real    :: length ! forecast duration in hours
 real    :: dt ! model time step
 integer :: npdt ! number of time steps
 real    :: t1, t2 ! execution time = t2 - t1
 
 complex, dimension(-mm:mm) :: xin, xout    ! prognostic variable (input/output)
 complex, dimension(-mm:mm) :: xin5, xout5  ! prognostic variable (input/output) - Trajectory
 
 complex, dimension(-mm:mm) :: xin0, xin50
 
 real :: scal1, scal2, xv 
!
! Initialisation for FFT991
! 
 print *,'number of longitudes',nlon
 
 call set99(trigs,ifax,nlon)
! 
! Define initial conditions in physical space
!
 do i=1,nlon
   ii = i - nlon/2
   u5_i(i) = -20.0*sin(2.*pi*(float(ii)-1.0)/float(nlon))
!   u_i(i)   = 0.01*u5_i(i)
 enddo 
! 
! Define initial conditions in spectral space
! 
 call fft_d(u5_i,u5_m(:,1))  
! call fft_d(u_i,u_m(:,1))  
 
 dt= 600.0
 length = 48.0
 npdt = int(length*3600.0/dt) - 1
 
 print *,'Courant number=',20.0*dt/(2*pi*a/float(nlon)),2*pi*a/nlon/1000
 
 xin5(:) = u5_m(:,1)
! xin50(:) = xin5(:)
! xin(:)  = u_m(:,1)
! 
! Non-linear model  
! 
 call simpleburgers2(xin5,xout5,dt,npdt)
 
 call fft_i(xout5,u5_f)
 
 do i=1,nlon
   ii = i - nlon/2
   xv = 2.0*pi*(float(ii)-1.0)/float(nlon)*a/1.E3 ! horizontal coordinate in km
   write (31,*) xv,u5_i(i),u5_f(i)
 enddo 
 
end program master
