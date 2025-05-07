program master

 use const
 use mod_vars
 use spectral_vars
 use fft99_mod
 
 implicit none
 
 integer :: i
 
 real    :: length ! forecast duration in hours
 real    :: dt ! model time step
 integer :: npdt ! number of time steps
 real    :: t1, t2 ! execution time = t2 - t1
 
 complex, dimension(-mm:mm) :: xin, xout    ! prognostic variable (input/output)
 complex, dimension(-mm:mm) :: xin5, xout5  ! prognostic variable (input/output) - Trajectory
!
! Initialisation for FFT991
! 
 print *,'number of longitudes',nlon
 
 call set99(trigs,ifax,nlon)
! 
! Define initial conditions in physical space
!
 do i=1,nlon
   u5_i(i) = -20.0*cos(2.0*pi*(float(i)-1.0)/float(nlon))
   u_i(i)   = 0.01*u5_i(i)
 enddo 
! 
! Define initial conditions in spectral space
! 
 call fft_d(u5_i,u5_m(:,1))  
 call fft_d(u_i,u_m(:,1))  
 call fft_i(u_m(:,1),u_i)
 
 dt= 180.0 
 length = 72.0
 npdt = int(length*3600.0/dt) - 1
  
 xin5(:) = u5_m(:,1)
 xin(:)  = u_m(:,1)
 
 call cpu_time(time=t1)
! 
! Non-linear model  
! 
 call simpleburgers(xin5,xout5,dt,npdt)
 
 call fft_i(xout5,u5_f)
 
 do i=1,nlon
   write (40,*) i,u5_i(i),u5_f(i)
 enddo 
 
 stop
! 
! Tangent-linear model  
!  
 xout5 = (0.,0.)
 xout  = (0.,0.)
 
 call burgers_tl(xin5,xin,xout5,xout,dt,npdt) 
 
 call fft_i(xout,u_f)
 call fft_i(xout5,u5_f)
 
 do i=1,nlon
   write (50,*) i,u5_i(i),u5_f(i),u_i(i),u_f(i)
 enddo 
 
 call cpu_time(time=t2)
 
 print *,'execution of Burgers model',t2-t1
 
end program master
