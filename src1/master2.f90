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
 complex, dimension(-mm:mm) :: xin6, xout6  ! prognostic variable (input/output) - Trajectory
 
 complex, dimension(-mm:mm) :: xin0, xin50, eta_m
 
 real, dimension(nlon) :: eta, sigmab
 
 real :: scal1, scal2, xv, Lb, zsum, c1 
!
 sigmab(:) = 2.0
 
 Lb = 208.0E3

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
 enddo 
! 
! Define initial conditions in spectral space
! 
 call fft_d(u5_i,u5_m(:,1))  

 do ii=1,nlon
   call gasdev(eta(ii))
   eta(ii) = eta(ii)*sigmab(ii) 
 enddo 
 
! Scaled perturbed state vector in spectral space
  
 call fft_d(eta,eta_m)

! Perturbation scaled by correlation matrix in spectral space

 zsum = 0.0
 do m=1,mm
   zsum = zsum + 1./(1.0 + (float(m)*Lb/a)**2)**2
 enddo
!
! Normalisation factor 
! 
 c1 = (1.0 + 2.0*float(mm))/(1.0 + 2.0*zsum)
 
 do m=-mm,mm
   eta_m(m) = eta_m(m)*sqrt(c1)/(1.0 + (float(m)*Lb/a)**2)
 enddo 
 
 dt= 600.0
 length = 24.0
 npdt = int(length*3600.0/dt) - 1
 
 print *,'Courant number=',20.0*dt/(2*pi*a/float(nlon)),2*pi*a/nlon/1000
 
 xin5(:) = u5_m(:,1)
 
 xin6(:) = u5_m(:,1) + eta_m(:)
 
 read (101,*) xin6  ! read perturbed field in order to have the same IC in terms of random numbers
 
 xin(:) = xin6(:) - xin5(:)
 
 call fft_i(xin,u_i)
!
! Tangent-linear model
! 
 call simpleburgers_tl(xin5,xin,xout5,xout,dt,npdt)
 
 call fft_i(xout,u_f)
 
 call fft_i(xout5,u5_f)
 
 call fft_i(xin6,u6_i)
 
! xin50(:) = xin5(:)
! xin(:)  = u_m(:,1)
! 
! Non-linear model  
! 
 call simpleburgers(xin6,xout6,dt,npdt)
 
 call fft_i(xout6,u6_f)
 
 do i=1,nlon
   ii = i - nlon/2
   xv = 2.0*pi*(float(ii)-1.0)/float(nlon)*a/1.E3 ! horizontal coordinate in km
   write (40,*) xv,u5_i(i),u5_f(i)
   write (50,*) xv,u6_i(i),u6_f(i)
   write (60,*) xv,u_i(i),u_f(i)
 enddo 
 
end program master
