program master_tlad3

 use const
 use mod_vars
 use spectral_vars
 use fft99_mod
 
 implicit none
 
 integer :: i, m
 
 real    :: length ! forecast duration in hours
 real    :: dt ! model time step
 integer :: npdt ! number of time steps
 real    :: t1, t2 ! execution time = t2 - t1
 
 integer, parameter :: islot = 1
 
 complex, dimension(-mm:mm) :: xin, xout    ! prognostic variable (input/output)
 complex, dimension(-mm:mm) :: xin5, xout5  ! prognostic variable (input/output) - Trajectory
 
 complex, dimension(-mm:mm) :: xin0, xin50
 
 real, dimension(nobs,nslots) :: yo5, yo 
 
 real, dimension(nlon)        :: sigmab
 
 real :: scal1, scal2, zsum, c1
 
!
! Background error covariance matrix
! 
 sigmab(:) = 2.0
!
! Preliminary computations for spectral correlations
! 
 zsum = 0.0
 do m=1,mm
   zsum = zsum + 1./(1.0 + (float(m)*L/a)**2)**2
 enddo
!
! Normalisation factor 
! 
 c1 = (1.0 + 2.0*float(mm))/(1.0 + 2.0*zsum)  
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
 
 dt= 180.0*2
 length = 48.0
 npdt = int(length*3600.0/dt) - 1
  
 xin5(:) = u5_m(:,1)
 xin50(:) = xin5(:)
 xin(:)  = u_m(:,1)
 
 call cpu_time(time=t1)
! 
! Non-linear model  
! 
 call burgers(xin5,xout5,dt,npdt)
 
 call fft_i(xout5,u5_f)
 
 do i=1,nlon
   write (40,*) i,u5_i(i),u5_f(i)
 enddo 
!
! Tangent-linear model  
!  
 xout5 = (0.,0.)
 xout  = (0.,0.)
 
 xin0(:) = xin(:) 
 
 call hopL_tl(xin5,xin,yo5,yo,sigmab,c1,dt,npdt,islot) 
 
 scal1 = 0.0
 do i=1,nobs
   scal1 = scal1 + yo(i,1)**2
 enddo
! 
! Adjoint model
!
 call hopL_ad(xin5,xin,yo5,yo,sigmab,c1,dt,npdt,islot) 
 
 scal2 = 0.0
 do m=-mm,mm
   scal2 = scal2 + real(xin0(m))*real(xin(m)) + aimag(xin0(m))*aimag(xin(m))
 enddo
 
 print *,'Test adjoint HOPL'
 
 print *,'scal1=',scal1,'scal2=',scal2,'ratio=',scal1/scal2
 
 call fft_i(xout,u_f)
 call fft_i(xout5,u5_f)
 
 do i=1,nlon
   write (50,*) i,u5_i(i),u5_f(i),u_i(i),u_f(i)
 enddo 
 
 call cpu_time(time=t2)
 
 print *,'execution of HOPL model',t2-t1
 
end program master_tlad3
