program master_tlad

 use const
 use mod_vars
 use spectral_vars
 use fft99_mod
 
 implicit none
 
 integer :: i, m, ii
 
 real    :: length ! forecast duration in hours
 real    :: dt ! model time step
 integer :: npdt ! number of time steps
 real    :: t1, t2 ! execution time = t2 - t1
 
 complex, dimension(-mm:mm) :: xin, xout    ! prognostic variable (input/output)
 complex, dimension(-mm:mm) :: xin5, xout5  ! prognostic variable (input/output) - Trajectory
 
 complex, dimension(-mm:mm) :: xin0, xin50
 
 complex, dimension(-mm:mm) :: eta_m
 
 real, dimension(nlon) :: sigmab, eta
 
 real :: scal1, scal2, zsum, Lb, c1
 
 Lb = 208.0E3
 
 sigmab(:) = 2.0 
!
! Initialisation for FFT991
! 
 call set99(trigs,ifax,nlon)
! 
! Define initial conditions in physical space
!
 do i=1,nlon
   ii = i - nlon/2
   u5_i(i) = -20.0*sin(2.*pi*(float(ii)-1.0)/float(nlon)) 
 enddo 
 
 do i=1,nlon
   call gasdev(eta(i))
   eta(i) = eta(i)*sigmab(i) 
 enddo  
 
 write (100,*) eta
 
! Scaled perturbed state vector in spectral space
  
 call fft_d(eta,eta_m)

! Perturbation scaled by correlation matrix in spectral space

 zsum = 0.0
 do m=1,mm
   zsum = zsum + 1.0/(1.0 + (float(m)*Lb/a)**2)**2
 enddo
!
! Normalisation factor 
! 
 c1 = (1.0 + 2.0*float(mm))/(1.0 + 2.0*zsum)
 
 do m=-mm,mm
   eta_m(m) = eta_m(m)*sqrt(c1)/(1.0 + (float(m)*Lb/a)**2)
 enddo
 
 u_m(:,1) = eta_m(:)  
 
! 
! Define initial conditions in spectral space
! 
 call fft_d(u5_i,u5_m(:,1))  
 
 dt= 600.0
 length = 24.0
 npdt = int(length*3600.0/dt) - 1
 
 xin5(:) = u5_m(:,1)
 xin50(:) = xin5(:)
 xin(:)  = u_m(:,1)
 
 call cpu_time(time=t1)
! 
! Non-linear model  
! 
 call simpleburgers(xin5,xout5,dt,npdt)
 
 call fft_i(xout5,u5_f)
 
 do i=1,nlon
   write (45,*) i,u5_i(i),u5_f(i)
 enddo 
 
 call burgers(xin5,xout5,dt,npdt)
 
 call fft_i(xout5,u5_f)
 
 do i=1,nlon
   write (46,*) i,u5_i(i),u5_f(i)
 enddo 
! 
! Tangent-linear model  
!  
 xout5 = (0.,0.)
 xout  = (0.,0.)
 
! xin(:) = (10.,0.0)
 
 xin0(:) = xin(:) 
 
 call simpleburgers_tl(xin5,xin,xout5,xout,dt,npdt) 
 
 call fft_i(xout5,u5_f)
 
 do i=1,nlon
   write (44,*) i,u5_i(i),u5_f(i)
 enddo 
 
 scal1 = 0.0
 do m=-mm,mm
  ! scal1 = scal1 + real(xout(m))**2 + aimag(xout(m))**2
    scal1 = scal1 + 2.0*xout(m)*conjg(xout(m))
 enddo
! 
! Adjoint model
! 
 xin5(:) = xin50(:) ! should not be necessary - input array not changed after routine call
 
 call simpleburgers_ad(xin5,xin,xout5,xout,dt,npdt)
 
 call fft_i(xout5,u5_f)
 
 call fft_i(xin,u_i)  ! output from adjoint (in physical space)
 
 do i=1,nlon
   write (43,*) i,u5_i(i),u5_f(i)
 enddo 
 
 scal2 = 0.0
 do m=-mm,mm
   !scal2 = scal2 + real(xin0(m))*real(xin(m)) + aimag(xin0(m))*aimag(xin(m))
   scal2 = scal2 + xin0(m)*conjg(xin(m)) + conjg(xin0(m))*xin(m)
 enddo
 
 print *,'Test adjoint BURGERS'
 
 print *,'scal1=',scal1,'scal2=',scal2,'ratio=',scal1/scal2
 
 call fft_i(xout,u_f)
 call fft_i(xout5,u5_f)
 
 do i=1,nlon
   write (70,*) i,u5_i(i),u5_f(i),u_i(i),u_f(i)
 enddo 
 
 call cpu_time(time=t2)
 
 print *,'execution of Burgers model',t2-t1
 
end program master_tlad
