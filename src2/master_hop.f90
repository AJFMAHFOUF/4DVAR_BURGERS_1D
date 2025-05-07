program master_hop

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
 integer :: islot ! current time slot
 
 complex, dimension(-mm:mm) :: xin, xout    ! prognostic variable (input/output)
 complex, dimension(-mm:mm) :: xin5, xout5  ! prognostic variable (input/output) - Trajectory
 complex, dimension(-mm:mm) :: xin0
 complex, dimension(-mm:mm) :: chi          ! control vector in spectral space
 complex, dimension(-mm:mm) :: gradient     ! gradient of cost-function is spectral space
 
 real, dimension(nobs,nslots) :: yo, yo5
 
 real, dimension(nlon)        :: zvar, sigmab, chi_pdg, innov
 
 real :: scal1, scal2, zsum, c1, sigmao 

! B matrix specifications 
 
 sigmab(:) = 2.0 
!
! R matrix specification
!
 sigmao = 1.0 
! 
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
!
! Arbitrary set-up of chi
! 
 chi(:) = u_m(:,1)
 
 dt= 180.0 
 length = 6.0 ! where an observation set is available
 islot = 1 ! test with one time slot
 
 npdt = int(length*3600.0/dt) - 1
  
 xin5(:) = u5_m(:,1)
 xin(:)  = u_m(:,1)
 
 call cpu_time(time=t1) 

 call hopt(xin5,yo,dt,npdt,islot)
 
! do i=1,nlon/5
!   print *,'output from observation operator first time slot',yo(i,1),i
! enddo 
  
! length = 12.0 ! where an observation set is available
! islot = 2 ! test with one time slot
! npdt = int(length*3600.0/dt) - 1
 
! call hopt(xin5,yo,dt,npdt,islot)
 
! do i=1,nlon/5
!   print *,'output from observation operator 2nd time slot ',yo(i,2),i
! enddo  
 
! Non-linear model  

print *,'xin5=',xin5(10)

 call burgers(xin5,xout5,dt,npdt)
 
print *,'xout5=',xout5(10) 
 
 call fft_i(xout5,u5_f)
 
! Call non-linear observation operator

 print *,'islot=',islot
 
 call hop(u5_f,yo,islot)
 
! do i=1,nlon/5
!   print *,'output from observation operator',yo(i,2),i
! enddo  
 
 do i=1,nlon
   write (40,*) i,u5_i(i),u5_f(i)
 enddo 
! 
! Tangent-linear model and adjoint
!  
 xout5 = (0.,0.)
 xout  = (0.,0.)
 
 xin0(:) = xin(:)
!
! Apply spectral correlation to spectral coefficients
! 
 do m=-mm,mm
   chi(m) = chi(m)*sqrt(c1)/(1.0 + (float(m)*L/a)**2) ! horizontal correlations in spectral space
 enddo
!
! Back in physical space 
! 
 call fft_i(chi,chi_pdg)
!
! Scale by standard deviation of background error
! 
 do i=1,nlon
   chi_pdg(i) = chi_pdg(i)*sigmab(i)
 enddo
!
! Back to spectral space
!   
 call fft_d(chi_pdg,xin)
!
! Time integration of TL model up to observation time
!  

 print *,'xin=',xin(1)
 call hopt_tl(xin5,xin,yo5,yo,dt,npdt,islot)
 
 print*,'after TL',yo5(20,1),yo(20,1)
!
! Apply R-1 to the output  
! 
! yo(:,islot) = yo(:,islot)/sigmao**2

 print *,'nobs=',nobs,nlon/5,islot
 
 scal1 = 0.0
 do i=1,nobs
   scal1 = scal1 + yo(i,islot)**2
 enddo
!
!  Time integration of AD model from obsevation time to initial time 
!
 call hopt_ad(xin5,xin,yo5,yo,dt,npdt,islot)
 print*,'after AD',yo5(20,1)
 
 scal2 = 0.0
 do m=-mm,mm
   scal2 = scal2 + real(xin0(m))*real(xin(m)) + aimag(xin0(m))*aimag(xin(m))
 enddo
 
 print *,'scal1',scal1,'scal2=',scal2,'ratio',scal1/scal2
 
 stop

! Back to physical space 
 
 call fft_i(xin,zvar)
! 
!  Scaling with  standard deviation of background errors
! 
 do i=1,nlon
   zvar(i) = zvar(i)*sigmab(i)  
 enddo
!
! Back to spectral space
!
 call fft_d(zvar,xin)  
!
! Horizontal correlations in spectral space
! 
 do m=-mm,mm
   xin(m) = xin(m)*sqrt(c1)/(1.0 + (float(m)*L/a)**2) 
 enddo 
! 
!First part of the gradient 
!
 gradient(:) = xin(:) 
 
!------------------------------------
 
 innov(:) = innov(:)/sigmao**2
 
 scal1 = 0.0
 do i=1,nlon/5
   scal1 = scal1 + yo(i,islot)**2
 enddo
!
!  Time integration of AD model from obsevation time to initial time 
!
 call hopt_ad(xin5,xin,yo5,innov,dt,npdt,islot)
 
 scal2 = 0.0
 do m=-mm,mm
   scal2 = scal2 + real(xin0(m))*real(xin(m)) + aimag(xin0(m))*aimag(xin(m))
 enddo
 
 print *,'scal1',scal1,'scal2=',scal2,'ratio',scal1/scal2

stop

! Back to physical space 
 
 call fft_i(xin,zvar)
! 
!  Scaling with  standard deviation of background errors
! 
 do i=1,nlon
   zvar(i) = zvar(i)*sigmab(i)  
 enddo
!
! Back to spectral space
!
 call fft_d(zvar,xin)  
!
! Horizontal correlations in spectral space
! 
 do m=-mm,mm
   xin(m) = xin(m)*sqrt(c1)/(1.0 + (float(m)*L/a)**2) 
 enddo  
!
! Second part of the gradient
! 
 gradient(:) = chi(:) + gradient(:) - xin(:)
 
!-------------------------------------------------
 
 print *,'After AD - xin=',xin(10)
 
 do i=1,nlon/5
   print *,'output from AD observation operator',yo(i,islot),yo5(i,islot),i
 enddo
 
 call burgers_tl(xin5,xin,xout5,xout,dt,npdt) 
 
 call fft_i(xout,u_f)
 call fft_i(xout5,u5_f)
 
 do i=1,nlon
   write (50,*) i,u5_i(i),u5_f(i),u_i(i),u_f(i)
 enddo 
 
 call cpu_time(time=t2)
 
 print *,'execution of Burgers model',t2-t1
 
end program master_hop
