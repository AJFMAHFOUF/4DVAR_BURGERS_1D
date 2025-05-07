program master_assim

 use const
 use mod_vars
 use spectral_vars
 use fft99_mod
 
 implicit none
 
 integer :: i, ii, m, k, npdt0
 
 real       :: length ! forecast duration in hours
 real       :: dt     ! model time step
 integer    :: islot  ! current time slot
 
 complex, dimension(-mm:mm) :: xin, xout    ! prognostic variable (input/output)
 complex, dimension(-mm:mm) :: xin5, xout5  ! prognostic variable (input/output) - initial trajectory
 complex, dimension(-mm:mm) :: xin6, xout6  ! prognostic variable (input/output) - final trajectory
 complex, dimension(-mm:mm) :: xout7, xout8, xout9  ! prognostic variable (input/output) - model prediction
 complex, dimension(-mm:mm) :: xint, xoutt  ! prognostic variable (true value) 
 complex, dimension(-mm:mm) :: xin0
 complex, dimension(-mm:mm) :: chi          ! control vector in spectral space
 complex, dimension(-mm:mm) :: dx           ! increment in spectral space 
 complex, dimension(-mm:mm) :: eta_m        ! random number in spectral space 
 complex, dimension(-mm:mm) :: grad1, grad2 ! parts of gradient
 complex, dimension(-mm:mm) :: gradient, gradientm ! gradient of cost-function is spectral space
 complex, dimension(-mm:mm) :: dk, dkm, zdkt     ! for descent direction
 complex, dimension(-mm:mm,nslots) :: adk        ! for descent direction according to time slot
 
 real, dimension(nobs,nslots) :: yo5, yot, yo, d0, d0s, hdx, eta_o
 
 real, dimension(nlon)        :: zvar, sigmab, chi_pdg, innov, dx_pdg, eta
 
 integer, dimension(nslots)   :: npdt ! number of time steps (to reach different time slots)
 
 real :: scal1, scal2, zsum, c1, sigmao, xj, xjb, xjo, alphak, betak, zgr1, zgr2, t1, t2, xv 

 call cpu_time(time=t1)  
 
! B matrix specifications 
 
 sigmab(:) = 2.0 
!
! R matrix specification
!
 sigmao = 1.0
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
 call set99(trigs,ifax,nlon)
! 
! Define true and perturbed initial conditions in physical space
!
 do i=1,nlon
   ii = i - nlon/2
   ut_i(i) = -20.0*sin(2.*pi*(float(ii)-1.0)/float(nlon)) 
 enddo 
 
 call fft_d(ut_i,ut_m(:,1))  
 
 do i=1,nlon
   call gasdev(eta(i))
   eta(i) = eta(i)*sigmab(i) 
 enddo  
 
! write (100,*) eta !- to reproduce the same random numbers => read the stored file
 
! Scaled perturbed state vector in spectral space
  
 call fft_d(eta,eta_m)
 
 do m=-mm,mm
   eta_m(m) = eta_m(m)*sqrt(c1)/(1.0 + (float(m)*L/a)**2)
 enddo
 
 u5_m(:,1) = ut_m(:,1) + eta_m(:)  
 
! Read perturbed state with random numbers for reproductibility 
 
 read (101,*) u5_m(:,1)
 
! call fft_i(u5_m(:,1),u5_i)  - not necessary
!
 dt= 600.0 
 length = 0.0
 
! Set a number of time steps for a final 24h prediction 
 
 npdt0 = int(24*3600.0/dt) - 1
 
 do islot=1,nslots
   length = length + 6.0 ! observations available every 6 hours
   npdt(islot) = int(length*3600.0/dt) - 1
 enddo
  
 xint(:) = ut_m(:,1)
 xin5(:) = u5_m(:,1)
 
 call cpu_time(time=t1) 
!
! Call model for initial and true trajectories
!
 call simpleburgers(xin5,xout5,dt,npdt(nslots))
 call simpleburgers(xint,xoutt,dt,npdt(nslots))
!
! Back in physical space
! 
 call fft_i(xout5,u5_f)
 call fft_i(xoutt,ut_f) 
 
! Set-up initial gradient to zero 
 
 gradientm(:) = (0.,0.)
 
! Define random vector for noisy observations
 
 do islot=1,nslots
   do ii=1,nobs
     call gasdev(eta_o(ii,islot))
     eta_o(ii,islot) = eta_o(ii,islot)*sigmao 
   enddo
 enddo
    
! Save random vector for reproductibility
 
 print *,'nobs=',nobs,'nslots=',nslots
 read (200,*) eta_o
 
 do islot=1,nslots ! loop over time slots

! Define true and perturbed observations

   call hopt(xint,yot,dt,npdt(islot),islot)
   call hopt(xin5,yo5,dt,npdt(islot),islot)
 
! Innovation vector 
      
   d0(:,islot) = yot(:,islot) + eta_o(:,islot) - yo5(:,islot)
 
! Initial gradient  
 
   d0s(:,islot) = d0(:,islot)/sigmao**2
 
   call hopL_ad(xin5,grad2,yo5,d0s,sigmab,c1,dt,npdt(islot),islot)
   
! Sum-up gradients over the various time slots   
 
   gradientm(:) = gradientm(:) + grad2(:)
 
 enddo 
 
 zgr1 = 0.0      
 do m=-mm,mm
   zgr1 = zgr1 + conjg(gradientm(m))*gradientm(m) 
 enddo  
 
 chi(:) = (0.,0.) 
 dkm(:) = (0.,0.)
 betak = 0.0
 
!=========================  
!  Compute cost-function
!=========================   
   
! 1) Background term

 xjb = 0.0
 do m=-mm,mm
   xjb = xjb + conjg(chi(m))*chi(m)
 enddo
 xjb = 0.5*xjb
 
! 2) Observation term
 
 xjo = 0.0
 do islot=1,nslots
   call hopL_tl(xin5,chi,yo5,hdx,sigmab,c1,dt,npdt(islot),islot)
   do i=1,nobs
     xjo = xjo + ((hdx(i,islot) - d0(i,islot))/sigmao)**2
   enddo
 enddo  
 xjo = 0.5*xjo
 
 xj = xjo + xjb   
   
 write (173,*) 0,xj,zgr1 
 
 do k=1,20  ! start iteration for conjugate gradient 
 
! Define a descent direction

   dk(:) = gradientm(:) + betak*dkm(:)  
  
! Define the step in the descent direction

   do islot=1,nslots

     call hopL_tl(xin5,dk,yo5,yo,sigmab,c1,dt,npdt(islot),islot) 
   
     yo(:,islot) = yo(:,islot)/sigmao**2
  
     call hopL_ad(xin5,adk(:,islot),yo5,yo,sigmab,c1,dt,npdt(islot),islot)
  
   enddo

   zsum = 0.0
   zgr1 = 0.0
   
   zdkt(:) = (0.0,0.0)
   do islot=1,nslots
     zdkt(:) = zdkt(:) + adk(:,islot)
   enddo
   
   do m=-mm,mm
     zsum = zsum + conjg(dk(m))*(dk(m) + zdkt(m))
     zgr1 = zgr1 + conjg(gradientm(m))*gradientm(m) 
   enddo  
   alphak =  zgr1/zsum   
   !print *,'alphak=',alphak 
   
! Define a new value of the state vector 
 
   chi(:) = chi(:) + alphak*dk(:)   
 
! New gradient and new beta factor for descent direction

   gradient(:) = gradientm(:) - alphak*(dk(:) + zdkt(:))  
  
   zgr2 = 0.0
   do m=-mm,mm
     zgr2 = zgr2 + conjg(gradient(m))*gradient(m) 
   enddo  
   betak = zgr2/zgr1
   !print *,'gradient norm',k,zgr1,zgr2,betak,alphak
   if (betak > 1.0) then
     !print *,'Warning: gradient norm is increasing !'
     !exit
     !betak = 0.0
   endif  
!
! Swapp gradient and direction 
!   
   gradientm(:) = gradient(:)
   dkm(:) = dk(:)
  
!=========================  
!  Compute cost-function
!=========================   
   
! 1) Background term

   xjb = 0.0
   do m=-mm,mm
     xjb = xjb + conjg(chi(m))*chi(m)
   enddo
   xjb = 0.5*xjb
 
! 2) Observation term

   xjo = 0.0
   do islot=1,nslots 
     call hopL_tl(xin5,chi,yo5,hdx,sigmab,c1,dt,npdt(islot),islot)
     do i=1,nobs
      xjo = xjo + ((hdx(i,islot) - d0(i,islot))/sigmao)**2
     enddo
   enddo  
   xjo = 0.5*xjo
 
   xj = xjo + xjb   
   
   print *,'value of cost-function and gradient at iteration',k,' ',xj,zgr2 
   write (173,*) k,xj,zgr2
 
 enddo  ! end of loop - conjugate gradient descent algorithm
!
! From control vector to spectral space : dx = L^{-1}chi
! 
 do m=-mm,mm
   chi(m) = chi(m)*sqrt(c1)/(1.0 + (float(m)*L/a)**2) 
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
! Back to spectral space and add background field
!   
 call fft_d(chi_pdg,xin6)
 xin6(:) = xin6(:) + xin5(:)
!
! Call model for final trajectory
!
 call simpleburgers(xin6,xout6,dt,npdt(nslots))
!
! Back in physical space
! 
 call fft_i(xout6,u6_f)
! 
! Subsequent 24-h predictions after 24-h 
!
 call simpleburgers(xout6,xout7,dt,npdt0) ! initial state from 4D-Var assimilation
 call simpleburgers(xoutt,xout8,dt,npdt0) ! true initial state
 call simpleburgers(xout5,xout9,dt,npdt0) ! initial state from background field
!
! Back in physical space
! 
 call fft_i(xout7,u7_f)
 call fft_i(xout8,u8_f)
 call fft_i(xout9,u9_f)

 do i=1,nlon
   ii = i - nlon/2
   xv = 2.0*pi*(float(ii)-1.0)/float(nlon)*a/1.E3 ! horizontal coordinate in km
   write (170,*) xv,u5_f(i),ut_f(i),u6_f(i)
   write (171,*) xv,u9_f(i),u8_f(i),u7_f(i)
 enddo 
 
 call cpu_time(time=t2)
 
 print *,'Total execution of 4D-Var assimilation',t2-t1
 
end program master_assim
