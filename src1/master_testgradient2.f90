program master_testgradient

 use const
 use mod_vars
 use spectral_vars
 use fft99_mod
 
 implicit none
 
 integer :: i, m, k, ii
 
 real    :: length ! forecast duration in hours
 real    :: dt ! model time step
 integer :: npdt ! number of time steps
 real    :: t1, t2 ! execution time = t2 - t1
 integer :: islot ! current time slot
 
 complex, dimension(-mm:mm) :: xin, xout    ! prognostic variable (input/output)
 complex, dimension(-mm:mm) :: xin5, xout5  ! prognostic variable (input/output) - initial trajectory
 complex, dimension(-mm:mm) :: xin6, xout6  ! prognostic variable (input/output) - final trajectory
 complex, dimension(-mm:mm) :: xint, xoutt  ! prognostic variable (true value) 
 complex, dimension(-mm:mm) :: xin0
 complex, dimension(-mm:mm) :: chi, chi2    ! control vector in spectral space
 complex, dimension(-mm:mm) :: dx           ! increment in spectral space 
 complex, dimension(-mm:mm) :: grad1, grad2 ! parts of gradient
 complex, dimension(-mm:mm) :: gradient, gradientm ! gradient of cost-function is spectral space
 complex, dimension(-mm:mm) :: dk, dkm, adk        ! descent directions
 complex, dimension(-mm:mm) :: eta_m
 
 real, dimension(nobs,nslots) :: yo5, yot, d0, hdx
 
 real, dimension(nobs)        :: d0s, yo
 
 real, dimension(nlon)        :: zvar, sigmab, chi_pdg, innov, dx_pdg, eta
 
 real :: scal1, scal2, zsum, c1, sigmao, xj1, xj2, xjb, xjo, alphak, betak, zgr1, zgr2, Lb 

! B matrix specifications 
 
 sigmab(:) = 2.0 
!
! R matrix specification
!
 sigmao = 1.0 
!
! Correlation length in physical space
! 
 Lb = 208.0E3
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
 print *,'number of longitudes',nlon,' number of observations',nobs
 
 call set99(trigs,ifax,nlon)
! 
! Define true and perturbed initial conditions in physical space 
! 
 do i=1,nlon
   ii = i - nlon/2
   ut_i(i) = -20.0*sin(2.*pi*(float(ii)-1.0)/float(nlon))
 enddo 
! 
! Define initial conditions in spectral space
! 
 call fft_d(ut_i,ut_m(:,1))  

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
 
 u5_m(:,1) = ut_m(:,1) + eta_m(:)
 
 u_m(:,1) = eta_m(:) 
!
 dt= 600.0 
 length = 24.0 ! where an observation set is available
 islot = 1 ! test with one time slot
 
 npdt = int(length*3600.0/dt) - 1
  
 xint(:) = ut_m(:,1)
 xin5(:) = u5_m(:,1)
 
 call cpu_time(time=t1) 
!
! Call model for initial and true trajectories
!
 call simpleburgers(xin5,xout5,dt,npdt)
 call simpleburgers(xint,xoutt,dt,npdt)
!
! Back in physical space
! 
 call fft_i(xout5,u5_f)
 call fft_i(xoutt,ut_f) 

! Define true and perturbed observations

 call hopt(xint,yot,dt,npdt,islot)
 call hopt(xin5,yo5,dt,npdt,islot)
 
! Innovation vector  
 
 d0(:,1) = yot(:,1) - yo5(:,1)
 
 do i=1,nobs
!   print *,'innovation vector for first time slot',d0(i,1),i
 enddo 
  
! length = 12.0 ! where an observation set is available
! islot = 2 ! test with one time slot
! npdt = int(length*3600.0/dt) - 1
 
! call hopt(xint,yot,dt,npdt,islot)
 
! Define initial increment in spectral space 

 dx(:) = u_m(:,1)
 
! Change of variable 

 call fft_i(dx,dx_pdg) 
 
 do i=1,nlon
   dx_pdg(i) = dx_pdg(i)/sigmab(i)
 enddo  
 
 call fft_d(dx_pdg,chi)
 
 do m=-mm,mm
   chi(m) = chi(m)*(1.0 + (float(m)*L/a)**2)/sqrt(c1)
 enddo
 
!
! Compute cost-function (1)
! 
 
! 1) Background term

   xjb = 0.0
   do m=-mm,mm
     xjb = xjb + conjg(chi(m))*chi(m)
   enddo
   xjb = 0.5*xjb
 
! 2) Observation term

   call hopL_tl(xin5,chi,yo5,hdx(:,1),sigmab,c1,dt,npdt,islot)
 
   xjo = 0.0
   do i=1,nobs
     xjo = xjo + ((hdx(i,1) - d0(i,1))/sigmao)**2
   enddo
   xjo = 0.5*xjo
 
   xj1 = xjo + xjb  
 
! Compute gradient of cost function : A*chi + b

! 1) Term A*chi

   hdx(:,1) = hdx(:,1)/sigmao**2

   call hopL_ad(xin5,grad1,yo5,hdx(:,1),sigmab,c1,dt,npdt,islot)
 
! 2) Term b

     d0s(:) = d0(:,1)/sigmao**2
     call hopL_ad(xin5,grad2,yo5,d0s,sigmab,c1,dt,npdt,islot)
     
! Total gradient in chi-space

   gradient(:) =  -grad2(:) + chi(:) + grad1(:) !- grad2(:) !+ chi(:)
   
   zgr1 = 0.0
   do m=-mm,mm
     zgr1 = zgr1 + conjg(gradient(m))*gradient(m)
   enddo
   
   do k=0,15 
   
   scal2 = 10.0**(-k)
   
   chi2(:) = chi(:) + scal2*gradient(:)  
!
! Compute cost-function (2)
! 
 
! 1) Background term

   xjb = 0.0
   do m=-mm,mm
     xjb = xjb + conjg(chi2(m))*chi2(m)
   enddo
   xjb = 0.5*xjb
 
! 2) Observation term

   call hopL_tl(xin5,chi2,yo5,hdx(:,1),sigmab,c1,dt,npdt,islot)
 
   xjo = 0.0
   do i=1,nobs
     xjo = xjo + ((hdx(i,1) - d0(i,1))/sigmao)**2
   enddo
   xjo = 0.5*xjo
 
   xj2 = xjo + xjb   
   
   scal1 = (xj2 -xj1)/(scal2*zgr1)
   
   print *,'ratio',scal1,scal2,'FD x1=',(xj2 -xj1),'GR x2=',scal2*zgr1!,'chi=',chi(1)

 enddo
 
 stop
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
! Back to spectral space
!   
 call fft_d(chi_pdg,xin6)
!
! Call model for final trajectory
!
 call simpleburgers(xin6,xout6,dt,npdt)
!
! Back in physical space
! 
 call fft_i(xout6,u6_f)

 do i=1,nlon
   write (80,*) i,u5_f(i),ut_f(i),u6_f(i)
 enddo 
 
 call cpu_time(time=t2)
 
 print *,'Total execution of 4D-Var assimilation',t2-t1
 
end program master_testgradient
