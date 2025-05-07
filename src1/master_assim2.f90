program master_assim

 use const
 use mod_vars
 use spectral_vars
 use fft99_mod
 
 implicit none
 
 integer :: i, m, k
 
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
 complex, dimension(-mm:mm) :: chi, chis    ! control vector in spectral space
 complex, dimension(-mm:mm) :: dx           ! increment in spectral space 
 complex, dimension(-mm:mm) :: grad1, grad2 ! parts of gradient
 complex, dimension(-mm:mm) :: gradient, gradientm ! gradient of cost-function is spectral space
 complex, dimension(-mm:mm) :: dk, dkm, adk        ! descent directions
 
 real, dimension(nobs,nslots) :: yo5, yot, d0, hdx
 
 real, dimension(nobs)        :: d0s, yo
 
 
 real, dimension(nlon)        :: zvar, sigmab, chi_pdg, innov, dx_pdg
 
 real :: scal1, scal2, zsum, c1, sigmao, xj, xjb, xjo, alphak, betak, zgr1, zgr2 

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
   ut_i(i) = -20.0*cos(2.0*pi*(float(i)-1.0)/float(nlon))
   u5_i(i) = -10.0*cos(2.0*pi*(float(i)-1.0)/float(nlon))
 enddo 
! 
! Define initial conditions in spectral space
! 
 call fft_d(ut_i,ut_m(:,1))  
 call fft_d(u5_i,u5_m(:,1))  
!
 dt= 600.0 
 length = 12.0 ! where an observation set is available
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
 
! do i=1,nobs
!   print *,'innovation vector for first time slot',d0(i,1),i
! enddo 
  
! length = 12.0 ! where an observation set is available
! islot = 2 ! test with one time slot
! npdt = int(length*3600.0/dt) - 1
 
! call hopt(xint,yot,dt,npdt,islot)
 
! Define initial increment in spectral space 

 chi(:) = (0.,0.) 
 
! Initial gradient  
 
 d0s(:) = d0(:,1)/sigmao**2
 
 call hopL_ad(xin5,grad2,yo5,d0s,sigmab,c1,dt,npdt,islot)
 
 gradientm(:) = grad2(:)
 
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

   call hopL_tl(xin5,chi,yo5,hdx(:,1),sigmab,c1,dt,npdt,islot)
 
   xjo = 0.0
   do i=1,nobs
     xjo = xjo + ((hdx(i,1) - d0(i,1))/sigmao)**2
   enddo
   xjo = 0.5*xjo
 
   xj = xjo + xjb   
   
   print *,'value of cost-function before iterations',xj,xjo,xjb 
 
 
 do k=1,20  ! start iteration for conjugate gradient 
 
! Define a descent direction

  dk(:) = gradientm(:) + betak*dkm(:)  
  
! Define the step in the descent direction

  call hopL_tl(xin5,dk,yo5,yo,sigmab,c1,dt,npdt,islot) 
  
  yo(:) = yo(:)/sigmao**2
  
  call hopL_ad(xin5,adk,yo5,yo,sigmab,c1,dt,npdt,islot)

   zsum = 0.0
   zgr1 = 0.0
   do m=-mm,mm
     zsum = zsum + conjg(dk(m))*(dk(m) + adk(m))
     zgr1 = zgr1 + conjg(gradientm(m))*gradientm(m) 
   enddo  
   alphak =  zgr1/zsum    
   
! Define a new value of the state vector 
 
   chi(:) = chi(:) + alphak*dk(:)   
 
! New gradient and new beta factor for descent direction

   gradient(:) = gradientm(:) - alphak*(dk(:) + adk(:))
   
   zgr1 = 0.0
   zgr2 = 0.0
   do m=-mm,mm
     zgr1 = zgr1 + conjg(gradient(m))*gradient(m) 
     zgr2 = zgr2 + conjg(gradientm(m))*gradientm(m) 
   enddo  
   betak = zgr1/zgr2
!   print *,'gradient norm',k,zgr1,zgr2,betak,alphak
   if (betak > 1.0) then
     print *,'Exit from loop since gradient norm is increasing'
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

   call hopL_tl(xin5,chi,yo5,hdx(:,1),sigmab,c1,dt,npdt,islot)
 
   xjo = 0.0
   do i=1,nobs
     xjo = xjo + ((hdx(i,1) - d0(i,1))/sigmao)**2
   enddo
   xjo = 0.5*xjo
 
   xj = xjo + xjb   
   
   print *,'value of cost-function at iteration',k,' ',xj,xjo,xjb 
 
   do m=-mm,mm
     chis(m) = chi(m)*sqrt(c1)/(1.0 + (float(m)*L/a)**2)
   enddo
!
! Back in physical space 
! 
  call fft_i(chis,chi_pdg)
!
! Scale by standard deviation of background error and add background field
! 
   do i=1,nlon
    chi_pdg(i) = chi_pdg(i)*sigmab(i) + u5_i(i) 
   enddo 
 
 enddo
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
! Scale by standard deviation of background error and add background field
! 
 do i=1,nlon
   chi_pdg(i) = chi_pdg(i)*sigmab(i) + u5_i(i) 
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
 
end program master_assim
