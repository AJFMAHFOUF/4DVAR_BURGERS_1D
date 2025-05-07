program master_assim

 use const
 use mod_vars
 use spectral_vars
 use fft99_mod
 
 implicit none
 
 integer :: i, ii, m, k, npdt0
 
 real       :: length ! forecast duration in hours
 integer    :: islot  ! current time slot
 
 complex, dimension(-mm:mm) :: xin, xout    ! prognostic variable (input/output)
 complex, dimension(-mm:mm) :: xin5, xout5  ! prognostic variable (input/output) - initial trajectory
 complex, dimension(-mm:mm) :: xin6, xout6  ! prognostic variable (input/output) - final trajectory
 complex, dimension(-mm:mm) :: xout7, xout8, xout9  ! prognostic variable (input/output) - model prediction
 complex, dimension(-mm:mm) :: xint, xoutt  ! prognostic variable (true value) 
 complex, dimension(-mm:mm) :: xin0, xadd
 complex, dimension(-mm:mm) :: chi          ! control vector in spectral space
 complex, dimension(-mm:mm) :: dx           ! increment in spectral space 
 complex, dimension(-mm:mm) :: eta_m        ! random number in spectral space 
 complex, dimension(-mm:mm) :: grad1, grad2 ! parts of gradient
 complex, dimension(-mm:mm) :: gradient, gradientm ! gradient of cost-function is spectral space
 complex, dimension(-mm:mm) :: dk, dkm, zdkt     ! for descent direction
 complex, dimension(-mm:mm) :: zvar, zvar5
 !complex, dimension(-mm:mm,nslots) :: xin55      ! initial conditions for trajectory at each time slot
 
 !real, dimension(nobs,nslots) :: yo5, yot, d0, d0s, hdx, eta_o, yo_save
 
 !real, dimension(nobs)        :: yo
 complex, allocatable    :: xin55(:,:)
 real, allocatable       :: yo5(:,:), yot(:,:), d0(:,:), d0s(:,:), yo_save(:,:) 
 
 real :: zsum, xj, xjb, xjo, alphak, betak, zgr1, zgr2, t1, t2, xv 
!
! Define various experimental set-ups (B and R matrices) + initial conditions
!
 call init
 
! Allocate arrays depending upon number of observations and time-slots 
 
 allocate (yo5(nobs,nslots))
 allocate (yot(nobs,nslots))
 allocate (d0(nobs,nslots))
 allocate (d0s(nobs,nslots))
 allocate (yo_save(nobs,nslots))
 allocate (xin55(-mm:mm,nslots))
!
! Index "t" is the truth and index "5" is the background 
! 
 xint(:) = ut_m(:)
 xin5(:) = u5_m(:)
 
 npdt0 = int(nlength/dt) - 1  
!
! Call model for initial and true trajectories (over the whole assimilation window)
!
 call burgers(xin5,xout5,dt,npdt0)
 call burgers(xint,xoutt,dt,npdt0)
!
! Back in physical space
! 
 call fft_i(xout5,u5_f)
 call fft_i(xoutt,ut_f) 
 
! Set-up initial gradient to zero 

 gradientm(:) = (0.,0.)
 
! Create the necessary observations 
 
 call create_obs(xin5,yo5,.false.)
 call create_obs(xint,yot,.true.)
 
! Save initial conditions of trajectory for the various time slots
 
 xin55(:,1) = xin5(:)
 
 do islot=1,nslots
   call burgers(xin5,xout5,dt,npdt)
   xin5(:) = xout5(:)
   if (islot /= nslots) xin55(:,islot+1) = xin5(:)
 enddo 
!
! Restore initial conditions 
! 
 xin5(:) = xin55(:,1) 
!
! Initial conditions for adjoint integrations
! 
 xin(:)  = (0.,0.)
 xadd(:) = (0.0,0.0)
 
! Compute the initial gradient  

 do islot=nslots,1,-1  
   d0(:,islot) = yot(:,islot) - yo5(:,islot)
   d0s(:,islot) = d0(:,islot)/sigmao**2 
   call hopt_ad(xin55(:,islot),xin,yo5(:,islot),d0s(:,islot))
   xin(:) = xin(:) + xadd(:) 
   call burgers_ad(xin55(:,islot),xout,xout5,xin,npdt)
   xadd(:) = xout(:)
 enddo 
 
 call chavarin_ad(gradientm,xout)
 
 zgr1 = 0.0      
 do m=-mm,mm
   zgr1 = zgr1 + conjg(gradientm(m))*gradientm(m) 
 enddo  
 
 chi(:) = (0.,0.) 
 dkm(:) = (0.,0.)
 betak = 0.0
 
 call cost_function(chi,d0,xin5,xj,xjo,xjb)
   
 write (173,*) 0,xj,zgr1 
!
! Start minimisation algorithm - conjugate gradient descent method
!  
 do k=1,20  ! start iteration for conjugate gradient 
 
! Define a descent direction

   dk(:) = gradientm(:) + betak*dkm(:)  
  
! Define the step in the descent direction

   zvar5(:) = xin5(:)  
   
   call chavarin(dk,zvar)
   
   do islot=1,nslots  
     call burgers_tl(zvar5,zvar,xout5,xout,npdt)
     call hopt_tl(xout5,xout,yo5(:,islot),yo_save(:,islot))
     zvar5(:) = xout5(:) 
     zvar(:)  = xout(:)             
     yo_save(:,islot) = yo_save(:,islot)/sigmao**2   
   enddo
   
   xin(:) = (0.0,0.0)
   xadd(:) = (0.0,0.0)
     
   do islot=nslots,1,-1     
     call hopt_ad(xin55(:,islot),xin,yo5(:,islot),yo_save(:,islot))
     xin(:) = xin(:) + xadd(:) 
     call burgers_ad(xin55(:,islot),xout,xout5,xin,npdt)
     xadd(:) = xout(:)  
   enddo
   
   call chavarin_ad(zdkt,xout)

   zsum = 0.0
   zgr1 = 0.0   
   
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
  
   call cost_function(chi,d0,xin5,xj,xjo,xjb)  

   write (173,*) k,xj,zgr2
 
 enddo  ! end of loop - conjugate gradient descent algorithm
!
! From control vector to spectral space : dx = L^{-1}chi
! 
 call chavarin(chi,xin6)

 xin6(:) = xin6(:) + xin5(:)
!
! Call model for final trajectory
!
 call burgers(xin6,xout6,dt,npdt0)
!
! Back in physical space
! 
 call fft_i(xout6,u6_f)
! 
! Subsequent 24-h predictions after 24-h 
!
 call burgers(xout6,xout7,dt,npdt0) ! initial state from 4D-Var assimilation
 call burgers(xoutt,xout8,dt,npdt0) ! true initial state
 call burgers(xout5,xout9,dt,npdt0) ! initial state from background field
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

! Deallocate arrays  
 
 deallocate (yo5)
 deallocate (yot)
 deallocate (d0)
 deallocate (d0s)
 deallocate (yo_save)
 deallocate (xin55)
 
 call cpu_time(time=t2)
 
 print *,'Total execution of 4D-Var assimilation',t2-t1
 
end program master_assim
