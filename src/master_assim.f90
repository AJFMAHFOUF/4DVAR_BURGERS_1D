program master_assim
!--------------------------------------------------------------------------------------------
! Main program to perform a 4D-Var incremental variational assimilation
! with a 1D Burgers model in spectral space (advection-diffusion equation)
!--------------------------------------------------------------------------------------------
! Simulated observations are taken from a reference state (with a spatial sampling)
! The default assimilation window is 24 h and a subsequent 24-h forecast is also performed
! The background error covariance matrix is a SOAR model (default)
! 
! Experimental set-up taken from Rabier and Liu (2003) (ECMWF Annual Seminar Proceedings)
!
!                                                                  J.-F. Mahfouf (05/2025)
!--------------------------------------------------------------------------------------------
 
 use const
 use mod_vars
 use spectral_vars
 use fft99_mod
 
 implicit none
 
 integer :: i, ii, m, npdt0

 integer    :: islot  ! current time slot
 
 complex, dimension(-mm:mm) :: xin, xout    ! prognostic variable (input/output)
 complex, dimension(-mm:mm) :: xin5, xout5  ! prognostic variable (input/output) - initial trajectory
 complex, dimension(-mm:mm) :: xin6, xout6  ! prognostic variable (input/output) - final trajectory
 complex, dimension(-mm:mm) :: xout7, xout8 
 complex, dimension(-mm:mm) :: xout9        ! prognostic variable (input/output) - model prediction
 complex, dimension(-mm:mm) :: xint, xoutt  ! prognostic variable (true value) 
 complex, dimension(-mm:mm) :: xadd
 complex, dimension(-mm:mm) :: chi          ! control vector in spectral space
 complex, dimension(-mm:mm) :: gradientm    ! gradient of cost-function is spectral space
 complex, dimension(-mm:mm) :: dkm          ! for descent direction
 
 complex, allocatable       :: xin55(:,:)
 real, allocatable          :: yo5(:,:), yot(:,:), d0(:,:), d0s(:,:)
 
 real :: xj, xjb, xjo, betak, zgr1, t1, t2, xv 
 
 call  cpu_time(time=t1)
!
! Define various experimental set-ups (B and R matrices) + initial conditions
!
 call init
 
! Allocate arrays depending upon number of observations and time-slots 
 
 allocate (yo5(nobs,0:nslots))
 allocate (yot(nobs,0:nslots))
 allocate (d0(nobs,0:nslots))
 allocate (d0s(nobs,0:nslots))
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
 call burgers(xin5,xout5,npdt0)
 call burgers(xint,xoutt,npdt0)
!
! Back in physical space
! 
 call fft_i(xout5,u5_f)
 call fft_i(xoutt,ut_f) 
 
! Create the necessary observations 
 
 call create_obs(xin5,yo5,.false.)
 call create_obs(xint,yot,.true.)
 
! Save initial conditions of trajectory for the various time slots
 
 xin55(:,1) = xin5(:)
 
 do islot=1,nslots
   call burgers(xin5,xout5,npdt)
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
! 
! Compute the initial gradient  
!
 do islot=nslots,1,-1  
   d0(:,islot) = yot(:,islot) - yo5(:,islot)  ! innovation vector
   d0s(:,islot) = d0(:,islot)/sigmao**2       ! scaled innovation
   call hopt_ad(xin55(:,islot),xin,yo5(:,islot),d0s(:,islot))
   xin(:) = xin(:) + xadd(:) 
   call burgers_ad(xin55(:,islot),xout,xout5,xin,npdt)
   xadd(:) = xout(:)
 enddo 
 
 d0(:,0) = yot(:,0) - yo5(:,0)      ! innovation vector
 d0s(:,0) = d0(:,0)/sigmao**2       ! scaled innovation
 call hopt_ad(xin55(:,1),xin,yo5(:,0),d0s(:,0))
 xout(:) = xout(:) + xin(:) 
 
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
! Minimisation algorithm - conjugate gradient descent method
!  
 call congrad(chi,gradientm,xin55,d0)
!
! From control vector to spectral space : dx = L^{-1}chi - final increment
! 
 call chavarin(chi,xin6)

 xin6(:) = xin6(:) + xin5(:)
!
! Initial states in physical space
!
 call fft_i(xin6,u6_i)
 call fft_i(xin5,u5_i)
 call fft_i(xint,ut_i) 
!
! Call model for final trajectory
!
 call burgers(xin6,xout6,npdt0)
!
! Back in physical space
! 
 call fft_i(xout6,u6_f)
! 
! Need to recover output from background state: xout5 - overwritten in DA process
!
 call fft_d(u5_f,xout5)  
! 
! Subsequent 24-h predictions after 24-h 
! 
 call burgers(xout6,xout7,npdt0) ! initial state from 4D-Var assimilation
 call burgers(xoutt,xout8,npdt0) ! true initial state
 call burgers(xout5,xout9,npdt0) ! initial state from background field
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
   write (172,*) xv,u5_i(i),ut_i(i),u6_i(i)
 enddo 

! Deallocate arrays  
 
 deallocate (yo5)
 deallocate (yot)
 deallocate (d0)
 deallocate (d0s)
 deallocate (xin55)
 
 call cpu_time(time=t2)
 
 print *,'Total execution of 4D-Var assimilation',t2-t1
 
end program master_assim
