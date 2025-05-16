program master_assim
!--------------------------------------------------------------------------------------------
! Main program to perform a 3D-Var incremental variational assimilation
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
 
 real, allocatable          :: yo5(:), yot(:,:), d0(:), d0s(:)
 
 real :: xj, xjb, xjo, betak, zgr1, t1, t2, xv 
 
 call  cpu_time(time=t1)
!
! Define various experimental set-ups (B and R matrices) + initial conditions
!
 call init
 
! Allocate arrays depending upon number of observations and time-slots 
 
 allocate (yo5(nobs))
 allocate (yot(nobs,0:nslots))
 allocate (d0(nobs))
 allocate (d0s(nobs))
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
 
 call create_obs(xint,yot,.true.)
!
! Initial conditions for adjoint integrations
! 
 xin(:)  = (0.,0.)
 xadd(:) = (0.0,0.0)
 
 do islot=0,nslots 
! 
! Compute the initial gradient  
!
   call hopt(xin5,yo5)
   
   d0(:) = yot(:,islot) - yo5(:)     ! innovation vector
   d0s(:) = d0(:)/sigmao**2          ! scaled innovation
   
   call hopt_ad(xin5,xout,yo5,d0s)
 
   call chavarin_ad(gradientm,xout)
 
   zgr1 = 0.0      
   do m=-mm,mm
     zgr1 = zgr1 + conjg(gradientm(m))*gradientm(m) 
   enddo  
 
   chi(:) = (0.,0.) 
   dkm(:) = (0.,0.)
   betak = 0.0
 
   call cost_function(chi,d0,xin5,xj,xjo,xjb)
   
   write (173,*) islot,0,xj,zgr1 
!
! Minimisation algorithm - conjugate gradient descent method
!  
   call congrad(chi,gradientm,xin5,d0,islot)
!
! From control vector to spectral space : dx = L^{-1}chi - final increment
! 
   call chavarin(chi,xin6)

   xin6(:) = xin6(:) + xin5(:)
!
! Call model for the next time slot
!
   call burgers(xin6,xout6,npdt)
!   
! Resulting forecast = initial conditions for next time slot
! 
   xin5(:) = xout6(:)
 
 enddo
!
! Back in physical space
! 
 call fft_i(xout6,u6_f)
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
 
 call cpu_time(time=t2)
 
 print *,'Total execution of 3D-Var assimilation',t2-t1
 
end program master_assim
