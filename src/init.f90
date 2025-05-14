subroutine init
!----------------------------------------------------------------------------
! Definition of experimental set-ups:
! * Background error variances (correlation model defined in correl.f90)
! * Background correlation length 
! * Observation error variances (assumption of uncorrelated observations)
! * Model time step
! * Length of assimilation window
! * Number of time slots over the assimilation window
! * Number of observations per time slot
! * Set-up FFT arrays
! * Initial state and perturbed state (according to B^{1/2} matrix)
!
! External routines 
!  - gasdev : random noise according to Gaussian distribution
!  - correl : correlation model in spectral space (function)
!
!                                                  J.-F. Mahfouf (05/2025)
!-----------------------------------------------------------------------------

 use const
 use mod_vars
 use spectral_vars
 use fft99_mod
 
 implicit none
 
 complex, dimension(-mm:mm) :: eta_m
 real, dimension(nlon)      :: eta
 real                       :: zsum, correl
 integer                    :: m, i, ii
!
! Background error standard deviation 
!
 sigmab(:) = 2.0 
!
! Background error correlation length
!
 Lb = 208.0E3
!
! Observation error standard deviation
!
 sigmao = 1.0
!
! Model time step (s)
! 
 dt = 600.0
!
! Length of assimilation window (s)
! 
 nlength = 24*3600
!
! Number of observations (according to a sampling w.r.t. the model grid)
!
 nsampling = 4
! 
 nobs = nlon/nsampling
!
! Temporal frequency of observations (s)
!
 nfreq = 12*3600 
!
! Number of time steps over a time slot
! 
 npdt = int(nfreq/dt) - 1 
!
! Number of time slots over the assimilation window
! 
 nslots = nlength/nfreq
! 
 print *,'number of time slots =',nslots,' with assim window =',nlength/3600,'hours'
!
! Preliminary computations for spectral correlations
! 
 zsum = 0.0
 do m=1,mm
   zsum = zsum + correl(m)
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
 
 call fft_d(ut_i,ut_m(:))  
 
 do i=1,nlon
   call gasdev(eta(i))
   eta(i) = eta(i)*sigmab(i) 
 enddo  
 
! write (100,*) eta !- to reproduce the same random numbers => read the stored file
 
! Scaled perturbed state vector in spectral space
  
 call fft_d(eta,eta_m)
 
 do m=-mm,mm
   eta_m(m) = eta_m(m)*sqrt(c1*correl(m))
 enddo
 
 u5_m(:) = ut_m(:) + eta_m(:)  
 
! Read perturbed state with random numbers for reproductibility 
 
 read (101,*) u5_m(:)
 
 return
 
end subroutine init
