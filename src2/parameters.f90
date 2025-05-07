module params

 implicit none

 integer, parameter :: mm = 100               ! maximum wave number
 integer, parameter :: nlon = 2*mm + 1        ! number of longitudes
 integer, parameter :: mmax = (mm+1)*(mm+4)/2 ! number of stored wavenumbers
 integer, parameter :: nfft = 1               ! number of FFT to be done
 integer, parameter :: nobs = nlon/4          ! number of obserations per time slot
 integer, parameter :: nslots = 1             ! number of timeslots over the assimilation window
 
 complex, parameter :: j = (0,1)              ! square root of -1 
 real, parameter    :: a = 6371.22E3          ! Earth radius
 real, parameter    :: pi = acos(-1.0)        ! Pi constant
 real, parameter    :: g = 9.80616            ! Earth gravitational acceleration
 real, parameter    :: omega = 2.0*pi/86164.1 ! Earth angular speed (stellar day)
 real, parameter    :: nu = 0.02, wk = 0.53   ! tunable parameters for 2*dt filter 
 real, parameter    :: kdiff = 1.0E3          ! Coefficient for horizontal diffusion
 real, parameter    :: dt  = 180.0            ! model time step
 integer, parameter :: nhtot = 72             ! number of hours of model integration
 integer, parameter :: npdt = nhtot*3600/dt   ! number of model time steps
 integer, parameter :: nfreq = 24*3600/dt     ! hourly output archiving frequency
 integer, parameter :: ndfi_win = 12          ! time window in hours for digital filter initialisation
 character(len=3)   :: expid='001'            ! experiment identifier
 character(len=8)   :: cdate='15012023'       ! DDMMYYY : date of initial conditions  
 logical            :: lreaduv=.true.         ! logical to use u v at initial time
 logical            :: lsemimp=.true.         ! semi-implicit scheme
 logical            :: linit=.false.          ! DFI initialisation 

 
end module params 
 
module model_vars  
    
 use params   
 
 implicit none
 
 real, dimension (nlon) :: u, u1, u2          ! storage of prognostic variable in physical space
 
 real, dimension (nlon,3) :: u_pdg            ! When is Burgers equation solved in physical space
 
 complex, dimension(-mm:mm,3) :: u_m          ! Prognostic variable in spectral space
 complex, dimension(-mm:mm)   :: u2_m
 
end module model_vars

module spectral_vars

 use params
 
 implicit none
  
! arrays for FFT991 

 integer, dimension(13)         :: ifax
 real, dimension(3*nlon/2+1)    :: trigs 
 real, dimension(nfft*(nlon+2)) :: acoef
 real, dimension(nfft*(nlon+1)) :: work
 
end module spectral_vars
