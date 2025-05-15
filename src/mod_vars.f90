module mod_vars

 implicit none

 integer, parameter :: mm = 42 ! 5*5*5         ! maximum wave number
 integer, parameter :: nlon = 3*mm + 2         ! number of longitudes 
 
 real, dimension (nlon) :: u_i, u_f, u, u2     ! storage of prognostic variable in physical space 
 
 complex, dimension(-mm:mm,3)   :: u_m           ! Prognostic variable in spectral space 
 complex, dimension(-mm:mm)     :: u2_m
  
 real, dimension(nlon) :: u5_i, u5_f, u5, u25  ! storage of prognostic variable in physical space (Trajectory)
 real, dimension(nlon) :: u6_i, u6_f
 real, dimension(nlon) :: u7_f, u8_f, u9_f     ! prognostic variables for subsequent prediction
 
 real, dimension(nlon)      :: ut_i, ut_f      ! True model state in physical space
 complex, dimension(-mm:mm) :: ut_m            ! True model state in spectral space
 
 complex, dimension(-mm:mm,3) :: u5_m            ! Prognostic variable in spectral space (Trajectory)
 complex, dimension(-mm:mm)   :: u25_m
 
 real, dimension(nlon) :: sigmab               ! standard deviation of background errors
 
 integer :: nobs, nslots, nfreq  
 integer :: nlength, npdt, nsampling
 real    :: dt, c1, Lb, sigmao  
 
 
end module mod_vars 
 
 
