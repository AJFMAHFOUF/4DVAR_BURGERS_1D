module const

 implicit none

 complex, parameter :: j = (0,1)              ! square root of -1 
 real, parameter    :: a = 6371.22E3          ! Earth radius
 real, parameter    :: pi = acos(-1.0)        ! Pi constant
 real, parameter    :: g = 9.80616            ! Earth gravitational acceleration
 real, parameter    :: omega = 2.0*pi/86164.1 ! Earth angular speed (stellar day)
 real, parameter    :: nu = 0.05, wk = 0.53   ! tunable parameters for 2*dt filter 
 real, parameter    :: kdiff = 1.0E3          ! Coefficient for horizontal diffusion 
 integer, parameter :: nfft = 1               ! number of FFT to be done
 real, parameter    :: L = 500.0E3            ! Correlation length for background errors
 
end module const
