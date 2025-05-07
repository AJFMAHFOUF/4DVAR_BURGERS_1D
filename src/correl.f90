function correl(m)
!----------------------------------------------------------
! Spatial correlation function expressed in spectral space
! as a function of wave number "m"
! => Homogeneous and isotropic in physical space 
!
!                             J.-F. Mahfouf (05/2025)
!----------------------------------------------------------
 use const
 use mod_vars
 implicit none
 real                :: correl
 integer,   intent(in)  :: m
!
 correl = 1.0/(1.0 + (float(m)*Lb/a)**2)**2 ! Second order autoregressive model
! correl = 1.0/(1.0 + (float(m)*Lb/a)**2) ! First order autoregressive model
! correl = exp(-0.5*(float(m)*Lb/a)**2)   ! Gaussian model 
!
end function correl
