subroutine hopt_ad(xin5,xin,yo5,yo,dt,npdt,islot)
 
 use mod_vars
 
 implicit none
 
 complex, dimension(-mm:mm), intent(in)    :: xin5 ! initial conditions in spectral space (traj)
 complex, dimension(-mm:mm), intent(out)   :: xin  ! output gradient in spectral space
 integer, intent(in)                       :: npdt, islot
 real, intent(in)                          :: dt
 real, dimension(nobs,nslots), intent(out) :: yo5    ! simulated observations at time slot : islot (traj)
 real, dimension(nobs,nslots), intent(in)  :: yo     ! input gradient in physical (obs) space
 
 real, dimension(nlon)                     :: x5_pdg, x_pdg
 complex, dimension(-mm:mm)                :: xout5, xout, zvar_m
 integer :: jj, i
 
!  TRAJECTORY COMPUTATIONS 
!  =======================

! 1) Model integration up from initial time to the observation time 
 
 call simpleburgers(xin5,xout5,dt,npdt) 
 
! 2) Back in physical space 
 
 call fft_i(xout5,x5_pdg)
 
! 3) Simplified observation operator - projection in obs space 
 
 jj = 1
 do i=1,nlon
   if (mod(i,4) == 0) then
     yo5(jj,islot) = x5_pdg(i)
     jj = jj + 1      
   endif
 enddo  
 
!  ADJOINT COMPUTATIONS   
!  ====================

 x_pdg(:) = 0.0
 xout(:) = (0.,0.)
 
! 1) Simplified observation operator - projection in obs space (adjoint) 
 
 jj = 1
 do i=1,nlon
   if (mod(i,4) == 0) then
     x_pdg(i) = x_pdg(i) + yo(jj,islot) 
!     yo(jj,islot) = 0.0 - not used after - therefore not modified by the routine
     jj = jj + 1      
   endif
 enddo   
 
! 2) Back in spectral space (adjoint of inverse FFT)
 
  call fft_d(x_pdg,zvar_m)
  xout(:) = xout(:) + zvar_m(:)*float(nlon)
  x_pdg(:) = 0.   
 
  
! 1) Adjoint model integration from observation time to initial time
 
 call simpleburgers_ad(xin5,xin,xout5,xout,dt,npdt) 

 return 

end subroutine hopt_ad
