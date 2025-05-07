subroutine hopt_ad(xin5,xin,yo5,yo)
!------------------------------------------------------------------------------------
! 
! Observation operator - adjoint version
!
! Input :  xin5 = model state in spectral space (trajectory)
!          yo   = gradient of Jo with respect to observation in physical space
! Output : yo5  = simulated observations in physical space (trajectory)
!          xin  = gradient of Jo with respect to model state in spectral space
!
! The model state and the simulated observation are available for the same time
! 
! Note: since the model state is given in spectral space 
!       an inverse Fourier transform needs to be applied before 
!       calling the "actual" observation operator
!
! The observation operator is a simple spatial sampling from the collocation grid       
!
!                                                           J.-F. Mahfouf (05/2025)
!-------------------------------------------------------------------------------------
 
 use mod_vars
 
 implicit none
 
 complex, dimension(-mm:mm), intent(in)  :: xin5 ! model state in spectral space (traj)
 complex, dimension(-mm:mm), intent(out) :: xin  ! output gradient in spectral (model) space
 real, dimension(nobs), intent(out)      :: yo5  ! simulated observations 
 real, dimension(nobs), intent(in)       :: yo   ! input gradient in physical (obs) space
 
 real, dimension(nlon)                   :: x5_pdg, x_pdg
 complex, dimension(-mm:mm)              :: zvar_m 
 integer :: jj, i
 
!  TRAJECTORY COMPUTATIONS 
!  =======================
 
! 1) Back in physical space 
 
 call fft_i(xin5,x5_pdg)
 
! 2) Simplified observation operator - projection in obs space 
 
 jj = 1
 do i=1,nlon
   if (mod(i,nsampling) == 0) then
     yo5(jj) = x5_pdg(i)
     jj = jj + 1      
   endif
 enddo  
 
!  ADJOINT COMPUTATIONS   
!  ====================

 x_pdg(:) = 0.0
 xin(:) = (0.0,0.0)
 
! 1) Simplified observation operator - projection in obs space (adjoint) 
 
 jj = 1
 do i=1,nlon
   if (mod(i,nsampling) == 0) then
     x_pdg(i) = x_pdg(i) + yo(jj) 
!     yo(jj) = 0.0 - not used after - therefore not modified by the routine
     jj = jj + 1      
   endif
 enddo   
 
! 2) Back in spectral space (adjoint of inverse FFT)
 
  call fft_d(x_pdg,zvar_m)
  xin(:) = xin(:) + zvar_m(:)*float(nlon)
  x_pdg(:) = 0.   
  
 return 

end subroutine hopt_ad
