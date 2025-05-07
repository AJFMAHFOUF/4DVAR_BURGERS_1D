subroutine hopt_tl(xin5,xin,yo5,yo)
!------------------------------------------------------------------------------------
! 
! Observation operator - tangent-linear version
!
! Input :  xin5 = model state in spectral space (trajectory)
!          xin  = model state in spectral space (perturbation)
! Output : yo5 = simulated observations in physical space (trajectory)
!          yo  = simulated observations in physical space (perturbation)
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
 
 complex, dimension(-mm:mm), intent(in) :: xin5, xin  ! model state in spectral space
 real, dimension(nobs), intent(out)     :: yo5, yo    ! simulated observations 
 
 real, dimension(nlon)                  :: x5_pdg, x_pdg
 integer :: jj, i
 
! 1) Back in physical space 
 
 call fft_i(xin5,x5_pdg)
 call fft_i(xin,x_pdg)

! 2) Simplified observation operator - projection in obs space 
 
 jj = 1
 do i=1,nlon
   if (mod(i,nsampling) == 0) then
     yo5(jj) = x5_pdg(i)
     yo(jj)  = x_pdg(i)
     jj = jj + 1      
   endif
 enddo    

 return 

end subroutine hopt_tl
