subroutine hopt(xin,yo)
!------------------------------------------------------------------------------------
! 
! Observation operator - non-linear version
!
! Input  : xin = model state in spectral space 
! Output : yo = simulated observations in physical space
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
 
 complex, dimension(-mm:mm), intent(in)  :: xin  ! model state in spectral space
 real, dimension(nobs),      intent(out) :: yo   ! simulated observations 
 
 real, dimension(nlon)                   :: x_pdg
 integer :: jj, i
 
! 1) Back in physical space 
 
 call fft_i(xin,x_pdg)
 
! 2) Simplified observation operator - projection in obs space 
 
 jj = 1
 do i=1,nlon
   if (mod(i,nsampling) == 0) then
     yo(jj) = x_pdg(i)
     jj = jj + 1      
   endif
 enddo    

 return 

end subroutine hopt
