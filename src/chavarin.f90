subroutine chavarin(xin,xout)
!--------------------------------------------------------------------------------
! Change of variable from chi space to model state
! dx = L^{-1}chi
! where :
! - dx = increment (departure from background) in spectral space (output)
! - chi = control vector to be minimised (input)
! - L = square-root of the B^{-1} matrix 
!
! Input  : xin = model state in spectral space (spectral field)
! Output : xout  = model state in chi-space (spectral field) 
! 
! L^{-1} = Sx(std)xS^{-1}xC^{1/2}
! 
! S = direct Fourier transform
! std = diagonal matrix of standard deviation background errors in physical space
! S^{-1} = inverse Fourier transform
! C = diagonal correlation matrix in spectral space 
!
!                                                     J.-F. Mahfouf (05/2025)
!--------------------------------------------------------------------------------
 
 use mod_vars
 use const
 
 implicit none
 
 complex, dimension(-mm:mm), intent(in)    :: xin   ! model state in chi-space (perturbation)
 complex, dimension(-mm:mm), intent(out)   :: xout  ! model state (perturbation)
 
 real, dimension(nlon)                     :: x_pdg
 complex, dimension(-mm:mm)                :: zvar
 real                                      :: correl
 integer :: i, m

 zvar(:) = xin(:)

 do m=-mm,mm
   zvar(m) = zvar(m)*sqrt(c1*correl(m))
 enddo
 
 call fft_i(zvar,x_pdg)
 
 do i=1,nlon
   x_pdg(i) = x_pdg(i)*sigmab(i)
 enddo 
   
 call fft_d(x_pdg,zvar)
 
 xout(:) = zvar(:)
 
 return
 
end subroutine chavarin
