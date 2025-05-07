subroutine chavarin_ad(xin,xout)
!--------------------------------------------------------------------------------------------
! Change of variable from chi space to model state - adjoint
! [chi]^* = [L^{-1}]^*[dx]^*
! where :
! - [dx]^*     = gradient with respect to the increment in spectral space (input)
! - [chi}^*    = gradient with respect to the control vector (output)
! - [L^{-1}]^* = transpose conjugate of the square-root of B matrix 
!
! Input  : xout = adjoint of model state in spectral space (spectral field)
! Output : xin  = adjoint of model state in chi-space (spectral field) 
! 
! [L^{-1}]^* = C^{1/2}xSx(std)xS^{-1} 
! 
! C = diagonal correlation matrix in spectral space 
! S = direct Fourier transform
! std = diagonal standard deviation background errors in physical space
! S^{-1} = inverse Fourier transform
!
!                                                                  J.-F. Mahfouf (05/2025)
!---------------------------------------------------------------------------------------------
 
 use mod_vars
 use const
 
 implicit none
 
 complex, dimension(-mm:mm), intent(in)    :: xout ! model state in spectral space (adjoint)
 complex, dimension(-mm:mm), intent(out)   :: xin  ! model state in chi-space (adjoint)
 
 real, dimension(nlon)                     :: x_pdg
 complex, dimension(-mm:mm)                :: zvar
 real                                      :: correl
 integer :: i, m
 
 zvar(:) = xout(:)

 call fft_i(zvar,x_pdg)
 
 do i=1,nlon
   x_pdg(i) = x_pdg(i)*sigmab(i)
 enddo 
 
 call fft_d(x_pdg,zvar)

 do m=-mm,mm
   zvar(m) = zvar(m)*sqrt(c1*correl(m))
 enddo   
 
 xin(:) = zvar(:) 
 
 return
 
end subroutine chavarin_ad
