subroutine fft_i(bm,b)
!-----------------------------------------------------------------------------------------
! Inverse Fourier Transform
!
! Input:
! - bm   : spectral coefficients (2*mm + 1 values)
!
! Output:
! - b    : field in physical space (nlon values)
! 
! It is possible to execute an explicit Fourier transform instead of a FFT (l_notfft)
!
!                                                           J.-F. Mahfouf (05/2025)
!-------------------------------------------------------------------------------------------  
  use fft99_mod
  use spectral_vars
  use mod_vars, only : nlon, mm
  use const, only : j, pi
  
  implicit none
  
  real, dimension(nlon), intent(out)      :: b
  complex, dimension(-mm:mm), intent(in)  :: bm 
  real, dimension(nfft*(nlon+2))               :: zcoef
  integer :: m, i1, j1
  logical :: l_notfft
  
  l_notfft=.false.

  zcoef(:) = 0.0
  do i1 = 1,2*mm+1,2
     zcoef(i1)   = 0.5*real(bm((i1-1)/2) + bm((1-i1)/2))
     zcoef(i1+1) = 0.5*aimag(bm((i1-1)/2) - bm((1-i1)/2))    
  enddo
  call fft991(zcoef,work,trigs,ifax,1,0,nlon,nfft,1)
  b(1:nlon) = zcoef(1:nlon)
  
  if (l_notfft) then
    do j1=1,nlon
      b(j1) = 0.0
      do m=-mm,mm
        b(j1) = b(j1) + real(bm(m)*cexp(2.0*j*pi*float(m)*(float(j1) - 1.0)/float(nlon)))
      enddo
    enddo   
  endif

  return
 
end subroutine fft_i 
