subroutine fft_d(b,bm)
!-----------------------------------------------------------------------------------------
! Direct Fourier Transform
!
! Input:
! - b    : field in physical space (nlon values)
!
! Output:
! - bm   : spectral coefficients (2*mm + 1 values)
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
  
  real, dimension(nlon), intent(in)       :: b
  complex, dimension(-mm:mm), intent(out) :: bm 
  real, dimension(nlon+2)                 :: zcoef
  integer :: m, i1, j1
  logical :: l_notfft
  
  l_notfft=.false.
  
  zcoef(:) = 0.0
  zcoef(1:nlon) = b(1:nlon)
  call fft991(zcoef,work,trigs,ifax,1,0,nlon,nfft,-1)
  do i1 = 1,2*mm+1,2
    bm((i1-1)/2) = zcoef(i1) + j*zcoef(i1+1)
    bm((1-i1)/2) = zcoef(i1) - j*zcoef(i1+1)  
  enddo
  
  if (l_notfft) then
    do m=0,mm
      bm(m) = 0.0
      do j1=1,nlon
        bm(m) =  bm(m) + b(j1)*cexp(-2.0*j*float(m)*pi*(float(j1)-1.0)/float(nlon)) 
      enddo
      bm(m) = bm(m)/float(nlon)
      bm(-m) = conjg(bm(m))
    enddo
  endif
   
return
 
end subroutine fft_d   
        
  
