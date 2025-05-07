module spectral_vars

 use const
 use mod_vars
 
 implicit none  
 
! arrays for FFT991 

 integer, dimension(13)         :: ifax
 real, dimension(3*nlon/2+1)    :: trigs 
 real, dimension(nfft*(nlon+2)) :: acoef
 real, dimension(nfft*(nlon+1)) :: work
 
end module spectral_vars
