subroutine hopL_ad(xin5,xin,yo5,yo,sigmab,c1,dt,npdt,islot)
 
 use mod_vars
 use const
 
 implicit none
 
 complex, dimension(-mm:mm), intent(in)    :: xin5
 complex, dimension(-mm:mm), intent(out)   :: xin  ! initial conditions in spectral space
 real, dimension(nlon), intent(in)         :: sigmab
 integer, intent(in)                       :: npdt, islot
 real, intent(in)                          :: dt, c1
 real, dimension(nobs,nslots), intent(out) :: yo5
 real, dimension(nobs,nslots), intent(in ) :: yo      ! simulated observations at time slot : islot
 
 real, dimension(nlon)                     :: x5_pdg, x_pdg
 complex, dimension(-mm:mm)                :: xout5, xout, zvar
 integer :: jj, i, m
 
 call hopt_ad(xin5,zvar,yo5,yo,dt,npdt,islot)

 call fft_i(zvar,x_pdg)
 
 do i=1,nlon
   x_pdg(i) = x_pdg(i)*sigmab(i)
 enddo 
 
 call fft_d(x_pdg,zvar)

 do m=-mm,mm
    zvar(m) = zvar(m)*sqrt(c1)/(1.0 + (float(m)*L/a)**2)
 enddo   
 

 
 xin(:) = zvar(:) 
 
 return
 
end subroutine hopL_ad
