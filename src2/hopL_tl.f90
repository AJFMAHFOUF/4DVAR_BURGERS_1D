subroutine hopL_tl(xin5,xin,yo5,yo,sigmab,c1,dt,npdt,islot)
 
 use mod_vars
 use const
 
 implicit none
 
 complex, dimension(-mm:mm), intent(in)    :: xin5
 complex, dimension(-mm:mm), intent(inout) :: xin  ! initial conditions in spectral space
 real, dimension(nlon), intent(in)         :: sigmab
 integer, intent(in)                       :: npdt, islot
 real, intent(in)                          :: dt, c1
 real, dimension(nobs,nslots), intent(out) :: yo5, yo    ! simulated observations at time slot : islot
 
 real, dimension(nlon)                     :: x5_pdg, x_pdg
 complex, dimension(-mm:mm)                :: xout5, xout
 integer :: jj, i, m

 !do m=-mm,mm
 !  xin(m) = xin(m)*sqrt(c1)/(1.0 + (float(m)*L/a)**2)
 !enddo
 
 !call fft_i(xin,x_pdg)
 
 !do i=1,nlon
 !  x_pdg(i) = x_pdg(i)*sigmab(i)
 !enddo 
   
 !call fft_d(x_pdg,xin)
 
 call hopt_tl(xin5,xin,yo5,yo,dt,npdt,islot)
 
 return
 
end subroutine hopL_tl
