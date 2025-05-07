subroutine hopt_tl(xin5,xin,yo5,yo,dt,npdt,islot)
 
 use mod_vars
 
 implicit none
 
 complex, dimension(-mm:mm), intent(in)    :: xin5, xin  ! initial conditions in spectral space
 integer, intent(in)                       :: npdt, islot
 real, intent(in)                          :: dt
 real, dimension(nobs,nslots), intent(out) :: yo5, yo    ! simulated observations at time slot : islot
 
 real, dimension(nlon)                     :: x5_pdg, x_pdg
 complex, dimension(-mm:mm)                :: xout5, xout
 integer :: jj, i

! 1) Model integration up from initial time to the observation time 

! print *, 'xin5 in TL',xin5(10),islot,dt
 
 call simpleburgers_tl(xin5,xin,xout5,xout,dt,npdt) 
 
! print *,'xout5 in TL',xout5(10)
 
! 2) Back in physical space 
 
 call fft_i(xout5,x5_pdg)
 call fft_i(xout,x_pdg)

! 3) Simplified observation operator - projection in obs space 
 
 jj = 1
 do i=1,nlon
   if (mod(i,4) == 0) then
     yo5(jj,islot) = x5_pdg(i)
     yo(jj,islot)  = x_pdg(i)
     jj = jj + 1      
   endif
 enddo    

 return 

end subroutine hopt_tl
