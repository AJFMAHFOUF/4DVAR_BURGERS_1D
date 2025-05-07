subroutine hopt(xin,yo,dt,npdt,islot)
 
 use mod_vars
 
 implicit none
 
 complex, dimension(-mm:mm), intent(in)    :: xin  ! initial conditions in spectral space
 integer, intent(in)                       :: npdt, islot
 real, intent(in)                          :: dt
 real, dimension(nobs,nslots), intent(out) :: yo   ! simulated observations at time slot : islot
 
 real, dimension(nlon)                     :: x_pdg
 complex, dimension(-mm:mm)                :: xout
 integer :: jj, i

! 1) Model integration up from initial time to the observation time 
 
 call simpleburgers(xin,xout,dt,npdt)
 
! 2) Back in physical space 
 
 call fft_i(xout,x_pdg)
 
! 3) Simplified observation operator - projection in obs space 
 
 jj = 1
 do i=1,nlon
   if (mod(i,4) == 0) then
     yo(jj,islot) = x_pdg(i)
     jj = jj + 1      
   endif
 enddo    

 return 

end subroutine hopt
