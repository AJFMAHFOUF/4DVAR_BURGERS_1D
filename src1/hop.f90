subroutine hop(x,yo,islot)
 
 use mod_vars
 
 implicit none
 
 real, dimension(nlon), intent(in)         :: x
 real, dimension(nobs,nslots), intent(out) :: yo
 integer :: jj, i, islot
 
 jj = 1
 do i=1,nlon
   if (mod(i,4) == 0) then
     yo(jj,islot) =  x(i)
     jj = jj + 1      
   endif
 enddo    

 return 

end subroutine hop
