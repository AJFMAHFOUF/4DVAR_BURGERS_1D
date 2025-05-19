subroutine create_obs(xin5,yo5,l_obspert)
!---------------------------------------------------------------------------------
! Create a series of pseudo-observations from a reference state 
!
! Call of :
! - burgers = propagates the initial state xin5 to the observation time
! - hopt    = converts the model state into the observation space 
! - gasdev  = computes random noise following a Gaussian distribution
!
! Inputs  : xin5 = model state at the beginning of the assimilation window
!           l_obspert = true for noisy observations
! Outputs : yo5 = time series of pseudo-observations for nslots time slots
!
! The number of observations "nobs" per time slots
! and the number of time slots "nslots" are defined in the subroutine init.f90
!
!                             J.-F. Mahfouf (05/2025)
!---------------------------------------------------------------------------------
 
 use mod_vars
 
 implicit none
 
 complex, dimension(-mm:mm), intent(in)      :: xin5
 logical, intent(in)                         :: l_obspert
 real, dimension(nobs,0:nslots), intent(out) :: yo5
 
 integer :: islot, ii
 
 complex, dimension(-mm:mm)   :: xout5, zvar
 real, dimension(nobs,0:nslots) :: eta_o
 
 zvar(:) = xin5(:)
 
 if (l_obspert) then
   do islot=0,nslots
     do ii=1,nobs
       call gasdev(eta_o(ii,islot))
       eta_o(ii,islot) = eta_o(ii,islot)*sigmao 
     enddo
   enddo
 endif
 
 if (l_obspert) read (211,*) eta_o

! The above random array can be saved for reproductibility of the results 
 
 call hopt(xin5,yo5(:,0))
 if (l_obspert) yo5(:,0) = yo5(:,0) + eta_o(:,0)
 
 do islot=1,nslots
   call burgers(zvar,xout5,npdt)
   call hopt(xout5,yo5(:,islot))
   if (l_obspert) yo5(:,islot) = yo5(:,islot) + eta_o(:,islot)
   zvar(:) = xout5(:)
 enddo 
 
 close(211)

return
end subroutine create_obs
