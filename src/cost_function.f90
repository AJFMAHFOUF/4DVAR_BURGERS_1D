subroutine cost_function(chi,d0,xin5,xj,xjo,xjb)
!--------------------------------------------------------------------------------------------------
! Estimation of the cost-function J (incremental formulation - linearized obs operator
!
! Call of:
!
! - burgers_tl    : Dynamical model for the time propagation of increments
! - hopt_tl       : Conversion of spectral increments into observation space
! 
! Inputs:
! - chi : control vector in transformed space (chi = Ldx where L is the square-root of B^{-1})
! - d0 : innovation vector = difference between simulated and real observations for each timeslot
! - xin5 : initial conditions in spectral space at the beginning of the assimilation window
!
! Outputs:
! - xj : total cost-function
! - xjo : cost-function (observation term)
! - xjb : cost-function (background term) 
!
!  The two called routines are necessary to compute the vector H(dx) where 
!  H is the linearized observation operator (including an inverse spectral transform) 
!  and dx the increment in spectral space 
!
!
!                                                           J.-F. Mahfouf (05/2025)
!----------------------------------------------------------------------------------------------------  
 use mod_vars
 
 implicit none
 
 complex, dimension(-mm:mm), intent(in)     :: chi, xin5
 real, dimension(nobs),      intent(in)     :: d0
 real, intent(out)                          :: xj, xjo, xjb
 
 complex, dimension(-mm:mm) :: zvar, zvar5, xout, xout5 
 real, dimension(nobs)      :: yo5, hdx
 integer :: m, ii

 xjb = 0.0
 do m=-mm,mm
   xjb = xjb + conjg(chi(m))*chi(m)
 enddo
 xjb = 0.5*xjb
 
! 2) Observation term
 
 zvar5(:) = xin5(:)
 
 xjo = 0.0
 
 call chavarin(chi,zvar)
 
 call hopt_tl(zvar5,zvar,yo5,hdx) 
 do ii=1,nobs
   xjo = xjo + ((hdx(ii) - d0(ii))/sigmao)**2
 enddo 
 
 xjo = 0.5*xjo
 
 xj = xjo + xjb 
 
return

end subroutine cost_function   
