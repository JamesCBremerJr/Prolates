program test_prolates
use prolates

implicit double precision (a-h,o-z)
double precision, allocatable :: coefsps(:), coefsqs(:)
double precision, allocatable :: xs(:),ys(:)
double complex                :: clambda




c   = 5000.0d0
n   = 4000

call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)
aa = sum(coefsps)

call prind("aa    = ",aa)
call prind("chi   = ",chi)


call prolates_integer_coefficients(c,n,ncoefs,coefsps,coefsqs,c1,c2)
call prind("c1    = ",c1)
call prind("c2    = ",c2)

! !
! !  Plot the phase functions
! !
! nn = 1000
! allocate(xs(nn),ys(nn))
! do i=1,nn
! x = -1.0d0 + (i+0.0d0)/(nn+0.0d0)  * 2
! call prolates_integer_evalder(n,ncoefs,coefsps,coefsqs,x,valps,derps,valqs,derqs)

! xs(i) = x
! ys(i) = 1/ ( (c1*valps)**2 + (c2*valqs)**2)
! end do


! call plot_function("plot.pdf","","",nn,xs,ys)



!
!  Check the Wronskian formula ... note that ps_n(x,c) and qs_n(x,c) are
!  a very badly conditioned basis when c is large and n is not
!

x = 0.5d0
call prolates_integer_evalder(n,ncoefs,coefsps,coefsqs,x,valps,derps,valqs,derqs)

W  = valps*derqs-derps*valqs
W0 = aa**2/(1-x**2)

call prind("Wronskian as computed = ",W)
call prind("Wronskian via formula = ",W0)

! errrel = abs(W-W0)/abs(W0)
! call prind("relative error in Wronskian = ",errrel)


!
!  Check the accuracy of the asymptotic estimate of chi
!

call prolates_integer_bound(c,n,chi1,chi2)
chi0   = (chi1+chi2)/2
errrel = (chi-chi0)/chi

call prind("chi0  = ",chi0)
call prin2("relative error in asymptotic estimate of chi = ",errrel)


end program
