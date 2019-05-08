program plot_phases
use prolates

implicit double precision (a-h,o-z)
double precision, allocatable :: coefsps(:), coefsqs(:)
double precision, allocatable :: xs(:),ys(:)
double complex                :: clambda

nn = 1000
allocate(xs(nn),ys(nn))


c   = 500.0d0
n   = 400
call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)
call prolates_integer_coefficients(c,n,ncoefs,coefsps,coefsqs,c1,c2)

aa = sum(coefsps)

call prind("aa    = ",aa)
call prind("chi   = ",chi)

call prind("c1    = ",c1)
call prind("c2    = ",c2)

!
!  Plot the first phase function
!
do i=1,nn
x = -1.0d0 + (i+0.0d0)/(nn+0.0d0)  * 2
call prolates_integer_evalder(n,ncoefs,coefsps,coefsqs,x,valps,derps,valqs,derqs)

xs(i) = x
ys(i) = 1/ ( (c1*valps)**2 + (c2*valqs)**2)
end do

call plot_function("phaseder1.pdf","","",nn,xs,ys)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c   = 5000.0d0
n   = 4000

call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)
aa = sum(coefsps)

call prind("aa    = ",aa)
call prind("chi   = ",chi)


call prolates_integer_coefficients(c,n,ncoefs,coefsps,coefsqs,c1,c2)
call prind("c1    = ",c1)
call prind("c2    = ",c2)

!
!  Plot the second phase function
!
do i=1,nn
x = -1.0d0 + (i+0.0d0)/(nn+0.0d0)  * 2
call prolates_integer_evalder(n,ncoefs,coefsps,coefsqs,x,valps,derps,valqs,derqs)

xs(i) = x
ys(i) = 1/ ( (c1*valps)**2 + (c2*valqs)**2)
end do

call plot_function("phaseder2.pdf","","",nn,xs,ys)


end program
