program plot_phases
use prolates

implicit double precision (a-h,o-z)
double precision, allocatable :: coefsps(:), coefsqs(:)
double precision, allocatable :: xs(:),ys(:)
double complex                :: clambda

c   = 100
nn  = c
allocate(xs(nn),ys(nn))


do n=1,nn
call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)
call prolates_integer_coefficients(c,n,ncoefs,coefsps,coefsqs,c1,c2)

aa = sum(coefsps)

xs(n) = n
ys(n) = aa

print *,c,n,aa
end do


call plot_function("Aplot.pdf","","",nn,xs,ys)

stop


end program
