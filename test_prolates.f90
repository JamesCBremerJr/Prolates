program test_prolates
use prolates

implicit double precision (a-h,o-z)
double precision, allocatable :: coefsps(:)
double complex                :: zlam

double precision, allocatable   ::  psi(:), xs(:), ys(:), zs(:,:)



eps0    = epsilon(0.0d0)
eps     = 2.0d0**(-53.0d0)
pi      = acos(-1.0d0)
gamma   =   50.0d0
n       =   17
x       =   0.5d0
inormal = 1

call prolates_dimension(gamma,eps,ndim)

call prin2("gamma     = ",gamma)
call prini("dimension = ",ndim)

! reference values
chi0   = 1578.36224927541989403207235632293111d0
valps0 = 0.0434046611865871801601265621207536732d0
derps0 = 3.35154420055253986102883382869371796d0

call elapsed(t1)
call prolates_integer(inormal,gamma,n,chi,n1,n2,coefsps)
call elapsed(t2)

call prin2("prolates_integer time = ",t2-t1)

errrel = abs(chi-chi0)/abs(chi0)
call prin2("error in eigenvalue = ",errrel)

call  prolates_integer_evalder(n,n1,n2,coefsps,x,valps,derps)
errrel = abs(valps-valps0)/abs(valps0)
call prin2("error in value = ",errrel)
errrel = abs(derps-derps0)/abs(derps0)
call prin2("error in derivative = ",errrel)


end program
