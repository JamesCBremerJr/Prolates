!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for constructing and evaluating expansions of the angular
!  prolate spheroidal wave functions of the first kind of order 0 and integer
!  characteristic values in Legendre polynomials.
!
!  These prolate functions are the solutions of the reduced spheroidal wave equation
!
!    (1-x^2) * y''(x) - 2x y'(x) + ( \chi - \gamma^2 x^2 ) y(x) = 0                       (1)
!
!  on the interval (-1,1) which satisfy the boundary conditions
!
!    \lim   y'(x) (1-x^2) = 0.                                                            (2)
!    x \to \pm 1
!
!  The equation (1) together with (2) is a singular self-adjoint Sturm-Liouville
!  problem.  Consequently, for fixed values of \gamma, there exist a countable number
!  of values of chi for which solutions of (1) satisfying (2) exist.  These values
!  chi are known as Sturm-Liouville eigenvalues, and we denote them by
!
!    \chi_0 (gamma^2) < \chi_1 (gamma^2) < \chi_2 (gamma^2) < ...
!
!  For each such eigenvalue, there is a corresponding one-dimensional eigenspace
!  of functions which satisfy (2) and a choice of normalization is required to
!  determine the angular spheroidal wave function of the first kind
!       
!    ps   (x; gamma^2)                                                                    (3)
!      n 
!
!  uniquely.
!
!  At the present time, this code provide two choices of normalization scheme:
!
!  (1) "L^2 normalization" in which (3) is chosen so that
!
!          1
!      \int  ps_n(x)^2  dx  = 1 
!         -1
!
!      and the sign of ps_n(x) is indeterminate.
!
!  (2) "Legendre normalization" in which (3) is chosen so that the value of (3) at 0 
!      agrees with the value of the Legendre polynomial of degree n at 0 if n is even
!      and the derivative of (3) at 0 agrees with the derivative of the Legendre 
!      polynomial of degree n at 0 n the case n is odd.
!
!  The function (3) admits an expansion of the form
!
!                       \infty      ~ 
!    ps (x,\gamma^2) =   \sum    a  P           (x),                                      (4)
!      n                 j = N    j   n + 2*j
!
!  where N = -Floor(n/2) and
!
!    ~ 
!    P  (x)  
!      n 
!
!  denotes the L^2(-1,1) normalized Legendre polynomial of degree n.  This code 
!  operates by constructing a truncated version (4) via the Osipov-Rokhlin method.
!
!  The following subroutines should be regarded as publicly callable:
!
!    prolates_integer - compute the the Sturm-Liouville eigenvalue chi_n and the 
!      coefficients in the associated Legendre expansions (4) of the function
!      ps_n(x; \gamma^2) using the Rokhlin-Osipov method
!
!    prolates_integer_eval - evaluate an expansion of the form (4) at a specified
!      point in (-1,1); the Legendre polynomials are calculated using the three-term 
!      recurrence relations
!
!    prolates_integer_evalder - evaluate an expansion of the form (4) and its derivative
!      at the a specified point in the interval (-1,1)
!
!    prolates_integer_lambda - the function ps_n(x; gamma^2) is also an eigenfunction of
!      the operator
!
!                        1
!        T[f](x) = \int      exp(i \gamma x t) f(t) dt.                                   (5)
!                       -1                       
!
!      This subroutine returns the eigenvalue of T corresponding to (4).
!
!    prolates_dimension - given a value of gamma and eps, this routine returns the
!      least integer n such that \lambda_n < eps, where \lambda_1,\lambda_2,...
!      are the eigenvalues of the integral operator T in (5)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module     prolates
use utils
contains


subroutine prolates_integer(inormal,gamma,n,chi,n1,n2,coefsps)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out)  :: coefsps(:)
!
!  Calculate the Sturm-Liouville eigenvalue chi_n(\gamma^2) and construct the 
!  coefficients in the expansion (4).
!                                                                                         
!  Input parameters:
!    inormal - an integer parameter which specifies the normalization scheme to use;
!
!       inormal = 0  means use "L^2 normalization"
!       inormal = 1  means use "Legendre normalization"
!
!    gamma - the "bandlimit" of the PSWF function to construct, which must be
!      a nonnegative real number
!    n - the characteristic exponent, which must be a nonnegative integer
!
!  Output parameters:
!    chi - the Sturm-Liouville eigenvalue
!    (n1, n2)  - the lower and upper indices in the sum (6)
!    coefsps - an array containing the coefficients in the expansion (6)
!
data pi / 3.14159265358979323846264338327950288d0 / 
double precision, allocatable :: as(:),bs(:),coefsps0(:)

eps    = epsilon(0.0d0)

!  Estimate the number of terms in the expansion
!mm = 50 + 2*(n) + 2*floor(sqrt(gamma*(n)))
mm = 50 + 2/pi*n + floor(sqrt(gamma*(n)))

! this is the quantity suggested in Osipov- Rokhlin-Xiao
!mm = 1000 + (1.1)*n + floor(sqrt(gamma))

1000 continue

!  Determine the limits of the index in the expansion
nn1 = max(-floor(n/2.0d0),-mm/2) 
nn2 = mm+nn1-1

!  Construct the symmetric tridiagonal matrix; the diagonal elements
!  are stored in the array as and the off-diagonal entrys in the array
!  bs
allocate(coefsps0(nn1:nn2),as(nn1:nn2),bs(nn1+1:nn2))

!  Form the appropriate symmetric tridiagonal matrix 

dn = n
dm = 0

do k=nn1,nn2
dk    = k
as(k) = (0.2d1*dk+dn)*(1+0.2d1*dk+dn)+gamma**2-(0.2d1*(-1+dm**2+  &
        (2*dk+dn)*(1+2*dk+dn))*gamma**2)/((-1+4*dk+2*dn)*(3+4*dk+  &
        2*dn))
end do

do k=nn1+1,nn2
dk    = k
bs(k) = (sqrt(((-0.1d1+0.2d1*dk-0.1d1*dm+dn)*(0.2d1*dk-0.1d1*dm+  &
        dn)*(-0.1d1+0.2d1*dk+dm+dn)*(0.2d1*dk+dm+dn))/((-0.3d1+  &
        0.4d1*dk+0.2d1*dn)*(0.1d1+0.4d1*dk+0.2d1*dn)))*gamma**2)/(-1  &
        +0.4d1*dk+0.2d1*dn)
end do


! Compute the specified eigenvalue using the Rokhlin-Osipov technique
nn  = nn2-nn1+1
idx = -nn1+1

call prolates_trieigen(nn,as,bs,idx,chi,coefsps0)

! Truncate the coefficient expansion so as to eliminate coefficients of 
! extremely small magnitude


dd  = abs(coefsps0(0))
eps = epsilon(0.0d0)**2

do i=nn1,0
if ( abs(coefsps0(i)) .ge. eps*dd) then
mm1 = i
exit
endif
end do

do i=nn2,0,-1
if ( abs(coefsps0(i)) .gt. eps*dd) then
mm2 = i
exit
endif
end do


!  If the number of coefficients is insufficient, repeat the process


if ( mm2 .eq. nn2) then
mm = mm*2
deallocate(coefsps0,as,bs)
goto 1000
endif

n1 = mm1
n2 = mm2


allocate(coefsps(n1:n2))
coefsps = coefsps0(n1:n2)

! this fixed the sign
dmax = 0
idx  = 0
do i=n1,n2
if (abs(coefsps(i)) .gt. dmax) then
dmax = abs(coefsps(i))
idx  = i
endif
end do

! normalize the coefficiens

if (inormal .eq. 1) then
call prolates_integer_evalder(n,n1,n2,coefsps,0.0d0,valps0,derps0)
call prolates_legen0(n,val0,der0)

if (mod(n,2) == 0) then
coefsps = coefsps * val0/valps0
else
coefsps = coefsps * der0/derps0
endif

else 

dnorm   = sum(coefsps**2)
coefsps = coefsps / sqrt(dnorm)


endif

end subroutine


subroutine prolates_integer2(inormal,gamma,n,chi,n1,n2,coefsps,mm)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out)  :: coefsps(:)
!
!  Calculate the Sturm-Liouville eigenvalue chi_n(\gamma^2) and construct the 
!  coefficients in the expansion (4).
!                                                                                         
!  Input parameters:
!    inormal - an integer parameter which specifies the normalization scheme to use;
!
!       inormal = 0  means use "L^2 normalization"
!       inormal = 1  means use "Legendre normalization"
!
!    gamma - the "bandlimit" of the PSWF function to construct, which must be
!      a nonnegative real number
!    n - the characteristic exponent, which must be a nonnegative integer
!
!  Output parameters:
!    chi - the Sturm-Liouville eigenvalue
!    (n1, n2)  - the lower and upper indices in the sum (6)
!    coefsps - an array containing the coefficients in the expansion (6)
!
data pi / 3.14159265358979323846264338327950288d0 / 
double precision, allocatable :: as(:),bs(:),coefsps0(:)

eps    = epsilon(0.0d0)

!  Estimate the number of terms in the expansion
!mm = 50 + 2*(n) + 2*floor(sqrt(gamma*(n)))
mm = 50 + 2/pi*n + floor(sqrt(gamma*(n)))

! this is the quantity suggested in Osipov- Rokhlin-Xiao
!mm = 1000 + (1.1)*n + floor(sqrt(gamma))

1000 continue

!  Determine the limits of the index in the expansion
nn1 = max(-floor(n/2.0d0),-mm/2) 
nn2 = mm+nn1-1

!  Construct the symmetric tridiagonal matrix; the diagonal elements
!  are stored in the array as and the off-diagonal entrys in the array
!  bs
allocate(coefsps0(nn1:nn2),as(nn1:nn2),bs(nn1+1:nn2))

!  Form the appropriate symmetric tridiagonal matrix 

dn = n
dm = 0

do k=nn1,nn2
dk    = k
as(k) = (0.2d1*dk+dn)*(1+0.2d1*dk+dn)+gamma**2-(0.2d1*(-1+dm**2+  &
        (2*dk+dn)*(1+2*dk+dn))*gamma**2)/((-1+4*dk+2*dn)*(3+4*dk+  &
        2*dn))
end do

do k=nn1+1,nn2
dk    = k
bs(k) = (sqrt(((-0.1d1+0.2d1*dk-0.1d1*dm+dn)*(0.2d1*dk-0.1d1*dm+  &
        dn)*(-0.1d1+0.2d1*dk+dm+dn)*(0.2d1*dk+dm+dn))/((-0.3d1+  &
        0.4d1*dk+0.2d1*dn)*(0.1d1+0.4d1*dk+0.2d1*dn)))*gamma**2)/(-1  &
        +0.4d1*dk+0.2d1*dn)
end do


! Compute the specified eigenvalue using the Rokhlin-Osipov technique
nn  = nn2-nn1+1
idx = -nn1+1

call prolates_trieigen(nn,as,bs,idx,chi,coefsps0)

! Truncate the coefficient expansion so as to eliminate coefficients of 
! extremely small magnitude


dd  = abs(coefsps0(0))
eps = epsilon(0.0d0)**2

do i=nn1,0
if ( abs(coefsps0(i)) .ge. eps*dd) then
mm1 = i
exit
endif
end do

do i=nn2,0,-1
if ( abs(coefsps0(i)) .gt. eps*dd) then
mm2 = i
exit
endif
end do


!  If the number of coefficients is insufficient, repeat the process


if ( mm2 .eq. nn2) then
mm = mm*2
deallocate(coefsps0,as,bs)
goto 1000
endif

n1 = mm1
n2 = mm2


allocate(coefsps(n1:n2))
coefsps = coefsps0(n1:n2)

! this fixed the sign
dmax = 0
idx  = 0
do i=n1,n2
if (abs(coefsps(i)) .gt. dmax) then
dmax = abs(coefsps(i))
idx  = i
endif
end do

! normalize the coefficiens

if (inormal .eq. 1) then
call prolates_integer_evalder(n,n1,n2,coefsps,0.0d0,valps0,derps0)
call prolates_legen0(n,val0,der0)

if (mod(n,2) == 0) then
coefsps = coefsps * val0/valps0
else
coefsps = coefsps * der0/derps0
endif

else 

dnorm   = sum(coefsps**2)
coefsps = coefsps / sqrt(dnorm)


endif

end subroutine


subroutine prolates_integer_eval(n,n1,n2,coefsps,x,valps)
implicit double precision (a-h,o-z)
double precision :: coefsps(n1:n2)
!
!  Evaluate the expansion (4) at a specified point x in the interval (-1,1).
!
!  Input parameters:
!    n - the characteristic exponent of the PSWF
!    n1 - the lower index in the sum in (4) as returned by prolates_integer
!    n2 - the upper index in the sum in (4) as returned by prolates_integer
!    coefsps - the coefficients in the expansion (4) as returned by prolates_integer
!    x - the point in (-1,1) at which to evaluate the expansion
!
!  Output parameters:
!    valps - the value of ps_n(x;\gamma^2)
!

double precision, allocatable :: pols(:)

! Determine how many Legendre polynomials we need to evaluate
m = n+2*n2
allocate(pols(0:m) )
call prolates_leges(m+1,x,pols)

! Sum the expansion
valps = 0
idx   = n+2*n1
do i=n1,n2
valps = valps + pols(idx)*coefsps(i)
idx   = idx + 2
end do

end subroutine


subroutine prolates_integer_evalder(n,n1,n2,coefsps,x,valps,derps)
implicit double precision (a-h,o-z)
double precision :: coefsps(n1:n2)
!
!  Evaluate the expansion (4) and its derivative at a specified point x in the
!  interval (-1,1).
!
!  Input parameters:
!    n - the characteristic exponent of the PSWF
!    n1 - the lower index in the sum in (4) as returned by prolates_integer
!    n2 - the upper index in the sum in (4) as returned by prolates_integer
!    coefsps - the coefficients in the expansion (4) as returned by prolates_integer
!    x - the point in (-1,1) at which to evaluate the expansion
!
!  Output parameters:
!    valps - the value of ps_n(x;\gamma^2)
!    derps - the derivative of ps_n(x;\gamma^2) at the specified point
!

double precision, allocatable :: pols(:), ders(:)
!double precision, allocatable :: alfs(:)

! Determine how many Legendre polynomials we need to evaluate
m = n+2*n2
allocate(pols(0:m), ders(0:m) )
call prolates_legeders(m+1,x,pols,ders)

! Sum the expansion
valps = 0
derps = 0
idx   = n+2*n1
do i=n1,n2
valps = valps + pols(idx)*coefsps(i)
derps = derps + ders(idx)*coefsps(i)
idx   = idx + 2
end do

end subroutine



subroutine prolates_integer_lambda(gamma,n,n1,n2,coefsps,zlam)
implicit double precision (a-h,o-z)
double precision  :: coefsps(n1:n2)
double complex    :: zlam
!
!  Calculate the eigenvalue of the integral operator
!
!                        1
!        T[f](x) = \int      exp(i \gamma x t) f(t) dt
!                       -1
!
!  corresponding to ps_n(x).
!                                                                                         
!  Input parameters:
!    gamma - the "bandlimit" of the PSWF function 
!    (n1, n2)  - the lower and upper indices in the sum (6)
!    coefsps - an array containing the coefficients in the expansion (6)
!
!  Output parameters:
!    zlam - the eigenvalue
!

x = 0.0d0
call  prolates_integer_evalder(n,n1,n2,coefsps,x,valps,derps)

dd = n+2*n1
if (mod(n,2) == 0) then
val=coefsps(n1)*sqrt(2.0d0/(2*dd+1.0d0))
zlam = val/valps
else
der = gamma*coefsps(n1)*sqrt(2.0d0/(2*dd+1.0d0))
zlam = der/derps * (0.0d0,1.0d0)
endif

end subroutine

subroutine prolates_dimension(gamma,eps,ndim)
implicit double precision (a-h,o-z)
!
!  Return the least integer n such that \lambda_n < eps, where \lambda_1,\lambda_2,...
!  are the eigenvalues of the integral operator (5).  This is the numerical dimension of
!  the space of functions with bandlimit gamma to precision eps.
!
!  Input parameters:
!    gamma - the bandlimit
!    eps - the precision
!
!  Output parameters:
!    ndim - the dimension of the space
!

double precision, allocatable :: coefsps(:)
double complex                :: zlam, zlam0

pi   = acos(-1.0d0)
m0   = 2/pi*gamma
ndim = -1


dd = abs(zlam0)

do m=m0,10000000
call prolates_integer(inormal,gamma,m,chi,n1,n2,coefsps)
call prolates_integer_lambda(gamma,m,n1,n2,coefsps,zlam)
zlam = zlam*sqrt(gamma/(2*pi))
dd   = abs(zlam)
if (dd .lt. eps) then
ndim = m
return

endif
end do


end subroutine


subroutine prolates_leges(n,x,pols)
implicit double precision (a-h,o-z)
double precision           :: pols(:)
!
!  Return the values of the L^2 normalized Legendre polynomials of the
!  first kind of degrees 0 through n-1 at a specified point.
!
!  Input parameters:
!     n - the number of polynomials to evaluate
!     x - the point at which to evaluate them
!
!  Output parameters:
!     pols - the array of length n containing the values of the first n 
!      Legendre polynomials
!

if (x == 1.0d0) then
do i=1,n
pols(i) = 1.0d0
end do
goto 1000
endif

if (x == -1.0d0) then
dsign = 1.0d0
do i=1,n
pols(i) = dsign
dsign   = -dsign
end do
goto 1000
endif


pols(1) = 1.0d0
if (n == 1) goto 1000
pols(2) = x
if (n == 2) goto 1000

do j=2,n-1
pols(j+1) = ((2*j-1)*x*pols(j)-(j-1)*pols(j-1))/j
end do

!
!  Normalize the polynomials
!
1000 continue

do j=1,n
dd      = sqrt(j - 0.5d0)
pols(j) = pols(j) * dd
end do

end subroutine


subroutine prolates_legeders(n,x,pols,ders)
implicit double precision (a-h,o-z)
integer                       :: n
double precision              :: x
double precision, intent(out) :: pols(:),ders(:)
!
!  Evaluate the L^2 normalized Legendre polynomials of degree 0 through n-1 at the
!  point x using the standard 3-term recurrence relation.  Return the values
!  of their derivative at the point x as well.
!
!  Input parameters:
!    n - an integer specifying the number of polynomials to be evaluated
!    x - the point at which to evaluate the polynomials and their derivatives 
!
!  Output parameters:
!    pols - the ith entry of this user-supplied array of length n
!      will contain the value of normalized Legendre polynomial of degree
!      i-1 at x
!    ders - the ith entry of this user-supplied array will contain the
!      value of the derivative of the normalized Legendre polynomial at x
!
!



if ( x == 1.0d0) then
do i=1,n
pols(i) = 1.0d0
ders(i) = (i-1)*i/2
end do
goto 1000
endif

if ( x == -1.0d0) then
dsign = 1.0d0
do i=1,n
pols(i) = dsign
ders(i) = -dsign*(i-1)*i/2
dsign   = -dsign
end do
goto 1000
endif

pols(1)  = 1
ders(1)  = 0

if (n .gt. 1) then
pols(2)  = x
ders(2)  = 1
end if

!
!  Calculate the values of the unnormalized polynomials
!
do j=2,n-1
pols(j+1)=((2.0d0*j-1.0d0)*x*pols(j)-(j-1.0d0)*pols(j-1))/j
end do

!
!  Compute the derivatives of the unnormalized polynomials
!
d=x**2.0d0-1.0d0
do j=3,n
ders(j)=(j-1.0d0)*(x*pols(j)-pols(j-1))/d
end do

!
!  Normalize the polynomials and scale the derivatives
!
1000 continue

do i=1,n
dd      = sqrt(i-0.5d0)
pols(i) = pols(i) * dd
ders(i) = ders(i) * dd
end do

end subroutine



subroutine prolates_legen0(n,val,der)
implicit double precision (a-h,o-z)
!
!  Evaluate the Legendre function of the first kind of integer degree n
!  and its derivative at the point 0.
!
!  Input parameters:
!     n - the degree of the polynomial to evaluate
!
!  Output parameters:
!     val - the value of the polynomial at 0
!     der - the value of its derivative at 0 
!
data pi / 3.14159265358979323846264338327950288d0 /

val = 0
der = 0
!dnu = n + delta
dnu = n

! Compute cosval = cos(pi/2 * dnu) and sinval = sin(pi/2*dnu)

select case ( mod(n,4) )
case(-3)
cosvaln = 0
sinvaln = 1
case(-2)
cosvaln = -1
sinvaln = 0
case(-1)
cosvaln = 0
sinvaln = -1
case(0)
cosvaln = 1
sinvaln = 0
case(1)
cosvaln = 0
sinvaln = 1
case(2)
cosvaln = -1
sinvaln = 0
case(3)
cosvaln = 0
sinvaln = -1
end select 

! cosval = cosvaln*cos(pi/2*delta) - sinvaln*sin(pi/2*delta)
! sinval = cosvaln*sin(pi/2*delta) + sinvaln*cos(pi/2*delta)


! Compute the ratio of gamma functions by brute force if n is small, and use an asymptotic
! expansion if not.
gammaratio =  gamma(dnu/2+0.5d0)/gamma(dnu/2+1.0d0)


if (n .ge. 100) then
gammaratio = -0.104957025613225462481750387218029248d16/dnu**30-  &
0.3435274821674695538341969877678682d13/dnu**29+  &
0.137012508230061862467393693767518715d14/dnu**28+  &
0.52184449177744945865946227168534838d11/dnu**27-  &
0.20802210228419220675266715332806065d12/dnu**26-  &
0.93358732709156442705307398873761416d9/dnu**25+  &
0.37189644043206042515191168099049809d10/dnu**24+  &
0.19959724959938103205809132306364276d8/dnu**23-  &
0.79435712783055360854592176681049474d8/dnu**22-  &
0.51897923979290808045377229245471765d6/dnu**21+  &
0.20627383965581267895358219780481014d7/dnu**20+  &
0.16766118747003365058835187265431109d5/dnu**19-  &
0.66511566853934577951098546692891475d5/dnu**18-  &
0.6912167061162564516424966534975484d3/dnu**17+  &
0.27339931645089473537297049513979147d4/dnu**16+  &
0.37638891605797570836668516338219077d2/dnu**15-  &
0.14815719515134617540873450267244743d3/dnu**14-  &
0.28341526372214554549138737809124132d1/dnu**13+  &
0.11064038264157344432164622228873395d2/dnu**12+  &
0.31450601605351481364931293755092798d0/dnu**11-  &
0.12103480338672938011211980036999511d1/dnu**10-  &
0.5638860579751321228004007809251394d-1/dnu**9+  &
0.21215037666443619840288699752634574d0/dnu**8+  &
0.18752313014255059774912529012118952d-1/dnu**7-  &
0.6888076310874815972557053234370966d-1/dnu**6-  &
0.14501213286052244152751691019728349d-1/dnu**5+  &
0.5524271728019902534381596578944133d-1/dnu**4+  &
0.4419417382415922027505277263155306d-1/dnu**3-  &
0.3535533905932737622004221810524245d0/dnu**2+  &
0.14142135623730950488016887242096981d1/dnu
gammaratio = gammaratio *sqrt(dnu)
endif

if (n .le. -100) then
dnu = -dnu
gammaratio = 0.104957025613225462481750387218029248d16/dnu**30-  &
0.343527482167469553834196987767868205d13/dnu**29-  &
0.137012508230061862467393693767518715d14/dnu**28+  &
0.521844491777449458659462271685348384d11/dnu**27+  &
0.208022102284192206752667153328060649d12/dnu**26-  &
0.933587327091564427053073988737614164d9/dnu**25-  &
0.37189644043206042515191168099049809d10/dnu**24+  &
0.199597249599381032058091323063642763d8/dnu**23+  &
0.794357127830553608545921766810494744d8/dnu**22-  &
0.518979239792908080453772292454717655d6/dnu**21-  &
0.206273839655812678953582197804810138d7/dnu**20+  &
0.167661187470033650588351872654311091d5/dnu**19+  &
0.66511566853934577951098546692891475d5/dnu**18-  &
0.691216706116256451642496653497548362d3/dnu**17-  &
0.273399316450894735372970495139791469d4/dnu**16+  &
0.376388916057975708366685163382190772d2/dnu**15+  &
0.148157195151346175408734502672447426d3/dnu**14-  &
0.283415263722145545491387378091241324d1/dnu**13-  &
0.110640382641573444321646222288733953d2/dnu**12+  &
0.314506016053514813649312937550927984d0/dnu**11+  &
0.121034803386729380112119800369995106d1/dnu**10-  &
0.563886057975132122800400780925139408d-1/dnu**9-  &
0.212150376664436198402886997526345737d0/dnu**8+  &
0.187523130142550597749125290121189519d-1/dnu**7+  &
0.688807631087481597255705323437096598d-1/dnu**6-  &
0.145012132860522441527516910197283494d-1/dnu**5-  &
0.552427172801990253438159657894413312d-1/dnu**4+  &
0.44194173824159220275052772631553065d-1/dnu**3+  &
0.35355339059327376220042218105242452d0/dnu**2+  &
0.141421356237309504880168872420969808d1/dnu
!gammaratio = -gammaratio*sqrt(dnu) * sinval/cosval
gammaratio = -gammaratio*sqrt(dnu) * sinvaln/cosvaln

endif


! val        = 1.0d0/sqrt(pi) * cosval * gammaratio
! der        = 2.0d0/sqrt(pi) * sinval * 1.0d0/gammaratio

val        = 1.0d0/sqrt(pi) * cosvaln * gammaratio
der        = 2.0d0/sqrt(pi) * sinvaln * 1.0d0/gammaratio

end subroutine


subroutine prolates_lege0(dnu,val,der)
implicit double precision (a-h,o-z)
!
!  Evaluate the Legendre function of the first kind of degree dnu
!  and its derivative at the point 0.
!
!  Input parameters:
!     dnu - the degree of the polynomial to evaluate
!
!  Output parameters:
!     val - the value of the polynomial at 0
!     der - the value of its derivative at 0 
!
data pi / 3.14159265358979323846264338327950288d0 /

val = 0
der = 0

! Compute cosval = cos(pi/2 * dnu) and sinval = sin(pi/2*dnu)


cosval = cos(pi/2*dnu)
sinval = sin(pi/2*dnu)

! Compute the ratio of gamma functions by brute force if n is small, and use an asymptotic
! expansion if not.
gammaratio =  gamma(dnu/2+0.5d0)/gamma(dnu/2+1.0d0)


if (dnu .ge. 100) then
gammaratio = -0.104957025613225462481750387218029248d16/dnu**30-  &
0.3435274821674695538341969877678682d13/dnu**29+  &
0.137012508230061862467393693767518715d14/dnu**28+  &
0.52184449177744945865946227168534838d11/dnu**27-  &
0.20802210228419220675266715332806065d12/dnu**26-  &
0.93358732709156442705307398873761416d9/dnu**25+  &
0.37189644043206042515191168099049809d10/dnu**24+  &
0.19959724959938103205809132306364276d8/dnu**23-  &
0.79435712783055360854592176681049474d8/dnu**22-  &
0.51897923979290808045377229245471765d6/dnu**21+  &
0.20627383965581267895358219780481014d7/dnu**20+  &
0.16766118747003365058835187265431109d5/dnu**19-  &
0.66511566853934577951098546692891475d5/dnu**18-  &
0.6912167061162564516424966534975484d3/dnu**17+  &
0.27339931645089473537297049513979147d4/dnu**16+  &
0.37638891605797570836668516338219077d2/dnu**15-  &
0.14815719515134617540873450267244743d3/dnu**14-  &
0.28341526372214554549138737809124132d1/dnu**13+  &
0.11064038264157344432164622228873395d2/dnu**12+  &
0.31450601605351481364931293755092798d0/dnu**11-  &
0.12103480338672938011211980036999511d1/dnu**10-  &
0.5638860579751321228004007809251394d-1/dnu**9+  &
0.21215037666443619840288699752634574d0/dnu**8+  &
0.18752313014255059774912529012118952d-1/dnu**7-  &
0.6888076310874815972557053234370966d-1/dnu**6-  &
0.14501213286052244152751691019728349d-1/dnu**5+  &
0.5524271728019902534381596578944133d-1/dnu**4+  &
0.4419417382415922027505277263155306d-1/dnu**3-  &
0.3535533905932737622004221810524245d0/dnu**2+  &
0.14142135623730950488016887242096981d1/dnu
gammaratio = gammaratio *sqrt(dnu)
endif

if (n .le. -100) then
dnu = -dnu
gammaratio = 0.104957025613225462481750387218029248d16/dnu**30-  &
0.343527482167469553834196987767868205d13/dnu**29-  &
0.137012508230061862467393693767518715d14/dnu**28+  &
0.521844491777449458659462271685348384d11/dnu**27+  &
0.208022102284192206752667153328060649d12/dnu**26-  &
0.933587327091564427053073988737614164d9/dnu**25-  &
0.37189644043206042515191168099049809d10/dnu**24+  &
0.199597249599381032058091323063642763d8/dnu**23+  &
0.794357127830553608545921766810494744d8/dnu**22-  &
0.518979239792908080453772292454717655d6/dnu**21-  &
0.206273839655812678953582197804810138d7/dnu**20+  &
0.167661187470033650588351872654311091d5/dnu**19+  &
0.66511566853934577951098546692891475d5/dnu**18-  &
0.691216706116256451642496653497548362d3/dnu**17-  &
0.273399316450894735372970495139791469d4/dnu**16+  &
0.376388916057975708366685163382190772d2/dnu**15+  &
0.148157195151346175408734502672447426d3/dnu**14-  &
0.283415263722145545491387378091241324d1/dnu**13-  &
0.110640382641573444321646222288733953d2/dnu**12+  &
0.314506016053514813649312937550927984d0/dnu**11+  &
0.121034803386729380112119800369995106d1/dnu**10-  &
0.563886057975132122800400780925139408d-1/dnu**9-  &
0.212150376664436198402886997526345737d0/dnu**8+  &
0.187523130142550597749125290121189519d-1/dnu**7+  &
0.688807631087481597255705323437096598d-1/dnu**6-  &
0.145012132860522441527516910197283494d-1/dnu**5-  &
0.552427172801990253438159657894413312d-1/dnu**4+  &
0.44194173824159220275052772631553065d-1/dnu**3+  &
0.35355339059327376220042218105242452d0/dnu**2+  &
0.141421356237309504880168872420969808d1/dnu
gammaratio = -gammaratio*sqrt(dnu) * sinval/cosval

endif


val        = 1.0d0/sqrt(pi) * cosval * gammaratio
der        = 2.0d0/sqrt(pi) * sinval * 1.0d0/gammaratio

end subroutine

subroutine prolates_trieigen(n,as,bs,j,chi,y)
implicit double precision (a-h,o-z)
double precision             :: as(n),bs(n),y(n)
!
!  Calculate a single eigenvalue and a corresponding eigenvector of a
!  symmetric tridiagonal matrix via the Rokhlin-Osipov approach.
!
!  Input parameters:
!    n - the dimensions of the matrix 
!    as - the n diagonal entries of the matrix
!    bs - the n-1 subdiagonal entries of the matrix
!    j - the index of the eigenvalue to compute
!
!  Output parameters:
!    chi - the eigenvalue
!    y - the eigenvector
!

double precision, allocatable   :: as0(:),bs0(:),x(:)
!
!  Use sturm bisection to estimate  the eigenvalue
!

eps = epsilon(0.0d0)*10
call prolates_sturm_bisection(eps,n,as,bs,j,chi)
!
!  Apply the inverse power method, starting with a random right hand side
!
allocate(as0(n),bs0(n),x(n))

call random_seed()
do i=1,n
call random_number(y(i))
end do
dd = 0
do i=1,n
dd = dd + y(i)**2
end do 
dd = sqrt(dd)
y  = y / dd

chi0 = chi
do iter  = 1,6

! solve the system ( A - chi * I ) x = y 
as0 = as - chi0
bs0 = bs

! eliminate the subdiagonal entries
do i=2,n
dd     = bs0(i-1) / as0(i-1)
as0(i) = as0(i) - dd * bs0(i-1)
y(i)   = y(i)   - dd * y(i-1)
end do

x(n) = y(n) / as0(n)

do i=n-1,1,-1
x(i) = y(i) - bs0(i) * x(i+1)
x(i) = x(i) / as0(i)
end do

! normalize the eigenvector
dd   = sqrt(sum(x**2))
y    = x / dd

end do

end subroutine




subroutine prolates_sturm_bisection(eps,n,as,bs,j,dlambda)
implicit double precision (a-h,o-z)

integer, intent(in)            :: n,j
double precision, intent(in)   :: eps,as(n),bs(n)
double precision, intent(out)  :: dlambda
!
!  Use the classic Sturm bisection algorithm in order to approximate the
!  jth eigenvalue of a symmetric tridiagonal matrix.
!
!  Input parameters:
!
!    eps - the desired precision for the calculations
!    n - the dimensions of the input matrix
!    as - an array of length n speciying the diagonal entries of the input
!      matrix
!    bs - an array of length n-1 specifying the off-diagonal entries of the
!      input matrix
!    j - the index (between 1 and n) of the eigenvalue to compute
!
!
!  Output parameters:
!    dlambda - an estimate of the jth eigenvalue of the input matrix
!

double precision, allocatable :: pols(:),pols0(:)


allocate(pols(0:n))

!
!  Bracket the eigenvalue using the Gershgorin theorem
!

x1 = as(1) - abs(bs(1))
x2 = as(1) + abs(bs(1))

do i=2,n-1
x1  = min(x1, as(i) - abs(bs(i-1)) - abs(bs(i)))
x2  = max(x2, as(i) + abs(bs(i-1)) + abs(bs(i)))
end do

x1  = min(x1, as(i) - abs(bs(i-1)))
x2  = max(x2, as(i) + abs(bs(i-1)))


!
!  Perform repeated bisection
!

do while ( abs(x2-x1) .gt. eps * (abs(x1)+abs(x2)))
x       = (x1+x2)/2

if (abs(x1) + abs(x2) .lt. 1d-200) exit

!
!  Evaluate the characteristic polynomials of the principal submatrices
!  of the input matrix via the standard recurrence relation at the point x.
!
                           
pols(0) = 1.0d0
pols(1) = as(1)-x

do l=2,n

pols(l)  = (as(l)-x)*pols(l-1) - bs(l-1)**2 * pols(l-2)
dd       = pols(l)**2+pols(l-1)**2

if (dd .lt. 1.0d-13) then
pols(l)   = pols(l)   * 1.0d+24
pols(l-1) = pols(l-1) * 1.0d+24
endif

if (dd .gt. 1.0d13) then
pols(l)   = pols(l)   * 1.0d-24
pols(l-1) = pols(l-1) * 1.0d-24
endif

end do

!
!  Count the number of sign changes
!

nn = 0
do l=1,n 
if (pols(l)*pols(l-1) .lt. 0) nn = nn + 1
if (pols(l) .eq. 0)           nn = nn - 1
end do

!
!  There are nn eigenvalues less than x
!

if (j .le. nn) then
x2 = x
else
x1 = x
endif

end do

dlambda = x

end subroutine


end module 
