!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for the spheroidal wave functions ps_n(x,c) and qs_n(x,c)
!  in time independent of c and n.
!
!  It makes use of a nonstandard normalization (mainly because the
!  standard Meixner and Schafke normalization is numerically problematic).  More
!  explicitly, it evaluates
!
!    u(x) = c1 ps_n(x)   
!
!  and                                                                                    (1)
!
!    v(x) = c2 qs_n(x)
!
!  where the constants c1 and c2 are chosen so that the Wronskian of the pair is 1
!  and the function u(x)^2 + v(x)^2 is absolutely monotone on the interval
!  (-1,1).  This uniquely determines c1 and c2.
!
!  Given c and n, it solves an initial value problem for Appell's equation to construct 
!  a nonoscillatory phase function which represents ps_n(x,c) and qs_n(x,c).  The initial
!  correct initial values are computed via precomputed expansions which are contained
!  in the file prolexps.f90.   This procedure runs in O(1) time.  The functions
!  u and v can then be run in O(1) time.
!
!  Because of the use of preocmputed expansions, this code can only function in
!  the regime
!
!    2^8 <= c <= 2^20          and       200 <= n <= 3c.
!
!  The following subroutines should be regarded as publicly callable:
!
!    prolates_fast_phase - construct a nonoscillatory phase function which represents
!     the angular prolate spheroidal wave functions of order 0, specified bandlimit
!     c and characteristic exponent n
!
!    prolates_fast_eval - use the nonoscillatory phase function constructed by
!     prolates_fast_phase to evaluate the functions in (1) at a point
!
!    prolates_fast_eval0 - use the nonoscillatory phase function constructed by
!     prolates_fast_phase to evaluate the function ps_n(x,c) at a point in (0,1)
!
!    prolates_fast_evalder - use the nonoscillatory phase function constructed by
!     prolates_fast_phase to evaluate the functions in (1) and their derivatives
!     at a point
!
!    prolates_fast_legendre0 - return the values of the Legendre functions 
!     P_nu(x) and Q_nu(x) and their derivatives at the point 0; this routine is
!     useful if Flammer's normalization for the prolate functions is desired
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module prolates_fast

use utils
use chebyshev
use tensor
use prolates_window


contains


subroutine prolates_fast_phase(chebdata,c,n,chi,nints,ab,acoefs,apcoefs,appcoefs)
implicit double precision (a-h,o-z)
type(chebexps_data)                        :: chebdata
integer, intent(out)                       :: nints
double precision, allocatable, intent(out) :: ab(:,:)
double precision, allocatable, intent(out) :: acoefs(:,:),apcoefs(:,:),appcoefs(:,:)
!
!  Solve an initial value problem for Appell's equation
!
!    w'''(x) + 4 q(x) w'(x) + 2 q'(x) w(x) = 0,                                           (2)
!
!  where q(x) is the coefficient in the spheroidal wave equation, on the interval 
!  (0,\infty) using an adaptive Chebyshev spectral method in order to construct a 
!  nonoscillatory phase function representing the functions u and v defined in (1).    
!  The solution of (2) to the nonoscillatory phase via the formula
!
!    w(x) = 1 / alpha'(x),
!
!  The initial values are for (1) calculated using precomputed expansions found in 
!  prolexps.f90 and constructed by create_prolexps.f90.
!
!  The resulting phase function is represented as a piecewise Chebyshev expansion
!  on a collection of intervals.  The number of coefficients in the expansion
!  on each interval is specified by the user via the chebdata structure
!  (see chebexps in chebyshev.f90).
!
!  Input parameters:
!    chebdata - the structure returned by chebexps; the entry k in this structure
!      determine the number of Chebyshev nodes used in the piecewise 
!      discretizations scheme which represents the phase functions and its
!      derivative
!    c - the bandlimit in (1)
!    n - the characteristic exponent in (1)
!
!  Output parameters:
!    chi - the Sturm-Liouville eigenvalue of ps_n(x,c)
!    nints - the number of intervals used in the adaptive discretization of the
!        phase function and its derivative
!    ab - a (2,nints) array specifying the endpoints of the discretization
!       intervals
!    acoefs - a (k,nints) array giving the coefficients in the piecewise
!     Chebyshev expansion of the phase function
!    apcoefs - a (k,nints) array giving the coefficients in the piecewise
!     Chebyshev expansion of the derivative of the phase function
!    appcoefs - a (k,nints) array giving the coefficients in the piecewise
!     Chebyshev expansion of the second derivative of the phase function
!

double precision, allocatable     :: ab0(:,:)
double precision, allocatable     :: ws(:,:),wders(:,:),wder2s(:,:),wder3s(:,:)
double precision, allocatable     :: vals(:),vals2(:),roots(:)

double precision, allocatable     :: ab_r(:,:),coefsps(:),coefsqs(:)
double complex, allocatable       :: rs(:,:),rders(:,:),rvals(:,:)
type(c_ptr)                       :: userptr
type(prolates_qfun_data), pointer :: qdata
double complex                    :: ra,rb,zval,ima,clambda

data pi / 3.14159265358979323846264338327950288d0 /

eps = 1.0d-12
k   = chebdata%k

!
!  Find the value of chi and the phase function at 0 
!
call elapsed(t1)
call prolexps(c,n,chi,apval0,appval0,apppval0)

!
!  Compute the correct initial values for w(t) = 1/alpha'(t) and solve
!  Appell's equation going forward
!

isubst = 2
call qprolates_convert_from0(isubst,apval0,appval0,apppval0,apval,appval,apppval)
call qprolates_convert_to_w(apval,appval,apppval,wa,wpa,wppa)

call prolates_fast_appell(eps,chebdata,c,chi,isubst,wa,wpa,wppa, &
  nints,ab,ws,wders,wder2s,wder3s)

!
!  Construct the expansions of alpha(x) an dalpha'(x)
!

allocate(acoefs(k,nints),apcoefs(k,nints),appcoefs(k,nints))
!appcoefs(k,nints),apppcoefs(k,nints))
allocate(vals(k),vals2(k))

aval = -pi/2*(n+1)
do int=1,nints
a               = ab(1,int)
b               = ab(2,int)
vals            = 1.0d0/ws(:,int) 

apcoefs(:,int)  = matmul(chebdata%u,vals)
vals2           = aval + (b-a)/2 * matmul(chebdata%aintl,vals)

aval            = vals2(k)
acoefs(:,int)   = matmul(chebdata%u,vals2)

vals            = -wders(:,int)/ws(:,int)**2
appcoefs(:,int) = matmul(chebdata%u,vals)

! vals             = 2*wders(:,int)**2/ws(:,int)**3 - wder2s(:,int)/ws(:,int)**2
! apppcoefs(:,int) = matmul(chebdata%u,vals)

end do



end subroutine



subroutine prolates_fast_eval(k,nints,ab,acoefs,apcoefs,t,valps,valqs)
implicit double precision (a-h,o-z)
double precision :: ab(2,nints), acoefs(k,nints),apcoefs(k,nints)
!
!  Evaluate the functions u and v in (1) at a specified point given the
!  data generated by prolates_fast_phase.
!
!  Input parameters:
!    k - the size of the piecewise Chebyshev expansions used to represent
!      the phase function and its derivative
!    nints  - the number of intervals in the discretization scheme used
!      to describe the phase function and its derivatives
!    ab - the (2,nints) array giving the endpoints of the discretization
!      intervals
!    acoefs - the (k,nints) array specifying the piecewise expansion
!     coefficients of the phase function
!    apcoefs - the (k,nints) array specifying the piecewise expansion
!     coefficients of the derivative of the phase function
!    t - the point at which to evaluate the functions (1)
!    
!
!  Output parameters:
!    valps - the value of u at t
!    valqs - the value of v at t
!
!



x = -log(1-t)

call chebpw_eval22(nints,ab,k,acoefs,apcoefs,x,aval,apval)
dd     = 1.0d0/sqrt(apval) * 1.0d0/sqrt(1+t)

valps  = sin(aval)*dd
valqs  = cos(aval)*dd

end subroutine


subroutine prolates_fast_eval0(k,nints,ab,acoefs,apcoefs,t,valps)
implicit double precision (a-h,o-z)
double precision :: ab(2,nints), acoefs(k,nints),apcoefs(k,nints)
!
!  Evaluate only the function u in (1).
!
!  Input parameters:
!    k - the size of the piecewise Chebyshev expansions used to represent
!      the phase function and its derivative
!    nints  - the number of intervals in the discretization scheme used
!      to describe the phase function and its derivatives
!    ab - the (2,nints) array giving the endpoints of the discretization
!      intervals
!    acoefs - the (k,nints) array specifying the piecewise expansion
!     coefficients of the phase function
!    apcoefs - the (k,nints) array specifying the piecewise expansion
!     coefficients of the derivative of the phase function
!    t - the point at which to evaluate the functions (1)
!    
!
!  Output parameters:
!    valps - the value of u at t
!
x = -log(1-t)
call chebpw_eval22(nints,ab,k,acoefs,apcoefs,x,aval,apval)

dd     = 1.0d0/sqrt(apval) * 1/sqrt(1+t)
valps  = sin(aval)*dd

end subroutine


subroutine prolates_fast_evalder(k,nints,ab,acoefs,apcoefs,appcoefs, &
  t,valps,derps,valqs,derqs)
implicit double precision (a-h,o-z)

double precision :: ab(2,nints), acoefs(k,nints), apcoefs(k,nints),appcoefs(k,nints)

!
!  Evaluate the functions u and v in (1) and their derivatives at a specified 
!  point given the data generated by prolates_fast_phase.
!
!  Input parameters:
!    k - the size of the piecewise Chebyshev expansions used to represent
!      the phase function and its derivative
!    nints  - the number of intervals in the discretization scheme used
!      to describe the phase function and its derivatives
!    ab - the (2,nints) array giving the endpoints of the discretization
!      intervals
!    acoefs - the (k,nints) array specifying the piecewise expansion
!     coefficients of the phase function
!    apcoefs - the (k,nints) array specifying the piecewise expansion
!     coefficients of the derivative of the phase function
!    appcoefs - the (k,nints) array specifying the piecewise expansion
!     coefficients of the second derivative of the phase function
!    t - the point at which to evaluate the functions (1)
!    
!
!  Output parameters:
!    valps - the value of u at t
!    derps - the derivative of u at t
!    valqs - the value of v at t
!    derqs - the derivative of v at t
!
!

x = -log(1-t)
call chebpw_eval23(nints,ab,k,acoefs,apcoefs,appcoefs,x,aval,apval,appval)

!
!  Compute the values of the functions under the change of variables
!
valps0 = sin(aval)/sqrt(apval)
derps0 = cos(aval)*sqrt(apval) - sin(aval)*appval/(2*apval**1.5d0)

valqs0 = cos(aval)/sqrt(apval)
derqs0 = -sin(aval)*sqrt(apval) - cos(aval)*appval/(2*apval**1.5d0)

!
!  Undo the change of variables
!

dd    = sqrt(2.0d0 - exp(-x))
valps =  valps0 / dd 
derps = derps0/sqrt(1+t) - (1-t)/2 * valps

valqs = -valqs0 / dd
derqs = derqs0/sqrt(1+t) - (1-t)/2 * valqs

end subroutine



subroutine prolates_fast_appell(eps,chebdata,c,chi,isubst,wa,wpa,wppa, &
  nints,ab,ws,wders,wder2s,wder3s)
implicit double precision (a-h,o-z)
type(chebexps_data)                        :: chebdata
integer, intent(out)                       :: nints
double precision, allocatable, intent(out) :: ab(:,:), ws(:,:),wders(:,:),wder2s(:,:),wder3s(:,:)
double precision, intent(in)               :: wa,wpa,wppa
type(c_ptr)                                :: userptr
!
!  Solve an initial value problem for Appell's equation
!
!    w'''(x) + 4 q(x) w'(x) + 2 q'(x) w(x) = 0,
!
!  where q(x) is the coefficient in the spheroidal wave equation, on the interval (0,\infty) 
!  using an adaptive Chebyshev spectral method.
!
!  This is a rather specialized version of the solve3_ivp_adap routine in linearode.f90.
!  It assumes that w is increasing (which it is if the initial conditions are properly
!  set) and trys to extend the solution as far as possible before it becomes too big.
!
!  Input parameters:
!    eps - precision for the adaptive discretization procedures
!    chebdata - the structure returned by the chebexps routine; the number of
!      Chebyshev nodes k used per discretization interval is specified by this
!      data structure
!    (c,chi) - the parameters in (2)
!    isubst - an integer parameter specifying the transformation (3); see 
!       qprolates
!    wa - the value of w(0)
!    wpa - the value of w'(0)
!    wppa -  the value of w''(0)
!
!  Output parameters:
!    nints - the number of intervals used to discretize the solution
!    ab - a (2,nints) array specifying the intervals used to discretize the
!      solution
!    ws - the (k,nints) array of values of w at the discretization nodes
!    wders - the (k,nints) array of the value of w' at the discretization nodes
!    wder2s - the (k,nints) array of the value of w'' at the discretization nodes
!    wder3s - the (k,nints) array of the value of w''' at the discretization nodes
!

double precision, allocatable   :: ab0(:,:)               ! the list of intervals to process


! data for solving over a single interval
double precision, allocatable   :: ys0(:),yders0(:),yder2s0(:),yder3s0(:),coefs0(:)
double precision, allocatable   :: amatr(:,:),ts(:),ps(:),rs(:),qs(:),fs(:)

! output data
double precision, allocatable   :: about(:,:)
double precision, allocatable   :: ysout(:,:),ydersout(:,:),yder2sout(:,:)
double precision, allocatable   :: yder3sout(:,:)


ier = 0

!
!  Set algorithm parameters and allocate memory for the procedure.
!

maxints  = 1000
k        = chebdata%k
nn       = k/4

allocate(ys0(k),yders0(k),yder2s0(k),yder3s0(k),coefs0(k),amatr(k,k))
allocate(ab0(2,maxints),rs(k),qs(k),ps(k),ts(k),fs(k)) 
allocate(about(2,maxints))
allocate(ysout(k,maxints))
allocate(ydersout(k,maxints))
allocate(yder2sout(k,maxints))
allocate(yder3sout(k,maxints))

!
!  Form an initial set of intervals on a large interval
!

a = 0
b = 10

dd      = (b-a)
nints0  = 11
do i=1,nints0
int        = nints0-i+1
ab0(1,int) = a + (b-a)*(i-1.0d0)/(nints0+0.0d0)
ab0(2,int) = a + (b-a)*(i)/(nints0+0.0d0)
end do

! dd     = log(b-a)
! nints0 = 20
! dd     = dd/(nints0+0.0d0)

! do i=1,nints0
! int       = nints0-i+1
! a0        = dd * (i-1)
! b0        = dd * i
! ab0(1,int) = a0
! ab0(2,int) = b0
! end do

! ab0(:,1:nints0)        = a + exp(ab0(:,1:nints0))-

!call prin2("ab0 = ",ab0(:,1:nints0))
! stop

nintsout = 0

do while (nints0 > 0) 

ifsplit = 0
a0      = ab0(1,nints0)
b0      = ab0(2,nints0)
nints0  = nints0 - 1


ts = (a0+b0)/2 + (b0-a0)/2 * chebdata%xs

do i=1,k
call qprolates(isubst,c,chi,ts(i),qval,qder)
ps(i) = 4*qval
qs(i) = 2*qder
end do


if (nintsout .eq. 0) then
ys0(1)     = wa
yders0(1)  = wpa
yder2s0(1) = wppa
else
ys0(1)     = ysout(k,nintsout)
yders0(1)  = ydersout(k,nintsout)
yder2s0(1) = yder2sout(k,nintsout)
endif



call prolates_fast_int(a0,b0,k,chebdata%xs,chebdata%aintl,chebdata%aintl2,chebdata%aintl3,  &
 ps,qs,ys0,yders0,yder2s0,yder3s0)

if (abs(ys0(k)) .gt. 1d15) exit

ifsplit = 0

coefs0 = matmul(chebdata%u,ys0)
dd1 = maxval(abs(coefs0))
dd2 = maxval(abs(coefs0(k-nn+1:k)))
if (dd2 .gt. eps*dd1) ifsplit = 1

coefs0 = matmul(chebdata%u,1/ys0)
dd1 = maxval(abs(coefs0))
dd2 = maxval(abs(coefs0(k-nn+1:k)))
if (dd2 .gt. eps*dd1) ifsplit = 1

if (ifsplit .eq. 1) then
if (nints0+2 .gt. maxints) then
ier = 8
return
endif

nints0 = nints0+1
ab0(1,nints0) = (a0+b0)/2
ab0(2,nints0) = b0

nints0 = nints0+1
ab0(1,nints0) = a0
ab0(2,nints0) = (a0+b0)/2


else

if (nintsout+1 .gt. maxints) then
ier = 8
return
endif

nintsout          = nintsout+1
about(1,nintsout) = a0
about(2,nintsout) = b0

ysout(:,nintsout)     = ys0
ydersout(:,nintsout ) = yders0
yder2sout(:,nintsout) = yder2s0
yder3sout(:,nintsout) = yder3s0

endif


end do
!
!  Copy out the data
!
nints = nintsout
allocate(ab(2,nints))
allocate(ws(k,nints))
allocate(wders(k,nints))
allocate(wder2s(k,nints))
allocate(wder3s(k,nints))

ab    = about(:,1:nints)
ws    = ysout(:,1:nints)
wders = ydersout(:,1:nints)
wder2s = yder2sout(:,1:nints)
wder3s = yder3sout(:,1:nints)

end subroutine



subroutine prolates_fast_int(a,b,k,xscheb,chebint,chebint2,chebint3,ps,qs,ys,yders,yder2s,yder3s)
implicit double precision (a-h,o-z)

integer, intent(in)           :: k
double precision, intent(in)  :: xscheb(k),chebint(k,k),chebint2(k,k),chebint3(k,k)
double precision, intent(in)  :: ps(k),qs(k)
double precision, intent(out) :: ys(k),yders(k),yder2s(k),yder3s(k)

!
!  Solve an initial value for the ordinary differential equation
!
!     y'''(t) + p(t) y'(t) + q(t) y(t) = 0                                        (5)
!
!  on the interval [a,b] using a standard Chebyshev spectral method.
!
!  Input parameters:
!
!    (a,b) - the interval on which the ODE (5) is given
!    k - the number of Chebyshev points on the interval [a,b]
!    xscheb - the nodes of the k-point Chebyshev grid on [-1,1] as returned by
!     the subroutine chebexps (see below)
!    chebint? - the "left" Chebyshev spectral integration matrice as constructed
!      by the subroutine chebexps 
!    ps - an array specifying the values of the function p(t) appearing in (5)
!      at the k Chebyshev nodes on [a,b]
!    qs - an array specifying the values of the function q(t) appearing in (5)
!      at the k Chebyshev node on [a,b]
!
!    ys(1) - the value of y(a)
!    yders(1) - the value of y'(a)
!    yder2s(1) - the value of y''(a)
!
!  Output parameters:
!
!    ys - the values of the solution y of (3) at the nodes of the k-point
!      Chebyshev grid on [a,b]
!    yders - the values of the solution y' of (3) at the nodes of the k-point
!      Chebyshev grid on [a,b]
!    yder2s - the values of the solution y'' of (3) at the nodes of the k-point
!      Chebyshev grid on [a,b]
!    yder3s - the values of the solution y''' of (3) at the nodes of the k-point
!      Chebyshev grid on [a,b]

!

double precision, allocatable :: amatr(:,:),xs(:),sigma(:),rhs(:)

!
!  Allocate memory for the procedure and setup some parameters.
!


allocate(amatr(k,k),xs(k),sigma(k),rhs(k))


xs       = (b-a)/2 *xscheb + (b+a)/2 
alpha    = ys(1)
beta     = yders(1)
eta      = yder2s(1)


!
!  We represent the solution in the form
!
!      y(t) = alpha + beta (t-a) + eta(t-a)^2/2 + \int_a^t (t-s)^2/2 sigma(s) ds,
!
!  insert this representation into (3), and solve the resulting system of
!  linear equations in order to obtain the values of sigma.  
!    

amatr  = 0
do i=1,k
amatr(i,i) = 1.0d0
sigma(i)   = 0.0d0
end do

!
! Handle the r(t) * y''(t) term.
!
! do i=1,k
! amatr(i,:) = amatr(i,:) + rs(i) * chebint(i,:)*(b-a)/2
! sigma(i)   = sigma(i) - rs(i)*eta
! end do

!
! Handle the p(t) * y'(t) term.
!
do i=1,k
amatr(i,:) = amatr(i,:) + ps(i) * chebint2(i,:)*((b-a)/2)**2
sigma(i)   = sigma(i) - ps(i)*(beta + eta * (xs(i)-a))
end do

!
!  Handle the q(t) y(t) term
!

do i=1,k
amatr(i,:) = amatr(i,:) + qs(i) * chebint3(i,:)*((b-a)/2)**3
sigma(i)   = sigma(i) - qs(i) * (alpha + beta*(xs(i)-a) + eta * (xs(i)-a)**2/2)
end do


! !
! !  Form the right-hand side.
! !
! do i=1,k
! sigma(i) = sigma(i) + fs(i)
! end do

!
!  Use a QR decomposition to invert the linear system
!
call qrsolv(amatr,k,sigma)
!call gesv(amatr,sigma)

!
!  Calculate y(t), y'(t) and y''(t) from sigma.
!

yder3s = sigma
yder2s = (b-a)/2*matmul(chebint,sigma)
yders  = ((b-a)/2)**2*matmul(chebint2,sigma)
ys     = ((b-a)/2)**3*matmul(chebint3,sigma)


do i=1,k
ys(i)      = ys(i) + alpha + beta*(xs(i)-a) + eta * (xs(i)-a)**2/2
yders(i)   = yders(i) + beta   + eta * (xs(i)-a)
yder2s(i)  = yder2s(i) + eta
end do


end subroutine



subroutine prolates_fast_legendre0(n,valp,derp,valq,derq)
implicit double precision (a-h,o-z)
!
!  Evaluate the Legendre polynomial P_n(x) and its derivative,
!  and the Legendre function of the second kind Q_n(x) and it
!  derivative at 0.
!
!  Input parameters:
!    n - the degree of the Legendre polynomial
!
!  Output parameters:
!    valp - the value of P_n(0)
!    derp - the value of P_n'(0)
!    valq - the value of Q_n(0)
!    derq - the value of Q_n'(0)
!

data pi / 3.14159265358979323846264338327950288d0 /

x = n
k = n/2

if (n .le. 100) then


! even and small n
if (k*2 .eq. n) then

valp = sqrt(pi) / gamma(x/2+1) * 1/gamma(0.5d0-x/2)
derp = 0

valq = 0
derq = (-1)**k * sqrt(pi) * gamma(x/2+1)/gamma(x/2+0.5d0) 

! odd and small n
else

valp = 0
derp = -2*sqrt(pi)/gamma(x/2+0.5d0) * 1/gamma(-x/2)

valq = (-1)**k*0.5d0*sqrt(pi)*gamma(x/2+0.5d0)/gamma(x/2+1)
derq = 0

endif

return
endif


! even and n large
if (k*2 .eq. n) then

valp = -0.148431649088563762198660022422650598d16/x**30-  &
0.485821224329117000978403785616984106d13/x**29+  &
0.193764947353708790444903653063270543d14/x**28+  &
0.73799955772136410633236984018979611d11/x**27-  &
0.294187678323667813407880442715722451d12/x**26-  &
0.132029185963253723789708229216300673d10/x**25+  &
0.52594098985729771344549331741203939d10/x**24+  &
0.282273137395812467922073432224222556d8/x**23-  &
0.112339062354570726254759134923877362d9/x**22-  &
0.73394747950520926359471024502667547d6/x**21+  &
0.291715261604023450980605058902028759d7/x**20+  &
0.237108725203699616179697784446034348d5/x**19-  &
0.94061559899519089480435241057421081d5/x**18-  &
0.97752804032846775328380317660048604d3/x**17+  &
0.386645021268388973112450912594795227d4/x**16+  &
0.53229430981609766604378819465637207d2/x**15-  &
0.209525914746191119775176048278808594d3/x**14-  &
0.4008097097394056618213653564453125d1/x**13+  &
0.1564691296778619289398193359375d2/x**12+  &
0.444778673350811004638671875d0/x**11-  &
0.1711690604686737060546875d1/x**10-  &
0.797455310821533203125d-1/x**9+0.30002593994140625d0/x**8+  &
0.26519775390625d-1/x**7-0.97412109375d-1/x**6-  &
0.205078125d-1/x**5+0.78125d-1/x**4+0.625d-1/x**3-0.5d0/x**2  &
+0.2d1/x

valp = valp * (-1)**k * sqrt(x/(2*pi))
derp = 0

valq = 0

derq = 0.314159265358979323846264338327950288d1+  &
0.58434230005548153042140794152260281d15/x**30+  &
0.233155889168425022325822071642256577d16/x**29-  &
0.76312619465517644873105170657110965d13/x**28-  &
0.304365267564812292081756713693344732d14/x**27+  &
0.115924699444497702336208285815639048d12/x**26+  &
0.462108924499136031064547262831559196d12/x**25-  &
0.207390960340799273896746340987369073d10/x**24-  &
0.82614617497871522732112250861911224d10/x**23+  &
0.443393607374213395012735496796180248d8/x**22+  &
0.176461786502142546981392300590730238d9/x**21-  &
0.115288200486715537926777701497525336d7/x**20-  &
0.458225261397612378845394349552890901d7/x**19+  &
0.372449514600991882757348738207674927d5/x**18+  &
0.147751552782762730886485108642108393d6/x**17-  &
0.153549745508697071452859057156106626d4/x**16-  &
0.60734057918192007915281554504711202d4/x**15+  &
0.83612594663295089678764321799727409d2/x**14+  &
0.329122537251657674709769935918117915d3/x**13-  &
0.62959041980238711423954008194179508d1/x**12-  &
0.245781134154779863708940998223948068d2/x**11+  &
0.69865670633616109895954676238959888d0/x**10+  &
0.26887173144512620304705748105919643d1/x**9-  &
0.125263987302154692771586128581664358d0/x**8-  &
0.471279644403147198114560894193521312d0/x**7+  &
0.416571657710194446445318161021406251d-1/x**6+  &
0.153014583591592712664379237051821491d0/x**5-  &
0.322135965455984658240798393793308401d-1/x**4-  &
0.122718463030851298377447007159355581d0/x**3+  &
0.98174770424681038701957605727484465d-1/x**2+  &
0.78539816339744830961566084581987572d0/x
derq = derq * (-1)**k * sqrt(x/(2*pi))

! odd and n large 
else

valp = 0
derp = 0.2d1+0.372003862046069566717328990006299221d15/x**30+  &
0.148431649088563762198660022422650598d16/x**29-  &
0.485821224329117000978403785616984106d13/x**28-  &
0.193764947353708790444903653063270543d14/x**27+  &
0.73799955772136410633236984018979611d11/x**26+  &
0.294187678323667813407880442715722451d12/x**25-  &
0.132029185963253723789708229216300673d10/x**24-  &
0.52594098985729771344549331741203939d10/x**23+  &
0.282273137395812467922073432224222556d8/x**22+  &
0.112339062354570726254759134923877362d9/x**21-  &
0.73394747950520926359471024502667547d6/x**20-  &
0.291715261604023450980605058902028759d7/x**19+  &
0.237108725203699616179697784446034348d5/x**18+  &
0.94061559899519089480435241057421081d5/x**17-  &
0.97752804032846775328380317660048604d3/x**16-  &
0.386645021268388973112450912594795227d4/x**15+  &
0.53229430981609766604378819465637207d2/x**14+  &
0.209525914746191119775176048278808594d3/x**13-  &
0.4008097097394056618213653564453125d1/x**12-  &
0.1564691296778619289398193359375d2/x**11+  &
0.444778673350811004638671875d0/x**10+  &
0.1711690604686737060546875d1/x**9-  &
0.797455310821533203125d-1/x**8-0.30002593994140625d0/x**7+  &
0.26519775390625d-1/x**6+0.97412109375d-1/x**5-  &
0.205078125d-1/x**4-0.78125d-1/x**3+0.625d-1/x**2+0.5d0/x
derp = (-1)**k * sqrt(x/(2*pi)) * derp

valq = 0.233155889168425022325822071642256577d16/x**30+  &
0.76312619465517644873105170657110965d13/x**29-  &
0.304365267564812292081756713693344732d14/x**28-  &
0.115924699444497702336208285815639048d12/x**27+  &
0.462108924499136031064547262831559196d12/x**26+  &
0.207390960340799273896746340987369073d10/x**25-  &
0.82614617497871522732112250861911224d10/x**24-  &
0.443393607374213395012735496796180248d8/x**23+  &
0.176461786502142546981392300590730238d9/x**22+  &
0.115288200486715537926777701497525336d7/x**21-  &
0.458225261397612378845394349552890901d7/x**20-  &
0.372449514600991882757348738207674927d5/x**19+  &
0.147751552782762730886485108642108393d6/x**18+  &
0.153549745508697071452859057156106626d4/x**17-  &
0.60734057918192007915281554504711202d4/x**16-  &
0.83612594663295089678764321799727409d2/x**15+  &
0.329122537251657674709769935918117915d3/x**14+  &
0.62959041980238711423954008194179508d1/x**13-  &
0.245781134154779863708940998223948068d2/x**12-  &
0.69865670633616109895954676238959888d0/x**11+  &
0.26887173144512620304705748105919643d1/x**10+  &
0.125263987302154692771586128581664358d0/x**9-  &
0.471279644403147198114560894193521312d0/x**8-  &
0.416571657710194446445318161021406251d-1/x**7+  &
0.153014583591592712664379237051821491d0/x**6+  &
0.322135965455984658240798393793308401d-1/x**5-  &
0.122718463030851298377447007159355581d0/x**4-  &
0.98174770424681038701957605727484465d-1/x**3+  &
0.78539816339744830961566084581987572d0/x**2-  &
0.314159265358979323846264338327950288d1/x
valq = valq*(-1)**k * sqrt(x/(2*pi))
derq = 0

endif

end subroutine


end module
