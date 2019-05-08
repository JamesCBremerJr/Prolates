!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This file contains code for evaluating the function
!
!    f(x) =  alpha(x,c,chi)
!
!  and its first three derivatives at 0, where alpha(x,c,chi) is a nonoscillatory
!  phase function for one of several differential equations of the form
!
!    z''(x) + q(x) z(x) = 0                                                               (1)
!
!  obtained from the spheroidal wave equation 
!
!    (1-t^2) y''(t) - 2 t y'(t) + (chi - c^2 t^2) y(t) = 0                                (2)
!
!  via transformations of the type
!
!    z(x) = y(\psi(x)) r(x).                                                              (3)
!
!  If chi_n(c) denotes the Sturm-Liouville eigenvalue of the prolate function
!  ps_n(x,c), then the value of alpha(0,c,chi_n(c)) is -pi/2 ( n + 1).
!
!  The code here makes use of the ``windowing'' algorithm and achieves 
!  double precision accuracy roughly in the regime
!
!    c >= 250, chi > chi_100(c).
!
!  The following subroutines should be regarded as publicly callable:
!
!    prolates_window_avals - compute the values of a nonoscillatory phase function
!      alpha(x) for one of the equations (1) and its first three derivatives at the 
!      point 0 given chi and c using the "windowing algorithm"
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module prolates_window

use utils
use chebyshev
use riccati
use linearode
use prolates
use tensor
use ieee_arithmetic

! a data structure which passes information to the solver for Riccati's equation
type         prolates_qfun_data
integer                    :: ifwindow, isubst
double precision           :: a,b,dlambda,c,chi,eps
double precision           :: chi1, chi2
end type     prolates_qfun_data

double precision, parameter, private     :: dtail = .25d0

contains


subroutine prolates_window_avals(c,chi,isubst,aval,apval,appval,apppval)
implicit double precision (a-h,o-z)
!
!  Given c and chi, compute the values of the nonoscillatory phase function
!  \alpha(x) for one of the equations (1) and its first three derivatives at the
!  point 0.  The transformation (3) is determined via the isubst parameter.
!
!  Input parameters:
!    eps - precision for the adaptive discretization procedures
!    (c,chi) - the values of the parameters in (2)
!    isubst - an integer parameter specifying the operative transformation/
!      change of variables
!
!      isubst = 0    means   \psi(x) = x,              r(x) = \sqrt(1-x^2)
!      isubst = 1    menas   \psi(x) = w/sqrt(1+w^2),  r(x) = (1+w^2)^(1/4)
!      isubst = 2    means   \psi(x) = 1 - exp(-w),    r(x) =  sqrt(2 - \exp(-x))
!
!  Output parameters:
!    aval - the value of \alpha(0)
!    apval - the value of \alpha'(0)
!    appval - the value of \alpha''(0)
!    apppval - the value of \alpha'''(0)
! 
!
type(chebexps_data)                    :: chebdata
type(c_ptr)                            :: userptr
type(prolates_qfun_data), pointer      :: qdata
double precision, allocatable          :: roots(:),ab(:,:)
double complex, allocatable            :: rs(:,:),rders(:,:)
double precision, allocatable          :: ab_r(:,:)
double complex                         :: ra, rb, ima
double precision, allocatable          :: ws(:,:), wders(:,:), wder2s(:,:), &
                                          wder3s(:,:), ab_w(:,:)

eps0 = epsilon(0.0d0)


if (eps0 .lt. 1.0d-20) then
eps   = 1.0d-20
elseif (eps0 .lt. 1.0d-17) then
eps   = 1.0d-15
else
eps   = 1.0d-13
endif


k    = 30
call chebexps(k,chebdata)

isubst1 = 0
a      =  0.0d0
b      =  1.0d0
call qprolates_find_roots(isubst1,c,chi,a,b,nroots,roots)

!call prin2("roots = ",roots)

if (nroots .eq. 0) then
tp = 0.90d0
else
tp = roots(1)
endif

if (tp .ge. 1.0d0) then
tp = 0.90d0
endif

a              = 0
b              = tp

ima            = (0.0d0,1.0d0)
pi             = acos(-1.0d0)
dlambda        = sqrt(chi + 0.25d0)
ifwindow       = 2

allocate(qdata)
qdata%a        = a
qdata%b        = b
qdata%dlambda  = dlambda
qdata%c        = c
qdata%chi      = chi
qdata%ifwindow = ifwindow
qdata%isubst   = isubst1

userptr        = c_loc(qdata)


!  Solve the windowed Riccati equation backward
!

rb = ima*dlambda 
call riccati_tvp_adap(ier,eps/50,chebdata,a,b,prolates_ricode,rb,nints_r,ab_r,rs,rders,userptr)

if (ier .ne. 0) then
call prini("after riccati_tvp_adap, ier = ",ier)
call prind("c       = ",c)
call prind("chi     = ",chi)
call prind("tp      = ",tp)
call prini("isubst1 = ",isubst1)
call prin2("eps     = ",eps)
stop
endif


! call prin2("ab_r = ",ab_r)
!
!  Compute the values of the phase function from the values of r at 0
!

apval0   = imag(rs(1,1))
appval0  = imag(rders(1,1))
apppval0 = -(real(rders(1,1)) *2*apval0-appval0**2/apval0)

!
!  Apply the change of variables
!

isubst2 = 2
call qprolates_convert_from0(isubst2,apval0,appval0,apppval0,apval,appval,apppval)

!
!  Convert the values of the first three derivatives of \alpha(x) to
!  the values of w(0), w'(0) and w''(0)
!

ifwindow       = 0

qdata%isubst   = isubst2
qdata%ifwindow = 0
qdata%a        = a 
qdata%b        = b

call qprolates_convert_to_w(apval,appval,apppval,wa,wpa,wppa)

!
!  Solve Appell's equation going forward using an adaptive Chebyshev spectral
!  method
!

call prolates_window_appell(eps,chebdata,c,chi,isubst2,wa,wpa,wppa, &
  nints_w,ab_w,ws,wders,wder2s,wder3s)

if (ier .ne. 0) then
call prini("after prolates_window_appell, ier = ",ier)
stop
endif

! call prin2("ab_w = ",ab_w)

!
!  Extract the values of alpha'(0), alpha''(0) and alpha'''(0)
!

wa   = ws(1,1)
wpa  = wders(1,1)
wppa = wder2s(1,1)

call qprolates_convert_from_w(wa,wpa,wppa,apval,appval,apppval)

!
!  Integrate \alpha'(x) to calculate the value of alpha(0)
!


! call chebpw_plot(ifshow,"w1.pdf",nints_w,ab_w,k,chebdata%xs,ws)
! call chebpw_plot(ifshow,"ap.pdf",nints_w,ab_w,k,chebdata%xs,1/ws)

aval = 0

do int = nints_w,1,-1
a0   = ab_w(1,int)
b0   = ab_w(2,int)
aval = aval + (b0-a0)/2*dot_product(chebdata%aintr(1,:),1.0d0/ws(:,int))
end do

!
!  Compute the values of the derivatives of the phase function alpha0(t)
!  with respect to the default change of variables and then to the 
!  change of variables desired by the user
!

call qprolates_convert_to0(isubst2,apval,appval,apppval,  &
  apval0,appval0,apppval0)

call qprolates_convert_from0(isubst,apval0,appval0,apppval0,  &
  apval,appval,apppval)

end subroutine



subroutine prolates_winfun(a,b,t,val,der)
implicit double precision (a-h,o-z)
!
!  Return the value of a smooth window function which is, up to double
!  precision accuracy, equal to 1 when t = a and 0 when t = b.
!
!  That is the graph of the function has the following shape:
!
!    
!   ------------\
!                ---
!                   \
!                    ---
!                       \--------------  
!                
!                     
!   a                                 b
!                             
!
!
!
!  Input parameters:
!    (a,b) the interval on which the window function is given
!    t - point at which to evaluate it 
!
!  Output parameters:
!    val - the value of the window function
!    der - the derivative of the window function
!
!

data dd / 0.564189583547756286948079451560772586d0 /

c     = (b+a)/2
alpha = 8d0/(b-c)
val   = ( 1-erf(alpha*(t-c)) ) / 2
der   = -dd*alpha*exp(-alpha**2*(t-c)**2)

end subroutine




subroutine qprolates(isubst,c,chi,x,val,der)
implicit double precision (a-h,o-z)
!
!  Evaluate the coefficient q(x) in one of several differential equations of the form
!
!     z''(x) + q(x) z(x) = 0                                                              (1)
!
!  obtained from the spheroidal wave equation
!
!    (1-t^2) * y''(t) - 2 t y'(t) + (chi - c^2 t^2) y(t) = 0,   0 < t < 1                 (2)
!
!  via a change of variables of the type
!
!    z(x) = y( \psi(x) )  r(x).                                                           (3)
!
!  The particular change of variables is specified via the parameter isubst.
!
!  Input parameters:
!    isubst - specify the operative change of variables
!
!      isubst = 0    means   \psi(x) = x,              r(x) = \sqrt(1-x^2)
!      isubst = 1    menas   \psi(x) = w/sqrt(1+w^2),  r(x) = (1+w^2)^(1/4)
!      isubst = 2    means   \psi(x) = 1 - exp(-w),    r(x) =  sqrt(2 - \exp(-x))
!
!    (c,chi) - the parameters in (2)
!    x - the point at which to evaluate q(x)
!
!  Output parameters:
!    val - the value of q(x) 
!    der - the value of q'(x)
!

select case (isubst)

case(0)

val   = (1+chi-0.1d1*(c**2+chi)*x**2+c**2*x**4)/(-1+x**2)**2
der   = (-0.2d1*x*(2+chi-chi*x**2+c**2*(-1+x**2)))/(-1+x**2)**3

case(1) 

val   = (-2+0.4d1*chi+(-1-0.4d1*c**2+0.4d1*chi)*x**2+x**4)/(0.4d1*(1  &
         +x**2)**3)

der   = (0.2d1*(-1-4*c**2+4*chi)*x+0.4d1*x**3)/(0.4d1*(1+x**2)**3)-  &
        (3*x*(-2+0.4d1*chi+(-1-0.4d1*c**2+0.4d1*chi)*x**2+  &
        x**4))/(0.2d1*(1+x**2)**4)

case(2)

z   = exp(x)
 if (z .gt. 10**7) then
 val = (1-0.1d1*c**2+chi)/(0.8d1*z**3)+(3+0.12d2*c**2+  &
 0.4d1*chi)/(0.16d2*z**2)+(0.5d0-0.1d1*c**2+chi)/(0.2d1*z)

 der = (3*(-1+c**2-0.1d1*chi))/(0.8d1*z**3)+(-3-0.12d2*c**2-  &
 0.4d1*chi)/(0.8d1*z**2)+(-1+0.2d1*c**2-0.2d1*chi)/(0.4d1*z)
 else

val = c**2/z**2 + 0.25d0 / (2*z-1)**2 + (0.5d0 - c**2 + chi) / (2*z-1)

der = (-1+0.2d1*c**2-0.2d1*chi)/(1-0.2d1*z)**2-(0.2d1*c**2)/z**3-  &
      0.1d1/(-1+2*z)**3

der = der * z
 endif


end select

end subroutine



subroutine qprolates_find_roots(isubst,c,chi,a,b,nroots,roots)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out) :: roots(:)
!
!  Attempt to find all simple roots of the function q(x) in a specified interval 
!  via a randomized procedure.
!
!  Input parameters:
!    isubst - the operative substitution
!    (a,b) - the interval on which to 
!
!  Output parameters:  
!    nroots - the numbers of roots
!    roots - an array containing the roots
!

double precision, allocatable  :: ts(:),vals(:)
double precision, allocatable  :: ab(:,:)

nbisect  = 100
nnewt    = 8

maxroots = 1000
nn       = maxroots*4

allocate(vals(nn),ts(nn),ab(2,maxroots))

!
!  Use either randomized brackets or equal sampled ones.
!

ts(1) = a
do i=2,nn-1
!call random_number(ts(i))
ts(i) = 1.0d0/(nn+0.0d0)*i
end do

ts(nn) = b

!call quicksort(nn,ts)

do i=1,nn
call qprolates(isubst,c,chi,ts(i),vals(i),der)
end do
m = 0
do i=1,nn-1
if (vals(i) * vals(i+1) .le. 0) then

m = m +1
if (m .gt. maxroots) then
call prina("find_roots failed")
stop
endif


ab(1,m) = ts(i)
ab(2,m) = ts(i+1)
endif
end do

!call prin2("ab = ",ab(:,1:m))

nroots = m
allocate(roots(nroots))

!
!  Use bisection to refine the brackets
!


! do iter=1,nbisect
! do int=1,m

! a0 = ab(1,int)
! b0 = ab(2,int)
! c0 = (a0+b0)/2

! call qprolates(isubst,c,chi,a0,val1,der)
! call qprolates(isubst,c,chi,c0,val2,der)

! if (val1*val2 <= 0) then
! ab(2,int) = c0
! else
! ab(1,int) = c0
! endif

! end do
! end do

!
!  Use Newton to refine the roots
!
do int=1,m

t = (ab(1,int)+ab(2,int))/2
do iter=1,nnewt
call qprolates(isubst,c,chi,t,val,der)
t = t - val/der
end do

roots(int) = t
end do


end subroutine


subroutine qprolates_convert_from0(isubst,apval0,appval0,apppval0, &
  apval,appval,apppval)
implicit double precision (a-h,o-z)
!
!  Given the values of the first three derivatives of a phase function \alpha0(t)
!  for the differential (1) corresponding to isubst = 0, compute the
!  first three derivatives of the phase function \alpha(x) corresponding to (1) 
!  under a different change of variables.
!
!  Input parameters:
!    isubst - integer parameter specifying the change of variables (see qprolates)
!    apval0 - the value of \alpha0'(0)
!    appval0 - the value of \alpha0''(0)
!    apppval0 - the value of \alpha0'''(0)
!
!  Output parameters:  
!    apval - the value of alpha'(0)
!    appval - the value of alpha''(0)
!    apppval - the value of alpha'''(0)
!

select case(isubst)

case(0)

apval   = apval0
appval  = appval0
apppval = apppval0

case(1)

apval   = apval0
appval  = appval0 
apppval = apppval0 + 3*apval0

case(2)

apval   = apval0
appval  = appval0 - apval0
apppval  = apppval0 - 3*appval - 2*apval

end select


end subroutine



subroutine qprolates_convert_to0(isubst,apval,appval,apppval, &
  apval0,appval0,apppval0)
implicit double precision (a-h,o-z)
!
!  Given the values of the first three derivatives of a phase function \alpha(x)
!  for the differential (1) corresponding to one of the substitutions used in
!  qprolates, find the first three derivatives of the phase function \alpha0(t) 
!  for (1) corresponding to the substitution isubst = 0.
!
!  Input parameters:
!    isubst - integer parameter specifying the change of variables (see qprolates)
!    apval - the value of \alpha'(0)
!    appval - the value of \alpha''(0)
!    apppval - the value of \alpha'''(0)
!
!  Output parameters:  
!    apval0 - the value of alpha'(0)
!    appval0 - the value of alpha''(0)
!    apppval0 - the value of alpha'''(0)
!

select case(isubst)

case(0)

apval0   = apval
appval0  = appval
apppval0 = apppval

case(1)

apval0   = apval
appval0  = appval
apppval0 = apppval - 3 * apval 

case(2)

apval0   = apval
appval0  = appval + apval 
apppval0 = apppval + 3*appval + 2*apval

end select


end subroutine



subroutine qprolates_convert_to_w(apval,appval,apppval,wval,wpval,wppval)
implicit double precision (a-h,o-z)
!
!  Given the values of the first three derivatives of the alpha'(x),
!  compute the value of 
!
!    w(x) = 1 / alpha'(x)
!
!  and its first two derivatives.
!
!  Input parameters:
!    apval - the value of alpha'(x)
!    appval - the value of alpha''(x)
!    apppval - the value of alpha'''(x)
!
!  Output parameters:  
!    wval - the value of w(0)
!    wpval - the value of w'(0)
!    wppval - the value of w''(0)
!
!

wval   = 1/apval
wpval  = -appval/apval**2
wppval = 2*apval**2/apval**3 - apppval / apval**2

end subroutine



subroutine qprolates_convert_from_w(wval,wpval,wppval,apval,appval,apppval)
implicit double precision (a-h,o-z)
!
!  Given the values of 
!
!    w(x) = 1 / alpha'(x)
!
!  and its first two derivatives compute the first three derivatives of the 
!  alpha'(x).
!
!  Input parameters:
!    wval - the value of w(0)
!    wpval - the value of w'(0)
!    wppval - the value of w''(0)
!
!  Output parameters:  
!    apval - the value of alpha'(x)
!    appval - the value of alpha''(x)
!    apppval - the value of alpha'''(x)
!
!

apval   = 1 / wval
appval  = -wpval*apval**2
apppval = -wppval*apval**2 + 2*apval**2/apval

end subroutine



subroutine prolates_qfun(ifwindow,isubst,dlambda,c,chi,a,b,t,val,der)
implicit double precision (a-h,o-z)

!
!  Evaluate either the function q(x) or a windowed version of it.
!

call qprolates(isubst,c,chi,t,val,der)

if (ifwindow .eq. 0) return

qval = val
qder = der
call prolates_winfun(a,b,t,bval,bder)


! q(a) is equal to dlambda^2 
if (ifwindow .eq. 1) then
val = dlambda**2 * bval + (1-bval) * qval
der = dlambda**2 * bder  - bder * qval + (1-bval)*qder
endif

! q(b) is equal to dlambda^2 
if (ifwindow .eq. 2) then
val = dlambda**2 * (1-bval) + bval * qval
der = -dlambda**2 * bder  + bder * vqal + bval*qder
endif

end subroutine


subroutine prolates_ricode(k,ts,qs,userptr)
use iso_c_binding
implicit double precision (a-h,o-z)
integer                            :: k
double precision                   :: ts(k)
double complex                     :: qs(k)
type(c_ptr)                        :: userptr
type(prolates_qfun_data), pointer  :: qdata

!
!  External subroutine for solving the Riccati equation r'(t) + r(t)^2 + q(t) =0
!  via the solver in riccati.f90.
!

call c_f_pointer(userptr,qdata)

a        = qdata%a
b        = qdata%b
c        = qdata%c
chi      = qdata%chi
ifwindow = qdata%ifwindow
dlambda  = qdata%dlambda
isubst   = qdata%isubst


do i=1,k
t = ts(i)
call prolates_qfun(ifwindow,isubst,dlambda,c,chi,a,b,t,val,der)
qs(i) = val
end do

end subroutine


subroutine prolates_window_appell(eps,chebdata,c,chi,isubst,wa,wpa,wppa, &
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
!  where q(x) one of the coefficients in (1), on the interval (0,\infty) using an
!  adaptive Chebyshev spectral method.  The function
!
!    w(x) = 1 / alpha'(x)
!
!  satisfies the Appell equation.
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

maxints  = 10000000
k        = chebdata%k
nn       = k * dtail

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
b = 1d15

dd     = log(b-a)
nints0 = ceiling(dd)*2
dd     = dd/(nints0+0.0d0)

do i=1,nints0
int       = nints0-i+1
a0        = dd * (i-1)
b0        = dd * i
ab0(1,int) = a0
ab0(2,int) = b0
end do

ab0(:,1:nints0)        = a + exp(ab0(:,1:nints0))-1

nintsout = 0

do while (nints0 > 0) 

ifsplit = 0
a0      = ab0(1,nints0)
b0      = ab0(2,nints0)
nints0  = nints0 - 1

if (b0-a0 .eq. 0) then
ier = 256
call prin2("in solve3_ivp_adap interval of length 0, a = ",a)
return
endif

ts = (a0+b0)/2 + (b0-a0)/2 * chebdata%xs

do i=1,k

call qprolates(isubst,c,chi,ts(i),qval,qder)

rs(i) = 0
ps(i) = 4*qval
qs(i) = 2*qder
fs(i) = 0
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

call solve3_ivp_int(a0,b0,k,chebdata%xs,chebdata%aintl,chebdata%aintl2,chebdata%aintl3,  &
 rs,ps,qs,fs,ys0,yders0,yder2s0,yder3s0)

ifstop  = 0

do i=1,k

if (.not. ieee_is_finite(ys0(i))) then
ifsplit = 1
call prin2("ps = ",ps)
call prin2("qs = ",qs)

print *,a0,b0
stop
endif

if (abs(ys0(i)) .gt. 1d50) then
ifstop = 1
exit
endif

end do

if (ifstop .eq. 1) exit

ifsplit = 0

coefs0 = matmul(chebdata%u,ys0)
dd1 = maxval(abs(coefs0))
dd2 = maxval(abs(coefs0(k-nn+1:k)))
if (dd2 .gt. eps*dd1) ifsplit = 1

! coefs0 = matmul(chebdata%u,yders0)
! dd1 = maxval(abs(coefs0))
! dd2 = maxval(abs(coefs0(k-nn+1:k)))
! if (dd2 .gt. eps*dd1) ifsplit = 1

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




end module
