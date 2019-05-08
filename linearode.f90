!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This file contains code for solving second order linear ordinary differential equations
!  of the form
!
!    y''(t) + p(t) y'(t) + q(t) y(t) = f(t)                                                 (1)
!
!  and third order linear ordinary differential equations of the form
!
!    y'''(t) + r(t) y''(t) + p(t) y'(t) + q(t) y(t) = f(t)                                  (2)
!
!  via a straightforward adaptive spectral method.  The resulting solutions  are
!  represented via piecewise Chebyshev expansions of a specified order on a 
!  collection of subintervals determined adaptively.
!
!  The following subroutines  should be regarded as publicly callable:
!
!    solve2_ivp_adap - adaptively solve an initial value problem for (1)
!
!    solve2_tvp_adap - adaptively solve a terminal value problem for (1)
!
!    solve3_ivp_adap - adaptively solve an initial value problem for (2)
!
!    solve3_tvp_adap - adaptively solve a terminal value problem for (2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module linearode

use utils
use chebyshev
use iso_c_binding

integer, parameter, private          :: maxiterstrap = 8
double precision, parameter, private :: dtail = 0.5d0

interface

subroutine solve3fun(k,ts,rs,ps,qs,fs,userptr)
import c_ptr
implicit double precision (a-h,o-z)
integer                       :: k
double precision              :: ts(k),rs(k),ps(k),qs(k),fs(k)
type(c_ptr)                   :: userptr
!
!  Return the value of the coefficients r(t), p(t) and q(t) in (2).
!
end subroutine


subroutine solve2fun(k,ts,ps,qs,fs,userptr)
import c_ptr
implicit double precision (a-h,o-z)
integer                       :: k
double precision              :: ts(k),ps(k),qs(k),fs(k)
type(c_ptr)                   :: userptr
!
!  Return the value of the coefficients p(t) and q(t) in (1).
!
end subroutine

end interface

contains


subroutine solve2_ivp_adap(ier,eps,a,b,chebdata,nints,ab,ys,yders,yder2s, &
  odefun,ya,ypa,userptr)
implicit double precision (a-h,o-z)
type(chebexps_data)                        :: chebdata
procedure(solve2fun)                       :: odefun
integer, intent(out)                       :: nints
double precision, allocatable, intent(out) :: ab(:,:), ys(:,:),yders(:,:),yder2s(:,:)
double precision, intent(in)               :: ya,ypa
type(c_ptr)                                :: userptr

!
!  Solve an initial value problem for the ODE (1) via an adaptive Chebyshev
!  spectral method.  The solution and its derivatives are represented via
!  their values at 
!
!  Input parameters:
!    eps - the desired precision for the solution
!    (a,b) - the interval on which to solve (1)
!    chebdata - the structure returned by the chebexps routine -- this structure
!       determines the number of Chebyshev nodes used on each interval
!    odefun - a user-specified external conforming to the solve2fun interface
!     which supplies the values of the coefficients in (1)
!    ya - the initial value of the solution at a 
!    ypa - the initial value of the derivative of the solution at a
!    userptr - a "void *" pointer which is passed to odefun and can be
!      used to pass whatever data is desired to that subroutine
!
!  Output parameters:
!    ier - an error return code; 
!       ier =    0   indicates successful execution
!       ier =    4   the maximum number of intervals was exceeded while
!                    discretizing the coefficients p, q and r
!       ier =    8   the maximum number of intervals was exceeded while
!                    adaptively discretizing the solution
!       ier =   256  an interval of length 0 was encountered -- its endpoint
!                    will be displayed
!
!    nints - the number of intervals used to discretize the solution
!    ab - a (2,nints) array specifying the intervals used to discretize the
!      solution
!    ys - a (k,nints) array giving the values of the solutions at the k-point
!     Chebyshev
!    yders - 
!    yder2s - 
!

double precision, allocatable   :: ab0(:,:)               ! the list of intervals to process


! data for solving over a single interval
double precision, allocatable   :: ys0(:),yders0(:),yder2s0(:),coefs0(:)
double precision, allocatable   :: amatr(:,:),ts(:),ps(:),qs(:),fs(:)

! output data
double precision, allocatable   :: about(:,:)
double precision, allocatable   :: ysout(:,:),ydersout(:,:),yder2sout(:,:)


ier = 0

!
!  Set algorithm parameters and allocate memory for the procedure.
!

epssq    = eps**2
maxints  = 1000000
k        = chebdata%k
ntail    = k*dtail

allocate(ys0(k),yders0(k),yder2s0(k),coefs0(k),amatr(k,k))
allocate(ab0(2,maxints),qs(k),ps(k),ts(k),fs(k))
allocate(about(2,maxints))
allocate(ysout(k,maxints))
allocate(ydersout(k,maxints))
allocate(yder2sout(k,maxints))

!
!  First adaptive discretize the coefficients
!

nintsout = 0
nints0   = 1

ab0(1,1) = a
ab0(2,1) = b

do while (nints0 > 0)

a0 = ab0(1,nints0)
b0 = ab0(2,nints0)
nints0 = nints0-1

if (b0-a0 .eq. 0) then
call prin2("in solve2_ivp_adap interval of length 0, a = ",a)
ier = 256
return
endif

ts = (b0-a0)/2 * chebdata%xs + (b0+a0)/2
call odefun(k,ts,ps,qs,fs,userptr)

ifsplit = 0

coefs0 = matmul(chebdata%u,ps)
dd3    = sum(coefs0**2)
dd4    = sum(coefs0(k-ntail+1:k)**2)

if (dd4 > dd3 * epssq) then
ifsplit = 1
goto 1000
endif

coefs0 = matmul(chebdata%u,qs)
dd5    = sum(coefs0**2)
dd6    = sum(coefs0(k-ntail+1:k)**2)

if (dd6 > dd5 * epssq) then
ifsplit = 1
goto 1000
endif

1000 continue

if (ifsplit .eq. 0) then

if (nintsout + 1 .gt. maxints) then
ier = 4
return
endif

nintsout          = nintsout +1 
about(1,nintsout) = a0
about(2,nintsout) = b0

else

if (nints0+2 .gt. maxints) then
ier = 4
return
endif


nints0 = nints0 + 1
ab0(1,nints0) = (a0+b0)/2
ab0(2,nints0) = b0

nints0 = nints0 + 1
ab0(1,nints0) = a0
ab0(2,nints0) = (a0+b0)/2

endif

end do


!
!  Initialize the list of intervals which need to be processed.
!

nints0   = nintsout
do i     = 1,nints0
int      = nints0-i+1
ab0(:,i) = about(:,int)
enddo

nintsout        = 0


do while (nints0 > 0) 

ifsplit = 0
a0      = ab0(1,nints0)
b0      = ab0(2,nints0)
nints0  = nints0 - 1

if (b0-a0 .eq. 0) then
ier = 256
call prin2("in solve2_ivp_adap interval of length 0, a = ",a)
return
endif

ts = (a0+b0)/2 + (b0-a0)/2 * chebdata%xs
call odefun(k,ts,ps,qs,fs,userptr)

if (nintsout .eq. 0) then
ys0(1)     = ya
yders0(1)  = ypa
else
ys0(1)     = ysout(k,nintsout)
yders0(1)  = ydersout(k,nintsout)
endif

call solve2_ivp_int(a0,b0,k,chebdata%xs,chebdata%aintl,chebdata%aintl2,  &
 ps,qs,fs,ys0,yders0,yder2s0)

coefs0 = matmul(chebdata%u,ys0)
coefs0 = coefs0/coefs0(1)

dd1 = sum(coefs0(1:k)**2)
dd2 = sum(coefs0(k-ntail+1:k)**2)

if (dd2 .gt. epssq*dd1) then

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

if (nintsout+2 .gt. maxints) then
ier = 8
return
endif

nintsout          = nintsout+1
about(1,nintsout) = a0
about(2,nintsout) = b0

ysout(:,nintsout)     = ys0
ydersout(:,nintsout ) = yders0
yder2sout(:,nintsout) = yder2s0

endif


end do


!
!  Copy out the data
!

nints = nintsout
allocate(ab(2,nints))
allocate(ys(k,nints))
allocate(yders(k,nints))
allocate(yder2s(k,nints))

ab = about(:,1:nints)
ys = ysout(:,1:nints)
yders = ydersout(:,1:nints)
yder2s = yder2sout(:,1:nints)

end subroutine



subroutine solve2_ivp_int(a,b,k,xscheb,chebint,chebint2,ps,qs,fs,ys,yders,yder2s)
implicit double precision (a-h,o-z)

integer, intent(in)           :: k
double precision, intent(in)  :: xscheb(k),chebint(k,k),chebint2(k,k)
double precision, intent(in)  :: ps(k),qs(k),fs(k)
double precision, intent(out) :: ys(k),yders(k),yder2s(k)

!
!  Solve an initial value for the ordinary differential equation
!
!     y''(t) + p(t) y'(t) + q(t) y(t) = f(t)                                                (3)
!
!  on the interval [a,b] using a standard Chebyshev spectral method.
!
!  Input parameters:
!
!    (a,b) - the interval on which the ODE (3) is given
!    k - the number of Chebyshev points on the interval [a,b]
!    xscheb - the nodes of the k-point Chebyshev grid on [-1,1] as returned by
!     the subroutine chebexps (see below)
!    chebint? - the "left" Chebyshev spectral integration matrice as constructed
!      by the subroutine chebexps 
!    ps - an array specifying the values of the function p(t) appearing in (3)
!      at the k Chebyshev nodes on [a,b]
!    qs - an array specifying the values of the function q(t) appearing in (3)
!      at the k Chebyshev node on [a,b]
!    fs - an array speciying the values of the function f(t) appearing in (3)
!      at the k Chebyshev nodes on [a,b]
!
!    ys(1) - the value of y(a)
!    yders(1) - the value of y'(a)
!
!  Output parameters:
!
!    ys - the values of the solution y of (3) at the nodes of the k-point
!      Chebyshev grid on [a,b]
!    yders - the values of the solution y' of (3) at the nodes of the k-point
!      Chebyshev grid on [a,b]
!    yder2s - the values of the solution y'' of (3) at the nodes of the k-point
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



!
!  We represent the solution in the form
!
!      y(t) = alpha + beta (t-a) + \int_a^t (t-s) sigma(s) ds,
!
!  insert this representation into (3), and solve the resulting system of
!  linear equations in order to obtain the values of sigma.  
!    

amatr  = 0
do i=1,k
amatr(i,i) = 1.0d0
sigma(i)   = 0.0d0
end do

!     y''(t) + p(t) y'(t) + q(t) y(t) = f(t)                                                (3)


!
! Handle the p(t) * y'(t) term.
!
do i=1,k
amatr(i,:) = amatr(i,:) + ps(i) * chebint(i,:)*((b-a)/2)
sigma(i)   = sigma(i) - ps(i)*beta
end do

!
!  Handle the q(t) y(t) term
!

do i=1,k
amatr(i,:) = amatr(i,:) + qs(i) * chebint2(i,:)*((b-a)/2)**2
sigma(i)   = sigma(i) - qs(i) * (alpha + beta*(xs(i)-a) )
end do


!
!  Form the right-hand side.
!
do i=1,k
sigma(i) = sigma(i) + fs(i)
end do

!
!  Use a QR decomposition to invert the linear system
!

call qrsolv(amatr,k,sigma)

!
!  Calculate y(t), y'(t) and y''(t) from sigma.
!

yder2s = sigma
yders  = (b-a)/2*matmul(chebint,sigma)
ys     = ((b-a)/2)**2*matmul(chebint2,sigma)

do i=1,k
ys(i)      = ys(i) + alpha + beta*(xs(i)-a) 
yders(i)   = yders(i) + beta
end do

end subroutine


subroutine solve2_tvp_adap(ier,eps,a,b,chebdata,nints,ab,ys,yders,yder2s, &
  odefun,yb,ypb,userptr)
implicit double precision (a-h,o-z)
type(chebexps_data)                        :: chebdata
procedure(solve2fun)                       :: odefun
integer, intent(out)                       :: nints
double precision, allocatable, intent(out) :: ab(:,:), ys(:,:),yders(:,:),yder2s(:,:)
double precision, intent(in)               :: yb,ypb
type(c_ptr)                                :: userptr

!
!  Solve a terminal value problem for the ODE (1) via an adaptive Chebyshev
!  spectral method.
!
!  Input parameters:
!    eps - the desired precision for the calculatins
!    (a,b) - the interval on which to solve (1)
!    chebdata - the structure returned by the chebexps routine 
!    odefun - 
!    ya -
!    ypa -
!    userptr - 
!  
!
!  Output parameters:
!    ier - an error return code; 
!       ier =    0   indicates successful execution
!       ier =    4   the maximum number of intervals was exceeded while
!                    discretizing the coefficients p, q and r
!       ier =    8   the maximum number of intervals was exceeded while
!                    adaptively discretizing the solution
!       ier =   256  an interval of length 0 was encountered -- its endpoint
!                    will be displayed
!
!    nints - the number of intervals used to discretize the solution
!    ab - a (2,nints) array specifying the intervals used to discretize the
!      solution
!    ys -
!    yders - 
!    yder2s - 
!

double precision, allocatable   :: ab0(:,:)               ! the list of intervals to process


! data for solving over a single interval
double precision, allocatable   :: ys0(:),yders0(:),yder2s0(:),coefs0(:)
double precision, allocatable   :: amatr(:,:),ts(:),ps(:),qs(:),fs(:)

! output data
double precision, allocatable   :: about(:,:)
double precision, allocatable   :: ysout(:,:),ydersout(:,:),yder2sout(:,:)

ier = 0

!
!  Set algorithm parameters and allocate memory for the procedure.
!

epssq    = eps**2
maxints  = 1000000
k        = chebdata%k
ntail    = k*dtail

allocate(ys0(k),yders0(k),yder2s0(k),coefs0(k),amatr(k,k))
allocate(ab0(2,maxints),qs(k),ps(k),ts(k),fs(k))
allocate(about(2,maxints))
allocate(ysout(k,maxints))
allocate(ydersout(k,maxints))
allocate(yder2sout(k,maxints))

!
!  First adaptive discretize the coefficients
!

nintsout = 0
nints0   = 1

ab0(1,1) = a
ab0(2,1) = b

do while (nints0 > 0)

a0 = ab0(1,nints0)
b0 = ab0(2,nints0)
nints0 = nints0-1


if (b0-a0 .eq. 0) then
call prin2("in solve2_tvp_adap interval of length 0, a = ",a)
ier = 256
return
endif

ts = (b0-a0)/2 * chebdata%xs + (b0+a0)/2
call odefun(k,ts,ps,qs,fs,userptr)

ifsplit = 0

coefs0 = matmul(chebdata%u,ps)
dd3    = sum(coefs0**2)
dd4    = sum(coefs0(k-ntail+1:k)**2)

if (dd4 > dd3 * epssq) then
ifsplit = 1
goto 1000
endif

coefs0 = matmul(chebdata%u,qs)
dd5    = sum(coefs0**2)
dd6    = sum(coefs0(k-ntail+1:k)**2)

if (dd6 > dd5 * epssq) then
ifsplit = 1
goto 1000
endif


1000 continue


if (ifsplit .eq. 0) then

if (nintsout + 1 .gt. maxints) then
ier = 4
return
endif

nintsout          = nintsout +1 
about(1,nintsout) = a0
about(2,nintsout) = b0

else

if (nints0+2 .gt. maxints) then
ier = 4
return
endif


nints0 = nints0 + 1
ab0(1,nints0) = a0
ab0(2,nints0) = (a0+b0)/2

nints0 = nints0 + 1
ab0(1,nints0) = (a0+b0)/2
ab0(2,nints0) = b0


endif

end do



!
!  Initialize the list of intervals which need to be processed.
!

nints0   = nintsout
do i     = 1,nints0
int      = nints0-i+1
ab0(:,i) = about(:,int)
enddo

nintsout        = 0



do while (nints0 > 0) 

ifsplit = 0
a0      = ab0(1,nints0)
b0      = ab0(2,nints0)
nints0  = nints0 - 1


if (b0-a0 .eq. 0) then
ier = 256
call prin2("in solve2_tvp_adap interval of length 0, a = ",a)
return
endif

ts = (a0+b0)/2 + (b0-a0)/2 * chebdata%xs
call odefun(k,ts,ps,qs,fs,userptr)

if (nintsout .eq. 0) then
ys0(k)     = yb
yders0(k)  = ypb
else
ys0(k)     = ysout(1,nintsout)
yders0(k)  = ydersout(1,nintsout)
endif

call solve2_tvp_int(a0,b0,k,chebdata%xs,chebdata%aintr,chebdata%aintr2,  &
 ps,qs,fs,ys0,yders0,yder2s0)

coefs0 = matmul(chebdata%u,ys0)
coefs0 = coefs0/coefs0(1)

dd1 = sum(coefs0(1:k)**2)
dd2 = sum(coefs0(k-ntail+1:k)**2)

if (dd2 .gt. epssq*dd1) then

if (nints0+2 .gt. maxints) then
ier = 8
return
endif


nints0 = nints0+1
ab0(1,nints0) = a0
ab0(2,nints0) = (a0+b0)/2

nints0 = nints0+1
ab0(1,nints0) = (a0+b0)/2
ab0(2,nints0) = b0



else

if (nintsout+2 .gt. maxints) then
ier = 8
return
endif

nintsout          = nintsout+1
about(1,nintsout) = a0
about(2,nintsout) = b0

ysout(:,nintsout)     = ys0
ydersout(:,nintsout ) = yders0
yder2sout(:,nintsout) = yder2s0

endif


end do


!
!  Copy out the data
!

nints = nintsout
allocate(ab(2,nints))
allocate(ys(k,nints))
allocate(yders(k,nints))
allocate(yder2s(k,nints))

do i=1,nints
int = nints-i+1
ab(:,i)= about(:,int)
ys(:,i) = ysout(:,int)
yders(:,i) = ydersout(:,int)
yder2s(:,i) = yder2sout(:,int)
end do

end subroutine



subroutine solve2_tvp_int(a,b,k,xscheb,chebint,chebint2,ps,qs,fs,ys,yders,yder2s)
implicit double precision (a-h,o-z)

integer, intent(in)           :: k
double precision, intent(in)  :: xscheb(k),chebint(k,k),chebint2(k,k)
double precision, intent(in)  :: ps(k),qs(k),fs(k)
double precision, intent(out) :: ys(k),yders(k),yder2s(k)

!
!  Solve a terminal value for the ordinary differential equation
!
!     y''(t) + p(t) y'(t) + q(t) y(t) = f(t)                                                (4)
!
!  on the interval [a,b] using a standard Chebyshev spectral method.
!
!  Input parameters:
!
!    (a,b) - the interval on which the ODE (4) is given
!    k - the number of Chebyshev points on the interval [a,b]
!    xscheb - the nodes of the k-point Chebyshev grid on [-1,1] as returned by
!     the subroutine chebexps (see below)
!    chebint? - the "left" Chebyshev spectral integration matrice as constructed
!      by the subroutine chebexps 
!    ps - an array specifying the values of the function p(t) appearing in (4)
!      at the k Chebyshev nodes on [a,b]
!    qs - an array specifying the values of the function q(t) appearing in (4)
!      at the k Chebyshev node on [a,b]
!    fs - an array speciying the values of the function f(t) appearing in (4)
!      at the k Chebyshev nodes on [a,b]
!
!    ys(k) - the value of y(b)
!    yders(k) - the value of y'(b)
!
!  Output parameters:
!
!    ys - the values of the solution y of (4) at the nodes of the k-point
!      Chebyshev grid on [a,b]
!    yders - the values of the solution y' of (4) at the nodes of the k-point
!      Chebyshev grid on [a,b]
!    yder2s - the values of the solution y'' of (4) at the nodes of the k-point
!      Chebyshev grid on [a,b]

!

double precision, allocatable :: amatr(:,:),xs(:),sigma(:),rhs(:)

!
!  Allocate memory for the procedure and setup some parameters.
!


allocate(amatr(k,k),xs(k),sigma(k),rhs(k))


xs       = (b-a)/2 *xscheb + (b+a)/2 
alpha    = ys(k)
beta     = yders(k)


!
!  We represent the solution in the form
!
!      y(t) = alpha + beta (t-b) + \int_b^t (t-s) sigma(s) ds,
!
!  insert this representation into (3), and solve the resulting system of
!  linear equations in order to obtain the values of sigma.  
!    

amatr  = 0
do i=1,k
amatr(i,i) = 1.0d0
sigma(i)   = 0.0d0
end do

!     y''(t) + p(t) y'(t) + q(t) y(t) = f(t)                                                (3)


!
! Handle the p(t) * y'(t) term.
!
do i=1,k
amatr(i,:) = amatr(i,:) + ps(i) * chebint(i,:)*((b-a)/2)
sigma(i)   = sigma(i) - ps(i)*beta
end do

!
!  Handle the q(t) y(t) term
!

do i=1,k
amatr(i,:) = amatr(i,:) + qs(i) * chebint2(i,:)*((b-a)/2)**2
sigma(i)   = sigma(i) - qs(i) * (alpha + beta*(xs(i)-b) )
end do


!
!  Form the right-hand side.
!
do i=1,k
sigma(i) = sigma(i) + fs(i)
end do


!
!  Use a QR decomposition to invert the linear system
!

call qrsolv(amatr,k,sigma)


!
!  Calculate y(t), y'(t) and y''(t) from sigma.
!

yder2s = sigma
yders  = (b-a)/2*matmul(chebint,sigma)
ys     = ((b-a)/2)**2*matmul(chebint2,sigma)




do i=1,k
ys(i)      = ys(i) + alpha + beta*(xs(i)-b) 
yders(i)   = yders(i) + beta
end do




end subroutine




subroutine solve3_ivp_adap(ier,eps,a,b,chebdata,nints,ab,ys,yders,yder2s,yder3s,&
  odefun,ya,ypa,yppa,userptr)
implicit double precision (a-h,o-z)
type(chebexps_data)                        :: chebdata
procedure(solve3fun)                       :: odefun
integer, intent(out)                       :: nints
double precision, allocatable, intent(out) :: ab(:,:), ys(:,:),yders(:,:),yder2s(:,:),yder3s(:,:)
double precision, intent(in)               :: ya,ypa,yppa
type(c_ptr)                                :: userptr

!
!  Solve an initial value problem for the ODE (2) via an adaptive Chebyshev
!  spectral method.
!
!  Input parameters:
!    eps - the desired precision for the calculatins
!    (a,b) - the interval on which to solve (1)
!    chebdata - the structure returned by the chebexps routine 
!    odefun - 
!    ya -
!    ypa -
!    yppa - 
!    userptr - 
!  
!
!  Output parameters:
!    ier - an error return code; 
!       ier =    0   indicates successful execution
!       ier =    4   the maximum number of intervals was exceeded while
!                    discretizing the coefficients p, q and r
!       ier =    8   the maximum number of intervals was exceeded while
!                    adaptively discretizing the solution
!       ier =   256  an interval of length 0 was encountered -- its endpoint
!                    will be displayed
!
!    nints - the number of intervals used to discretize the solution
!    ab - a (2,nints) array specifying the intervals used to discretize the
!      solution
!    ys -
!    yders - 
!    yder2s - 
!    yder3s - 
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

epssq    = eps**2
maxints  = 10000000
k        = chebdata%k
ntail    = k*dtail

allocate(ys0(k),yders0(k),yder2s0(k),yder3s0(k),coefs0(k),amatr(k,k))
allocate(ab0(2,maxints),rs(k),qs(k),ps(k),ts(k),fs(k)) 
allocate(about(2,maxints))
allocate(ysout(k,maxints))
allocate(ydersout(k,maxints))
allocate(yder2sout(k,maxints))
allocate(yder3sout(k,maxints))

!
!  First adaptive discretize the coefficients 
!

nints0   = 1
ab0(1,1) = a
ab0(2,1) = b


nintsout = 0


do while (nints0 > 0)

a0 = ab0(1,nints0)
b0 = ab0(2,nints0)
nints0 = nints0-1

if (b0-a0 .eq. 0) then
call prin2("in solve3_ivp_adap interval of length 0, a = ",a)
ier = 256
return
endif

ts = (b0-a0)/2 * chebdata%xs + (b0+a0)/2
call odefun(k,ts,rs,ps,qs,fs,userptr)

ifsplit = 0

coefs0 = matmul(chebdata%u,rs)
dd1    = sum(coefs0**2)
dd2    = sum(coefs0(k-ntail+1:k)**2)

if (dd2 > dd1 * epssq) then
ifsplit = 1
goto 1000
endif

coefs0 = matmul(chebdata%u,ps)
dd3    = sum(coefs0**2)
dd4    = sum(coefs0(k-ntail+1:k)**2)

if (dd4 > dd3 * epssq) then
ifsplit = 1
goto 1000
endif

coefs0 = matmul(chebdata%u,qs)
dd5    = sum(coefs0**2)
dd6    = sum(coefs0(k-ntail+1:k)**2)

if (dd6 > dd5 * epssq) then
ifsplit = 1
goto 1000
endif

1000 continue

if (ifsplit .eq. 0) then

if (nintsout + 1 .gt. maxints) then
ier = 4
return
endif

nintsout          = nintsout +1 
about(1,nintsout) = a0
about(2,nintsout) = b0

else

if (nints0+2 .gt. maxints) then
ier = 4
return
endif


nints0 = nints0 + 1
ab0(1,nints0) = (a0+b0)/2
ab0(2,nints0) = b0

nints0 = nints0 + 1
ab0(1,nints0) = a0
ab0(2,nints0) = (a0+b0)/2

endif

end do


!
!  Initialize the list of intervals which need to be processed.
!

nints0   = nintsout
do i     = 1,nints0
int      = nints0-i+1
ab0(:,i) = about(:,int)
enddo

nintsout        = 0


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
call odefun(k,ts,rs,ps,qs,fs,userptr)

if (nintsout .eq. 0) then
ys0(1)     = ya
yders0(1)  = ypa
yder2s0(1) = yppa
else
ys0(1)     = ysout(k,nintsout)
yders0(1)  = ydersout(k,nintsout)
yder2s0(1) = yder2sout(k,nintsout)
endif

ya0   = ya
ypa0  = ypa
yppa0 = yppa

call solve3_ivp_int(a0,b0,k,chebdata%xs,chebdata%aintl,chebdata%aintl2,chebdata%aintl3,  &
 rs,ps,qs,fs,ys0,yders0,yder2s0,yder3s0)

! do i=1,k

! if (ys0(i)+1 .eq. ys0(i)) then

! call prina("in solve3_ivp_adap, NaN or Inf")
! call prin2("a0   = ",a0)
! call prin2("b0   = ",b0)
! call prin2("rs   = ",rs)
! call prin2("ps   = ",ps)
! call prin2("qs   = ",qs)
! call prin2("fs   = ",fs)
! call prin2("y0   = ",ya0)
! call prin2("yp0  = ",ypa0)
! call prin2("ypp0 = ",yppa0)
! call prin2("ys0 = ",ys0)

! endif

! end do


coefs0 = matmul(chebdata%u,ys0)
coefs0 = coefs0/coefs0(1)

dd1 = sum(coefs0(1:k)**2)
dd2 = sum(coefs0(k/2+1:k)**2)


if (dd2 .gt. epssq*dd1) then

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
allocate(ys(k,nints))
allocate(yders(k,nints))
allocate(yder2s(k,nints))
allocate(yder3s(k,nints))

ab = about(:,1:nints)
ys = ysout(:,1:nints)
yders = ydersout(:,1:nints)
yder2s = yder2sout(:,1:nints)
yder3s = yder3sout(:,1:nints)

end subroutine




subroutine solve3_ivp_int(a,b,k,xscheb,chebint,chebint2,chebint3,rs,ps,qs,fs,ys,yders,yder2s,yder3s)
implicit double precision (a-h,o-z)

integer, intent(in)           :: k
double precision, intent(in)  :: xscheb(k),chebint(k,k),chebint2(k,k),chebint3(k,k)
double precision, intent(in)  :: rs(k),ps(k),qs(k),fs(k)
double precision, intent(out) :: ys(k),yders(k),yder2s(k),yder3s(k)

!
!  Solve an initial value for the ordinary differential equation
!
!     y'''(t) + r(t) y''(t) + p(t) y'(t) + q(t) y(t) = f(t)                                      (5)
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
!    rs - an array specifying the values of the function r(t) appearing in (5)
!      at the k Chebyshev nodes on [a,b]
!    ps - an array specifying the values of the function p(t) appearing in (5)
!      at the k Chebyshev nodes on [a,b]
!    qs - an array specifying the values of the function q(t) appearing in (5)
!      at the k Chebyshev node on [a,b]
!    fs - an array speciying the values of the function f(t) appearing in (5)
!      at the k Chebyshev nodes on [a,b]
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
do i=1,k
amatr(i,:) = amatr(i,:) + rs(i) * chebint(i,:)*(b-a)/2
sigma(i)   = sigma(i) - rs(i)*eta
end do

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


!
!  Form the right-hand side.
!
do i=1,k
sigma(i) = sigma(i) + fs(i)
end do

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




subroutine solve3_tvp_adap(ier,eps,a,b,chebdata,nints,ab,ys,yders,yder2s,yder3s,&
  odefun,yb,ypb,yppb,userptr)
implicit double precision (a-h,o-z)
type(chebexps_data)                        :: chebdata
procedure(solve3fun)                       :: odefun
integer, intent(out)                       :: nints
double precision, allocatable, intent(out) :: ab(:,:), ys(:,:),yders(:,:),yder2s(:,:),yder3s(:,:)
double precision, intent(in)               :: yb,ypb,yppb
type(c_ptr)                                :: userptr

!
!  Solve a terminal value problem for the ODE (2) via an adaptive Chebyshev
!  spectral method.
!
!  Input parameters:
!    eps - the desired precision for the calculatins
!    (a,b) - the interval on which to solve (1)
!    chebdata - the structure returned by the chebexps routine 
!    odefun - 
!    yb -
!    ypb -
!    yppb - 
!    userptr - 
!  
!
!  Output parameters:
!    ier - an error return code; 
!       ier =    0   indicates successful execution
!       ier =    4   the maximum number of intervals was exceeded while
!                    discretizing the coefficients p, q and r
!       ier =    8   the maximum number of intervals was exceeded while
!                    adaptively discretizing the solution
!       ier =   256  an interval of length 0 was encountered -- its endpoint
!                    will be displayed
!
!    nints - the number of intervals used to discretize the solution
!    ab - a (2,nints) array specifying the intervals used to discretize the
!      solution
!    ys -
!    yders - 
!    yder2s - 
!    yder3s - 
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

epssq    = eps**2
maxints  = 1000000
k        = chebdata%k
ntail    = k*dtail

allocate(ys0(k),yders0(k),yder2s0(k),yder3s0(k),coefs0(k),amatr(k,k))
allocate(ab0(2,maxints),rs(k),qs(k),ps(k),ts(k),fs(k))
allocate(about(2,maxints))
allocate(ysout(k,maxints))
allocate(ydersout(k,maxints))
allocate(yder2sout(k,maxints))
allocate(yder3sout(k,maxints))

!
!  First adaptive discretize the coefficients
!

nintsout = 0
nints0   = 1

ab0(1,1) = a
ab0(2,1) = b

do while (nints0 > 0)

a0 = ab0(1,nints0)
b0 = ab0(2,nints0)
nints0 = nints0-1

if (b0-a0 .eq. 0) then
call prin2("in solve3_tvp_adap interval of length 0, a = ",a)
ier = 256
return
endif

ts = (b0-a0)/2 * chebdata%xs + (b0+a0)/2
call odefun(k,ts,rs,ps,qs,fs,userptr)

ifsplit = 0

coefs0 = matmul(chebdata%u,rs)
dd1    = sum(coefs0**2)
dd2    = sum(coefs0(k-ntail+1:k)**2)


if (dd2 > dd1 * epssq) then
ifsplit = 1
goto 1000
endif


coefs0 = matmul(chebdata%u,ps)
dd3    = sum(coefs0**2)
dd4    = sum(coefs0(k-ntail+1:k)**2)

if (dd4 > dd3 * epssq) then
ifsplit = 1
goto 1000
endif

coefs0 = matmul(chebdata%u,qs)
dd5    = sum(coefs0**2)
dd6    = sum(coefs0(k-ntail+1:k)**2)

if (dd6 > dd5 * epssq) then
ifsplit = 1
goto 1000
endif

1000 continue

if (ifsplit .eq. 0) then

if (nintsout + 1 .gt. maxints) then
ier = 4
return
endif

nintsout          = nintsout +1 
about(1,nintsout) = a0
about(2,nintsout) = b0

else

if (nints0+2 .gt. maxints) then
ier = 4
return
endif


nints0 = nints0 + 1
ab0(1,nints0) = (a0+b0)/2
ab0(2,nints0) = b0

nints0 = nints0 + 1
ab0(1,nints0) = a0
ab0(2,nints0) = (a0+b0)/2

endif

end do


!
!  Initialize the list of intervals which need to be processed.
!

nints0   = nintsout
do i     = 1,nints0
int      = i
ab0(:,i) = about(:,int)
enddo
nintsout        = 0


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
call odefun(k,ts,rs,ps,qs,fs,userptr)

if (nintsout .eq. 0) then
ys0(k)     = yb
yders0(k)  = ypb
yder2s0(k) = yppb
else
ys0(k)     = ysout(1,nintsout)
yders0(k)  = ydersout(1,nintsout)
yder2s0(k) = yder2sout(1,nintsout)
endif


call solve3_tvp_int(a0,b0,k,chebdata%xs,chebdata%aintr,chebdata%aintr2,chebdata%aintr3,  &
 rs,ps,qs,fs,ys0,yders0,yder2s0,yder3s0)


coefs0 = matmul(chebdata%u,ys0)
coefs0 = coefs0/coefs0(1)

dd1 = sum(coefs0(1:k)**2)
dd2 = sum(coefs0(k-ntail+1:k)**2)

if (dd2 .gt. epssq*dd1) then

if (nints0+2 .gt. maxints) then
ier = 8
return
endif

nints0 = nints0+1
ab0(1,nints0) = a0
ab0(2,nints0) = (a0+b0)/2

nints0 = nints0+1
ab0(1,nints0) = (a0+b0)/2
ab0(2,nints0) = b0

else

if (nintsout+2 .gt. maxints) then
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
allocate(ys(k,nints))
allocate(yders(k,nints))
allocate(yder2s(k,nints))
allocate(yder3s(k,nints))

do i=1,nints
int = nints-i+1
ab(:,i)     = about(:,int)
ys(:,i)     = ysout(:,int)
yders(:,i)  = ydersout(:,int)
yder2s(:,i) = yder2sout(:,int)
yder3s(:,i) = yder3sout(:,int)
end do

end subroutine


subroutine solve3_tvp_int(a,b,k,xscheb,chebint,chebint2,chebint3,rs,ps,qs,fs,ys,yders,yder2s,yder3s)
implicit double precision (a-h,o-z)

integer, intent(in)           :: k
double precision, intent(in)  :: xscheb(k),chebint(k,k),chebint2(k,k),chebint3(k,k)
double precision, intent(in)  :: rs(k),ps(k),qs(k),fs(k)
double precision, intent(out) :: ys(k),yders(k),yder2s(k),yder3s(k)

!
!  Solve a terminal value for the ordinary differential equation
!
!     y'''(t) + r(t) y''(t) + p(t) y'(t) + q(t) y(t) = f(t)                                      (6)
!
!  on the interval [a,b] using a standard Chebyshev spectral method.
!
!  Input parameters:
!
!    (a,b) - the interval on which the ODE (6) is given
!    k - the number of Chebyshev points on the interval [a,b]
!    xscheb - the nodes of the k-point Chebyshev grid on [-1,1] as returned by
!     the subroutine chebexps (see below)
!    chebint? - the "right" Chebyshev spectral integration matrices as constructed
!      by the subroutine chebexps 
!    rs - an array specifying the values of the function r(t) appearing in (6)
!      at the k Chebyshev nodes on [a,b]
!    ps - an array specifying the values of the function p(t) appearing in (6)
!      at the k Chebyshev nodes on [a,b]
!    qs - an array specifying the values of the function q(t) appearing in (6)
!      at the k Chebyshev node on [a,b]
!    fs - an array speciying the values of the function f(t) appearing in (6)
!      at the k Chebyshev nodes on [a,b]
!
!    ys(k) - the value of y(b)
!    yders(k) - the value of y'(b)
!    yder2s(k) - the value of y''(b)
!
!  Output parameters:
!
!    ys - the values of the solution y of (6) at the nodes of the k-point
!      Chebyshev grid on [a,b]
!    yders - the values of the solution y' of (6) at the nodes of the k-point
!      Chebyshev grid on [a,b]
!    yder2s - the values of the solution y'' of (6) at the nodes of the k-point
!      Chebyshev grid on [a,b]
!    yder3s - the values of the solution y''' of (6) at the nodes of the k-point
!      Chebyshev grid on [a,b]

!

double precision, allocatable :: amatr(:,:),xs(:),sigma(:),rhs(:)

!
!  Allocate memory for the procedure and setup some parameters.
!

allocate(amatr(k,k),xs(k),sigma(k),rhs(k))


xs       = (b-a)/2 *xscheb + (b+a)/2 
alpha    = ys(k)
beta     = yders(k)
eta      = yder2s(k)

!
!  We represent the solution in the form
!
!      y(t) = alpha + beta (t-b) + eta(t-b)^2/2 + \int_b^t (t-s)^2/2 sigma(s) ds,
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
do i=1,k
amatr(i,:) = amatr(i,:) + rs(i) * chebint(i,:)*(b-a)/2
sigma(i)   = sigma(i) - rs(i)*eta
end do

!
! Handle the p(t) * y'(t) term.
!
do i=1,k
amatr(i,:) = amatr(i,:) + ps(i) * chebint2(i,:)*((b-a)/2)**2
sigma(i)   = sigma(i) - ps(i)*(beta + eta * (xs(i)-b))
end do

!
!  Handle the q(t) y(t) term
!

do i=1,k
amatr(i,:) = amatr(i,:) + qs(i) * chebint3(i,:)*((b-a)/2)**3
sigma(i)   = sigma(i) - qs(i) * (alpha + beta*(xs(i)-b) + eta * (xs(i)-b)**2/2)
end do


!
!  Form the right-hand side.
!
do i=1,k
sigma(i) = sigma(i) + fs(i)
end do

!
!  Use a QR decomposition to invert the linear system
!

call qrsolv(amatr,k,sigma)

!
!  Calculate y(t), y'(t) and y''(t) from sigma.
!

yder3s = sigma
yder2s = (b-a)/2*matmul(chebint,sigma)
yders  = ((b-a)/2)**2*matmul(chebint2,sigma)
ys     = ((b-a)/2)**3*matmul(chebint3,sigma)

do i=1,k
ys(i)      = ys(i) + alpha + beta*(xs(i)-b) + eta * (xs(i)-b)**2/2
yders(i)   = yders(i) + beta   + eta * (xs(i)-b)
yder2s(i)  = yder2s(i) + eta
end do

end subroutine


end module
