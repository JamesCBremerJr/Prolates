!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for constructing a nonoscillatory phase function
!  for the differential equation 
!
!     z''(x) + q(x) z(x) = 0                                                              (1)
!
!  obtained from  the spheroidal wave equation 
!
!    (1-t^2) * y''(t) - 2x y'(t) + ( chi - dmu^2/ (1-t^2)  - gamma^2 t^2 ) y(t) = 0       (2)
!
!  via the transformation
!
!    z(x) = y(1-exp(-x)) sqrt(2-exp(-x)).
!
!  The solutions of (2) are typically indexed by dmu, gamma and by an implicit
!  parameter dnu, called the characteristic exponent.  The relationship between
!  chi and the parameters gamma, dmu and dnu is fairly complicated (for instance,
!  chi has branch points at each half integer value of dnu) and the code here
!  provide no means of computing the desired value of chi as a function of
!  them.
!
!  See, however, prolates_exps.f90 which provides a mechanism for computing the values
!  of chi corresponding to the prolate spheroidal wave functions of order 0 and
!  integer characteristic exponents.
!
!  This code operates under the assumptions that gamma^2 >= 0, dmu >= 0 and
!  that the value of chi is chosen so that the corresponding value of dnu
!  is greater than or equal to dmu.
!
!  The phase function Psi is determined up to an additive constant by the fact
!  that it corresponds to a solution of (2) which is in the Hardy space H^2 of 
!  functions analytic on the upper half of the complex plane whose boundary values are
!  square integrable.  We fix the constant through the requirement that
!
!     \lim   Psi(x) = 0.
!    x \to 1-
!
!  This ensures that 
!
!                    sin (  Psi(t) )  
!    u(x) =       -------------------                                                     (3)
!                   sqrt(  Psi'(t) )
!
!  is related to angular prolate spheroidal wave function of the first kind via
!  the formula
!
!      dmu                       dmu             u(-log(1-t))
!    ps    (t; gamma^2) =      C   (gamma^2)  ----------------                            (4)
!      dnu                       dnu               sqrt(1+t)
!
!  with C_dnu^dmu(\gamma^2) a normalization constant.  In general, it is not the 
!  case that the angular function of the second kind is a multiple of
!
!                    cos (  Psi(t) )  
!    v(x) =       -------------------                                                     (5)
!                   sqrt{  Psi'(t) } 
!
!  (as one might expect).  However, when dnu and dmu are integers it is the
!  case that
!
!      m                         m             v(-log(1-t))
!    qs    (t; gamma^2) =      D   (gamma^2)  --------------
!      n                         n               sqrt(1+t)
!
!  for some normalization constant D_n^m(\gamma^2) which is NOT equal to
!  C_n^m(\gamma^2).
!
!  Psi and its first few derivatives are represented via piecewise Chebyshev
!  expansions on an interval of the form (0,b) with b chosen so that the
!  phase function is infintesmal at b.
!  
!  The following subroutines should be regarded as publicly callable:
!
!    spheroidal_phase1 - given the values of the parameters chi, dmu and gamma
!      in (2), find the values of the phase function Psi and its first few
!      derivatives at 0 and construct the phase function on an interval of the
!      form (0,b) with b chosen so that the derivatie of the phase function is
!      of extremely small magnitude at b
!
!    spheroidal_phase2 - given the values of the phase function and its first few
!      derivatives at 0, construct piecewise Chebyshev expansions of the phase
!      function and its first few deriviatives on an interval of the form (0,b)
!      where b is chosen so that either the phase function's derivative has
!      decayed to magnitude less than machine precision by b or b = 7 if that
!      does not happen before the point 7.
!
!    spheroidal_phase_eval - use the phase function Psi to evaluate the functions
!      (3) and (5) at a specified point in [0,1)
!
!    spheroidal_phase_eval_psi - evaluate the phase function Psi and its first few
!      derivatives 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module spheroidal_phase

use utils
use chebyshev

!  The structure spheroidal_phase_data stores the phase function and all necessary
!  related data.
type         spheroidal_phase_data

! The parameters in (1)
double precision              :: chi
double precision              :: gamma
double precision              :: dmu       
double precision              :: dnorm             ! the normalizing factor for ps


integer                       :: k                 ! order of the piecewise Chebyshev expansions
                                                   ! used to represent phase functions
!double precision, allocatable :: xs(:)             ! the nodes of the kth order Chebyshev grid

! Representation of the phase function for (1) on the interval (0,b)
integer                       :: nints
double precision, allocatable :: ab(:,:)

! Piecewise Chebyshev expansions of the phase function and its first few derivatives
double precision, allocatable :: acoefs(:,:)
double precision, allocatable :: apcoefs(:,:)
double precision, allocatable :: appcoefs(:,:)
double precision, allocatable :: apppcoefs(:,:)

! the values of the nonoscillatory phase function and its derivatives at 0
double precision              :: aval
double precision              :: apval
double precision              :: appval
double precision              :: apppval

end type     spheroidal_phase_data


! parameters for the nonlinear solver
integer, private, parameter              :: niters = 10
double precision, parameter, private     :: dtail = .33d0

contains


subroutine spheroidal_phase1(gamma,dmu,chi,phase)
implicit double precision (a-h,o-z)
!type(chebexps_data)                                        :: chebdata
type(spheroidal_phase_data), intent(out)                   :: phase
!
!  Compute the values of the phase function Psi and its first few derivatives
!  at 0 given the values of the parameters gamma, dmu and chi in (2).
!
!  This subroutine operates by first solving the Riccati equation associated with
!  (1) down the imaginary axis in order to compute the values of the first
!  few derivatives of Psi at 0.  To calculate the value of Psi at 0, we 
!  construct a phase function for  the equation
!
!    z''(x) + q(x) z(x)
!
!  obtained from (2) via the change of variables
!
!    z(x) = y( phi(x) ) r(x)                                                             
!
!  with 

!    phi(x) = 1 - exp(-x)     and    r(x)  = sqrt(2 - exp(-x) ).
!
!  This exponential change of variables allows us to get very close to the
!  singularity at 1 and ensure Psi is 0 there.
!
!  Input parameters:
!    chebdata - a structure returned by the chebexps routine which, among other
!      things, specifies the order k of the Chebyshev expansions
!    (gamma,dmu,chi) - the parameters in (1)
!
!  Output parameters:
!    phase - a data structure which will hold the values of the phase
!      function and its first few derivatives upon the completion
!      of this routine
!
type(chebexps_data)                      :: chebdata
double precision, allocatable            :: vals(:), vals2(:), vals3(:)

! for the solution of Appell's equation on (0,1)
double precision, allocatable            :: ws(:,:), wders(:,:), wder2s(:,:)
double precision, allocatable            :: wder3s(:,:), ab(:,:)

! for the solution of Riccati's equation on the imaginary axis
double precision, allocatable            :: ab_r(:,:)
double precision, allocatable            :: rs(:,:),rders(:,:),rints(:,:)
double complex                           :: ima


eps0           = epsilon(0.0d0)
eps            = eps0*10*gamma

ima            = (0.0d0,1.0d0)
pi             = acos(-1.0d0)
k              = 30
call chebexps(k,chebdata)

phase%gamma    = gamma
phase%dmu      = dmu
phase%k        = k

!allocate(phase%xs(k))
!phase%xs  = chebdata%xs

!  Solve Riccati's equation along the imaginary axis
call spheroidal_riccati_tvp(ier,eps,chebdata,gamma,chi,dmu,nints_r,ab_r,&
  rs,rders,rints)


if (ier .ne. 0) then
call prini("after riccati_tvp_adap, ier = ",ier)
call prind("gamma   = ",gamma)
call prin2("eps     = ",eps)
call prind("chi     = ",chi)
call prind("dmu     = ",dmu)
stop
endif

rval = rs(1,1)
rder = rders(1,1)

!  Compute the values of the phase function and its derivatives at 0 under
!  the change of variable isubst = 0

apval0   = -rval
appval0  = 0
apppval0 = -rder * 2 * rval


!  Apply the transformation to the exponential change of variables and then
!  convert the values of the first three derivatives of \alpha(x) to
!  the values of w(0), w'(0) and w''(0)

isubst         = 2
call spheroidal_convert_from0(isubst,apval0,appval0,apppval0,apval,appval,apppval)
call spheroidal_convert_to_w(apval,appval,apppval,wa,wpa,wppa)


!  Solve Appell's equation going forward using an adaptive Chebyshev spectral
!  method

call elapsed(t1)
call spheroidal_appell(ier,eps,chebdata,gamma,chi,dmu,wa,wpa,wppa, &
  nints,phase%ab,ws,wders,wder2s,wder3s)
call elapsed(t2)

if (ier .ne. 0) then
call prini("after spheroidal_appell, ier = ",ier)
call prind("gamma   = ",gamma)
call prin2("eps     = ",eps)
call prind("chi     = ",chi)
call prind("dmu     = ",dmu)
stop
endif

phase%nints  = nints
phase%k      = k


!  Construct the expansions of the derivatives of alpha

allocate(phase%acoefs(k,nints))
allocate(phase%apcoefs(k,nints))
allocate(phase%appcoefs(k,nints))
allocate(phase%apppcoefs(k,nints))

allocate(vals(k),vals2(k),vals3(k))

aval = 0.0d0

do int=nints,1,-1
a0       =  phase%ab(1,int)
b0       =  phase%ab(2,int)

vals     = 1.0d0/ws(:,int)
vals2    = aval + (b0-a0)/2 * matmul(chebdata%aintr,vals)
aval     = vals2(1)

phase%acoefs(:,int)   = matmul(chebdata%u,vals2)
phase%apcoefs(:,int)  = matmul(chebdata%u,vals)

vals     = 1.0d0/ws(:,int)**2
vals2    = 1.0d0/ws(:,int)**3
vals3    = wders(:,int)**2

phase%appcoefs(:,int)  = matmul(chebdata%u,-wders(:,int)/vals)
phase%apppcoefs(:,int) = matmul(chebdata%u,-wder2s(:,int)/vals + 2*vals3/vals2)

end do

phase%aval    = aval
phase%apval   = apval
phase%appval  = appval
phase%apppval = apppval
phase%chi     = chi


end subroutine


subroutine spheroidal_phase2(phase)
implicit double precision (a-h,o-z)
type(spheroidal_phase_data), intent(inout)                 :: phase
!
!  Construct the nonoscillatory phase function Psi given gamma, dmu and the values
!  of Psi and its first few derivatives at 0. 
!
!  Input parameters:
!    phase - the data structure whose entries gamma, dmu, chi specify
!     the parameters is (1) and whose entries aval, apval, appval anmd apppval
!     give the value of Psi and its first few derivatives at 0
!
!  Output parameters:
!    phase - the entries of the data structure which store the Chebyshev expansions
!      of the phase function and its derivatives will be set upon return
!
type(chebexps_data)                      :: chebdata
double precision, allocatable            :: vals(:),vals2(:)

! for the solution of Appell's equation
double precision, allocatable            :: ws(:,:), wders(:,:), wder2s(:,:)
double precision, allocatable            :: wder3s(:,:)


k        = phase%k
call chebexps(k,chebdata)

gamma     = phase%gamma
dmu       = phase%dmu
chi       = phase%chi

aval      = phase%aval
apval     = phase%apval
appval    = phase%appval
apppval   = phase%apppval

eps0           = epsilon(0.0d0)
eps            = eps0*10*gamma

!isubst    = 2
! call spheroidal_convert_from0(isubst,apval0,appval0,apppval0, &
!   apval,appval,apppval)

!  Convert the values of the first three derivatives of \alpha(x) to
!  the values of w(0), w'(0) and w''(0)
call spheroidal_convert_to_w(apval,appval,apppval,wa,wpa,wppa)

call elapsed(t1)
call spheroidal_appell2(ier,eps,chebdata,gamma,chi,dmu,wa,wpa,wppa, &
  nints,phase%ab,ws,wders,wder2s,wder3s)
call elapsed(t2)

if (ier .ne. 0) then
call prini("after spheroidal_appell2, ier = ",ier)
call prind("c       = ",c)
call prin2("eps     = ",eps)
call prind("chi     = ",chi)
call prind("dmu     = ",dmu)
stop
endif

phase%nints = nints
phase%k     = k

!  Construct the expansions of alpha and alpha' going forward
if (allocated(phase%acoefs))  deallocate(phase%acoefs)
if (allocated(phase%apcoefs)) deallocate(phase%apcoefs)

allocate(phase%acoefs(k,nints))
allocate(phase%apcoefs(k,nints))
allocate(vals(k),vals2(k))

do int=1,nints
a0       =  phase%ab(1,int)
b0       =  phase%ab(2,int)
vals     = 1.0d0/ws(:,int)

vals2    = aval + (b0-a0)/2 * matmul(chebdata%aintl,vals)
aval     = vals2(k)

phase%acoefs(:,int)  = matmul(chebdata%u,vals2)
phase%apcoefs(:,int) = matmul(chebdata%u,vals)

end do


end subroutine



subroutine spheroidal_phase_eval_psi(phase,x,aval,apval,appval)
implicit double precision (a-h,o-z)
type(spheroidal_phase_data), intent(in)                   :: phase
!
!  Evaluate the phase function constructed by one of the spheroidal_phase?
!  routines as well as its first few derivatives at a specified point.
!
!  Input parameters:
!    phase - the structure describing the phase function 
!    x - point at which to evaluate the phase function
!
!  Output parameters:
!    aval - the value of the phase function at the point x
!    apval - the derivative of the phase function at x
!

call chebpw_eval22(phase%nints,phase%ab,phase%k,phase%acoefs,phase%apcoefs,&
  x,aval,apval)

end subroutine



subroutine spheroidal_phase_eval(phase,t,valps,valqs)
implicit double precision (a-h,o-z)
type(spheroidal_phase_data), intent(in)                   :: phase
!
!  Evaluate the functions (3) and (5) defined above.  The first of these
!  is a multiple of angular spheroidal wave function of the first kind
!  and the second is a multiple of the angular function of the second
!  kind in the event that chi corresponds to integer values of nu and mu.
!
!  Input parameters:
!    phase - the structure describing the phase function 
!    t - a point in the interval (0,1)
!
!  Output parameters:
!    valps - the value of the function u
!    valqs - the value of the function v
!

! if (phase%isubst .eq. 0) then

! call  chebpw_eval22(phase%nints,phase%ab,phase%k,phase%acoefs,phase%apcoefs,&
!   t,aval,apval)

! dd     = 1.0d0/sqrt(1-t**2) * 1/sqrt(apval)
! valps  = sin(aval)*dd
! valqs  = cos(aval)*dd

! else

x = -log(1-t)
call  chebpw_eval22(phase%nints,phase%ab,phase%k,phase%acoefs,phase%apcoefs,&
  x,aval,apval)
dd   = 1/sqrt(1+t)
valps  = sin(aval)/sqrt(apval)*dd
valqs  = cos(aval)/sqrt(apval)*dd

! endif


end subroutine


subroutine spheroidal_appell(ier,eps,chebdata,gamma,chi,dmu,wa,wpa,wppa, &
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
!  where q(x) one of the coefficients in (1), on an of the form interval (0,b) using an
!  adaptive Chebyshev spectral method.  The function
!
!    w(x) = 1 / alpha'(x)
!
!  satisfies the Appell equation.  The value of b is chosen so that the magnitude of 
!  alpha'(x) is extremely small at b (close to the smallest value representable using
!  double precision numbers).
!
!  Input parameters:
!    eps - precision for the adaptive discretization procedures
!    chebdata - the structure returned by the chebexps routine; the number of
!      Chebyshev nodes k used per discretization interval is specified by this
!      data structure
!    (gamma,chi,dmu) - the parameters in (1)
!
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
b = 1d12

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
call qspheroidal(gamma,chi,dmu,k,ts,ps,qs)

rs = 0
ps = ps * 4
qs = qs * 2
fs = 0

if (nintsout .eq. 0) then
ys0(1)     = wa
yders0(1)  = wpa
yder2s0(1) = wppa
else
ys0(1)     = ysout(k,nintsout)
yders0(1)  = ydersout(k,nintsout)
yder2s0(1) = yder2sout(k,nintsout)
endif

call spheroidal_ivp_int(a0,b0,k,chebdata%xs,chebdata%aintl,chebdata%aintl2, &
  chebdata%aintl3,rs,ps,qs,fs,ys0,yders0,yder2s0,yder3s0)

ifstop  = 0

if (abs(ys0(k)) .gt. 1d250) ifstop=1

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



subroutine spheroidal_appell2(ier,eps,chebdata,gamma,chi,dmu,wa,wpa,wppa, &
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
!  where q(x) one of the coefficients in (1), on an of the form interval (0,b) using an
!  adaptive Chebyshev spectral method.  The function
!
!    w(x) = 1 / alpha'(x)
!
!  satisfies the Appell equation.  The value of b is chosen so that the magnitude of 
!  alpha'(x) is extremely small at b (close to the smallest value representable using
!  double precision numbers).
!
!  Input parameters:
!    eps - precision for the adaptive discretization procedures
!    chebdata - the structure returned by the chebexps routine; the number of
!      Chebyshev nodes k used per discretization interval is specified by this
!      data structure
!    (gamma,chi,dmu) - the parameters in (1)
!
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


ier  = 0
eps0 = epsilon(0.0d0)

!
!  Set algorithm parameters and allocate memory for the procedure.
!

maxints  = 1000
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
a = 0.0d0
b = 7.0d0

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
call qspheroidal(gamma,chi,dmu,k,ts,ps,qs)

rs = 0
ps = ps * 4
qs = qs * 2
fs = 0

if (nintsout .eq. 0) then
ys0(1)     = wa
yders0(1)  = wpa
yder2s0(1) = wppa
else
ys0(1)     = ysout(k,nintsout)
yders0(1)  = ydersout(k,nintsout)
yder2s0(1) = yder2sout(k,nintsout)
endif

call spheroidal_ivp_int(a0,b0,k,chebdata%xs,chebdata%aintl,chebdata%aintl2, &
  chebdata%aintl3,rs,ps,qs,fs,ys0,yders0,yder2s0,yder3s0)

ifstop  = 0

if (abs(ys0(k)) .gt. 1/eps0) ifstop=1

if (ifstop .eq. 1) exit

ifsplit = 0

coefs0 = matmul(chebdata%u,ys0)
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


! subroutine spheroidal_appell3(ier,eps,chebdata,gamma,chi,dmu,wa,wpa,wppa, &
!   nints,ab,ws,wders,wder2s,wder3s)
! implicit double precision (a-h,o-z)
! type(chebexps_data)                        :: chebdata
! integer, intent(out)                       :: nints
! double precision, allocatable, intent(out) :: ab(:,:), ws(:,:),wders(:,:),wder2s(:,:),wder3s(:,:)
! double precision, intent(in)               :: wa,wpa,wppa
! type(c_ptr)                                :: userptr
! !
! !  Solve an initial value problem for Appell's equation
! !
! !    w'''(x) + 4 q(x) w'(x) + 2 q'(x) w(x) = 0,
! !
! !  where q(x) one of the coefficients in (1), on an of the form interval (0,b) using an
! !  adaptive Chebyshev spectral method.  The function
! !
! !    w(x) = 1 / alpha'(x)
! !
! !  satisfies the Appell equation.  The value of b is chosen so that the magnitude of 
! !  alpha'(x) is extremely small at b (close to the smallest value representable using
! !  double precision numbers).
! !
! !  Input parameters:
! !    eps - precision for the adaptive discretization procedures
! !    chebdata - the structure returned by the chebexps routine; the number of
! !      Chebyshev nodes k used per discretization interval is specified by this
! !      data structure
! !    (gamma,chi,dmu) - the parameters in (1)
! !
! !    wa - the value of w(0)
! !    wpa - the value of w'(0)
! !    wppa -  the value of w''(0)
! !
! !  Output parameters:
! !    nints - the number of intervals used to discretize the solution
! !    ab - a (2,nints) array specifying the intervals used to discretize the
! !      solution
! !    ws - the (k,nints) array of values of w at the discretization nodes
! !    wders - the (k,nints) array of the value of w' at the discretization nodes
! !    wder2s - the (k,nints) array of the value of w'' at the discretization nodes
! !    wder3s - the (k,nints) array of the value of w''' at the discretization nodes
! !

! double precision, allocatable   :: ab0(:,:)               ! the list of intervals to process


! ! data for solving over a single interval
! double precision, allocatable   :: ys0(:),yders0(:),yder2s0(:),yder3s0(:),coefs0(:)
! double precision, allocatable   :: amatr(:,:),ts(:),ps(:),rs(:),qs(:),fs(:)

! ! output data
! double precision, allocatable   :: about(:,:)
! double precision, allocatable   :: ysout(:,:),ydersout(:,:),yder2sout(:,:)
! double precision, allocatable   :: yder3sout(:,:)


! ier  = 0
! eps0 = epsilon(0.0d0)

! !
! !  Set algorithm parameters and allocate memory for the procedure.
! !

! maxints  = 1000
! k        = chebdata%k
! nn       = k * dtail

! allocate(ys0(k),yders0(k),yder2s0(k),yder3s0(k),coefs0(k),amatr(k,k))
! allocate(ab0(2,maxints),rs(k),qs(k),ps(k),ts(k),fs(k)) 
! allocate(about(2,maxints))
! allocate(ysout(k,maxints))
! allocate(ydersout(k,maxints))
! allocate(yder2sout(k,maxints))
! allocate(yder3sout(k,maxints))

! !
! a = 0.0d0
! b = 7.0d0

! dd     = log(b-a)
! nints0 = ceiling(dd)*2
! dd     = dd/(nints0+0.0d0)

! do i=1,nints0
! int       = nints0-i+1
! a0        = dd * (i-1)
! b0        = dd * i
! ab0(1,int) = a0
! ab0(2,int) = b0
! end do

! ab0(:,1:nints0)        = a + exp(ab0(:,1:nints0))-1

! nintsout = 0

! do while (nints0 > 0) 

! ifsplit = 0
! a0      = ab0(1,nints0)
! b0      = ab0(2,nints0)
! nints0  = nints0 - 1

! if (b0-a0 .eq. 0) then
! ier = 256
! call prin2("in solve3_ivp_adap interval of length 0, a = ",a)
! return
! endif


! ts = (a0+b0)/2 + (b0-a0)/2 * chebdata%xs
! call qspheroidal(gamma,chi,dmu,k,ts,ps,qs)

! rs = 0
! ps = ps * 4
! qs = qs * 2
! fs = 0

! if (nintsout .eq. 0) then
! ys0(1)     = wa
! yders0(1)  = wpa
! yder2s0(1) = wppa
! else
! ys0(1)     = ysout(k,nintsout)
! yders0(1)  = ydersout(k,nintsout)
! yder2s0(1) = yder2sout(k,nintsout)
! endif

! call spheroidal_ivp_int(a0,b0,k,chebdata%xs,chebdata%aintl,chebdata%aintl2, &
!   chebdata%aintl3,rs,ps,qs,fs,ys0,yders0,yder2s0,yder3s0)

! ifstop  = 0

! if (abs(ys0(k)) .gt. 1/eps0) ifstop=1

! if (ifstop .eq. 1) exit

! ifsplit = 0

! coefs0 = matmul(chebdata%u,ys0)
! dd1 = maxval(abs(coefs0))
! dd2 = maxval(abs(coefs0(k-nn+1:k)))
! if (dd2 .gt. eps*dd1) ifsplit = 1


! if (ifsplit .eq. 1) then
! if (nints0+2 .gt. maxints) then
! ier = 8
! return
! endif

! nints0 = nints0+1
! ab0(1,nints0) = (a0+b0)/2
! ab0(2,nints0) = b0

! nints0 = nints0+1
! ab0(1,nints0) = a0
! ab0(2,nints0) = (a0+b0)/2


! else

! if (nintsout+1 .gt. maxints) then
! ier = 8
! return
! endif

! nintsout          = nintsout+1
! about(1,nintsout) = a0
! about(2,nintsout) = b0

! ysout(:,nintsout)     = ys0
! ydersout(:,nintsout ) = yders0
! yder2sout(:,nintsout) = yder2s0
! yder3sout(:,nintsout) = yder3s0

! endif


! end do
! !
! !  Copy out the data
! !
! nints = nintsout
! allocate(ab(2,nints))
! allocate(ws(k,nints))
! allocate(wders(k,nints))
! allocate(wder2s(k,nints))
! allocate(wder3s(k,nints))

! ab    = about(:,1:nints)
! ws    = ysout(:,1:nints)
! wders = ydersout(:,1:nints)
! wder2s = yder2sout(:,1:nints)
! wder3s = yder3sout(:,1:nints)

! end subroutine

! subroutine spheroidal_appell3(ier,eps,chebdata,gamma,chi,dmu,wa,wpa,wppa, &
!   nints,ab,ws,wders,wder2s,wder3s)
! implicit double precision (a-h,o-z)
! type(chebexps_data)                        :: chebdata
! integer, intent(out)                       :: nints
! double precision, allocatable, intent(out) :: ab(:,:), ws(:,:),wders(:,:),wder2s(:,:),wder3s(:,:)
! !double precision, intent(in)               :: wa,wpa,wppa
! type(c_ptr)                                :: userptr
! !
! !  Solve an initial value problem for Appell's equation
! !
! !    w'''(x) + 4 q(x) w'(x) + 2 q'(x) w(x) = 0,
! !
! !  where q(x) one of the coefficients in (1), on an of the form interval (0,b) using an
! !  adaptive Chebyshev spectral method.  The function
! !
! !    w(x) = 1 / alpha'(x)
! !
! !  satisfies the Appell equation.  The value of b is chosen so that the magnitude of 
! !  alpha'(x) is extremely small at b (close to the smallest value representable using
! !  double precision numbers).
! !
! !  Input parameters:
! !    eps - precision for the adaptive discretization procedures
! !    chebdata - the structure returned by the chebexps routine; the number of
! !      Chebyshev nodes k used per discretization interval is specified by this
! !      data structure
! !    (gamma,chi,dmu) - the parameters in (1)
! !
! !    wa - the value of w(0)
! !    wpa - the value of w'(0)
! !    wppa -  the value of w''(0)
! !
! !  Output parameters:
! !    nints - the number of intervals used to discretize the solution
! !    ab - a (2,nints) array specifying the intervals used to discretize the
! !      solution
! !    ws - the (k,nints) array of values of w at the discretization nodes
! !    wders - the (k,nints) array of the value of w' at the discretization nodes
! !    wder2s - the (k,nints) array of the value of w'' at the discretization nodes
! !    wder3s - the (k,nints) array of the value of w''' at the discretization nodes
! !

! double precision, allocatable   :: ab0(:,:)               ! the list of intervals to process


! ! data for solving over a single interval
! double precision, allocatable   :: ys0(:),yders0(:),yder2s0(:),yder3s0(:),coefs0(:)
! double precision, allocatable   :: amatr(:,:),ts(:),ps(:),rs(:),qs(:),fs(:)

! ! output data
! double precision, allocatable   :: about(:,:)
! double precision, allocatable   :: ysout(:,:),ydersout(:,:),yder2sout(:,:)
! double precision, allocatable   :: yder3sout(:,:)


! ier = 0

! !
! !  Set algorithm parameters and allocate memory for the procedure.
! !

! eps0     = epsilon(0.0d0)
! maxints  = 1000
! k        = chebdata%k
! nn       = k * dtail

! allocate(ys0(k),yders0(k),yder2s0(k),yder3s0(k),coefs0(k),amatr(k,k))
! allocate(ab0(2,maxints),rs(k),qs(k),ps(k),ts(k),fs(k)) 
! allocate(about(2,maxints))
! allocate(ysout(k,maxints))
! allocate(ydersout(k,maxints))
! allocate(yder2sout(k,maxints))
! allocate(yder3sout(k,maxints))

! !
! a = 0.00d0
! b = 0.99d0

! dd     = b-a
! nints0 = 4
! dd     = dd/(nints0+0.0d0)

! do i=1,nints0
! int       = nints0-i+1
! a0        = dd * (i-1)
! b0        = dd * i
! ab0(1,int) = a0
! ab0(2,int) = b0
! end do


! nintsout = 0

! do while (nints0 > 0) 

! ifsplit = 0
! a0      = ab0(1,nints0)
! b0      = ab0(2,nints0)
! nints0  = nints0 - 1


! if (b0-a0 .eq. 0) then
! ier = 256
! call prin2("in solve3_ivp_adap interval of length 0, a = ",a)
! return
! endif


! ts = (a0+b0)/2 + (b0-a0)/2 * chebdata%xs
! call qspheroidal2(gamma,chi,dmu,k,ts,ps,qs)


! rs = 0
! ps = ps * 4
! qs = qs * 2
! fs = 0


! if (nintsout .eq. 0) then
! ys0(1)     = wa
! yders0(1)  = wpa
! yder2s0(1) = wppa
! else
! ys0(1)     = ysout(k,nintsout)
! yders0(1)  = ydersout(k,nintsout)
! yder2s0(1) = yder2sout(k,nintsout)
! endif

! call spheroidal_ivp_int(a0,b0,k,chebdata%xs,chebdata%aintl,chebdata%aintl2, &
!   chebdata%aintl3,rs,ps,qs,fs,ys0,yders0,yder2s0,yder3s0)

! if (abs(ys0(k)) .gt. 1d15) exit

! ifsplit = 0

! coefs0 = matmul(chebdata%u,ys0)
! dd1 = maxval(abs(coefs0))
! dd2 = maxval(abs(coefs0(k-nn+1:k)))
! if (dd2 .gt. eps*dd1) ifsplit = 1

! ! coefs0 = matmul(chebdata%u,yders0)
! ! dd1 = maxval(abs(coefs0))
! ! dd2 = maxval(abs(coefs0(k-nn+1:k)))
! ! if (dd2 .gt. eps*dd1) ifsplit = 1


! if (ifsplit .eq. 1) then
! if (nints0+2 .gt. maxints) then
! ier = 8
! return
! endif

! nints0 = nints0+1
! ab0(1,nints0) = (a0+b0)/2
! ab0(2,nints0) = b0

! nints0 = nints0+1
! ab0(1,nints0) = a0
! ab0(2,nints0) = (a0+b0)/2

! else

! if (nintsout+1 .gt. maxints) then
! ier = 8
! return
! endif

! nintsout          = nintsout+1
! about(1,nintsout) = a0
! about(2,nintsout) = b0

! ysout(:,nintsout)     = ys0
! ydersout(:,nintsout ) = yders0
! yder2sout(:,nintsout) = yder2s0
! yder3sout(:,nintsout) = yder3s0

! endif


! end do

! !
! !  Copy out the data
! !
! nints = nintsout
! allocate(ab(2,nints))
! allocate(ws(k,nints))
! allocate(wders(k,nints))
! allocate(wder2s(k,nints))
! allocate(wder3s(k,nints))

! ab    = about(:,1:nints)
! ws    = ysout(:,1:nints)
! wders = ydersout(:,1:nints)
! wder2s = yder2sout(:,1:nints)
! wder3s = yder3sout(:,1:nints)

! end subroutine



subroutine qspheroidal(c,chi,dmu,k,xs,vals,ders)
implicit double precision (a-h,o-z)
double precision :: xs(k),vals(k),ders(k)
!
!  Evaluate the coefficient q(x) in the equation (1) and its derivative at
!  a collection of points.
!
!  Input parameters:
!    (c,chi,dmu) - the parameters in (1)
!    k - the number of points at which to evaluate q(x)
!    xs - the point at which to evaluate q(x)
!
!  Output parameters:
!    vals - the values of q(x) 
!    ders - the values of q'(x)
!


do i=1,k
x       = xs(i)
z       = exp(x)


if (z .gt. 1d7) then

val = -dmu**2/0.4d1+((-11*c**2)/0.256d3-(3*(-1+0.2d1*c**2-  &
0.2d1*chi))/0.512d3-(11*(1-0.2d2*c**2+0.4d1*chi))/0.4096d4-  &
(13*dmu**2)/0.4096d4)/(0.4d1*z**12)+((-5*c**2)/0.64d2-  &
(11*(-1+0.2d1*c**2-0.2d1*chi))/0.1024d4-(5*(1-0.2d2*c**2+  &
0.4d1*chi))/0.1024d4-(3*dmu**2)/0.512d3)/(0.4d1*z**11)+  &
((-9*c**2)/0.64d2-(5*(-1+0.2d1*c**2-0.2d1*chi))/0.256d3-  &
(9*(1-0.2d2*c**2+0.4d1*chi))/0.1024d4-  &
(11*dmu**2)/0.1024d4)/(0.4d1*z**10)+(-c**2/0.4d1+(-1+  &
0.2d2*c**2-0.4d1*chi)/0.64d2-(9*(-1+0.2d1*c**2-  &
0.2d1*chi))/0.256d3-(5*dmu**2)/0.256d3)/(0.4d1*z**9)+  &
((-7*c**2)/0.16d2+(1-0.2d1*c**2+0.2d1*chi)/0.16d2-(7*(1-  &
0.2d2*c**2+0.4d1*chi))/0.256d3-  &
(9*dmu**2)/0.256d3)/(0.4d1*z**8)+((-3*c**2)/0.4d1-(7*(-1+  &
0.2d1*c**2-0.2d1*chi))/0.64d2-(3*(1-0.2d2*c**2+  &
0.4d1*chi))/0.64d2-dmu**2/0.16d2)/(0.4d1*z**7)+  &
((-5*c**2)/0.4d1-(3*(-1+0.2d1*c**2-0.2d1*chi))/0.16d2-(5*(1-  &
0.2d2*c**2+0.4d1*chi))/0.64d2-  &
(7*dmu**2)/0.64d2)/(0.4d1*z**6)+(-0.2d1*c**2+(-1+0.2d2*c**2-  &
0.4d1*chi)/0.8d1-(5*(-1+0.2d1*c**2-0.2d1*chi))/0.16d2-  &
(3*dmu**2)/0.16d2)/(0.4d1*z**5)+(-0.3d1*c**2+(1-0.2d1*c**2+  &
0.2d1*chi)/0.2d1-(3*(1-0.2d2*c**2+0.4d1*chi))/0.16d2-  &
(5*dmu**2)/0.16d2)/(0.4d1*z**4)+(-0.4d1*c**2+(-1+0.2d2*c**2-  &
0.4d1*chi)/0.4d1-(3*(-1+0.2d1*c**2-0.2d1*chi))/0.4d1-  &
dmu**2/0.2d1)/(0.4d1*z**3)+(3+0.12d2*c**2+0.4d1*chi-  &
0.3d1*dmu**2)/(0.16d2*z**2)+(1-0.2d1*c**2+0.2d1*chi-  &
0.1d1*dmu**2)/(0.4d1*z)

der = (0.2d1*((33*c**2)/0.512d3+(33*(-9*c**2+chi))/0.4096d4+  &
(39*(-1+2*c**2-2*chi+dmu**2))/0.8192d4))/z**12+  &
(0.2d1*((55*c**2)/0.512d3+(55*(-9*c**2+chi))/0.4096d4+  &
(33*(-1+2*c**2-2*chi+dmu**2))/0.4096d4))/z**11+  &
(0.2d1*((45*c**2)/0.256d3+(45*(-9*c**2+chi))/0.2048d4+  &
(55*(-1+2*c**2-2*chi+dmu**2))/0.4096d4))/z**10+  &
(0.2d1*((9*c**2)/0.32d2+(9*(-9*c**2+chi))/0.256d3+(45*(-1+  &
2*c**2-2*chi+dmu**2))/0.2048d4))/z**9+  &
(0.2d1*((7*c**2)/0.16d2+(7*(-9*c**2+chi))/0.128d3+(9*(-1+  &
2*c**2-2*chi+dmu**2))/0.256d3))/z**8+  &
(0.2d1*((21*c**2)/0.32d2+(21*(-9*c**2+chi))/0.256d3+(7*(-1+  &
2*c**2-2*chi+dmu**2))/0.128d3))/z**7+  &
(0.2d1*((15*c**2)/0.16d2+(15*(-9*c**2+chi))/0.128d3+(21*(-1+  &
2*c**2-2*chi+dmu**2))/0.256d3))/z**6+(0.2d1*((5*c**2)/0.4d1+  &
(5*(-9*c**2+chi))/0.32d2+(15*(-1+2*c**2-2*chi+  &
dmu**2))/0.128d3))/z**5+(0.2d1*((3*c**2)/0.2d1+(3*(-9*c**2+  &
chi))/0.16d2+(5*(-1+2*c**2-2*chi+dmu**2))/0.32d2))/z**4+  &
(3*(-1+c**2-0.1d1*chi+dmu**2))/(0.8d1*z**3)+(-3-0.12d2*c**2-  &
0.4d1*chi+0.3d1*dmu**2)/(0.8d1*z**2)+(-1+0.2d1*c**2-  &
0.2d1*chi+dmu**2)/(0.4d1*z)


! val    = -dmu**2/0.4d1+(5-0.4d1*c**2+0.4d1*chi-  &
!           0.5d1*dmu**2)/(0.64d2*z**4)+(1-0.1d1*c**2+chi-  &
!           0.1d1*dmu**2)/(0.8d1*z**3)+(3+0.12d2*c**2+0.4d1*chi-  &
!           0.3d1*dmu**2)/(0.16d2*z**2)+(1-0.2d1*c**2+0.2d1*chi-  &
!           0.1d1*dmu**2)/(0.4d1*z)

! der    =  (-5+0.4d1*c**2-0.4d1*chi+0.5d1*dmu**2)/(0.16d2*z**5)-(3*(1-  &
!           0.1d1*c**2+chi-0.1d1*dmu**2))/(0.8d1*z**4)+(-3-0.12d2*c**2-  &
!           0.4d1*chi+0.3d1*dmu**2)/(0.8d1*z**3)+(-1+0.2d1*c**2-  &
!           0.2d1*chi+dmu**2)/(0.4d1*z**2)

else

val    = -(1+0.4d1*chi+0.4d1*c**2*exp(-2*x)*(-1+exp(x))**2*(-1+  &
          2*exp(x))+0.4d1*exp(x)*(-1-2*chi+dmu**2*exp(x)))/(0.4d1*(1-  &
          0.2d1*exp(x))**2)

der    = (-0.1d1*exp(x)*(1+4*chi+4*c**2*exp(-2*x)*(-1+exp(x))**2*(-1+  &
          2*exp(x))+4*exp(x)*(-1-2*chi+dmu**2*exp(x))))/(1-  &
          2*exp(x))**3-(0.8d1*c**2*exp(-x)*(-1+exp(x))**2+  &
          0.8d1*c**2*exp(-x)*(-1+exp(x))*(-1+2*exp(x))-  &
          0.8d1*c**2*exp(-2*x)*(-1+exp(x))**2*(-1+2*exp(x))+  &
          0.4d1*exp(x)*(-1-2*chi+dmu**2*exp(x))+  &
          0.4d1*dmu**2*exp(2*x))/(0.4d1*(1-0.2d1*exp(x))**2)

endif

vals(i) = val
ders(i) = der

end do



end subroutine


subroutine qspheroidal2(gamma,chi,dmu,k,xs,vals,ders)
implicit double precision (a-h,o-z)
double precision :: xs(k),vals(k),ders(k)
!
!  Evaluate the coefficient q(x) in the spheroidal wave equation.
!
!  Input parameters:
!    (gamma,chi,dmu) - the parameters in (1)
!    k - the number of points at which to evaluate q(x)
!    xs - the point at which to evaluate q(x)
!
!  Output parameters:
!    vals - the values of q(x) 
!    ders - the values of q'(x)
!

do i=1,k
x       = xs(i)
val = (1-0.1d1*dmu**2)/(1-0.1d1*x**2)**2+(chi-  &
      0.1d1*gamma**2*x**2)/(1-0.1d1*x**2)
der = (0.4d1*(1-dmu**2)*x)/(1-x**2)**3-(0.2d1*gamma**2*x)/(1-x**2)  &
      +(0.2d1*x*(chi-gamma**2*x**2))/(1-x**2)**2
vals(i) = val
ders(i) = der
end do

end subroutine


subroutine spheroidal_convert_from0(isubst,apval0,appval0,apppval0, &
  apval,appval,apppval)
implicit double precision (a-h,o-z)
!
!  Compute the values of the first three derivatives of a phase function \alpha at 0
!  for a differential equation of the form
!
!    z''(x) + q(x) z(x) = 0
!
!  arising from (1) via a change of variables of the form
!
!    z(x) = y( phi(x) ) r(x)
!
!  given the values of the first three derivatives of the phase function \alpha0
!  for the equation obtained by (1) via the change of variables
!
!    phi(x) = x,  r(x) = sqrt(1-x^2).
!
!  Input parameters:
!    isubst - integer parameter specifying the change of variables;
!      isubst = 0    means   \psi(x) = x,              r(x) = \sqrt(1-x^2)
!      isubst = 1    menas   \psi(x) = w/sqrt(1+w^2),  r(x) = (1+w^2)^(1/4)
!      isubst = 2    means   \psi(x) = 1 - exp(-w),    r(x) =  sqrt(2 - \exp(-x))
!  
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


subroutine spheroidal_convert_to0(isubst,apval,appval,apppval, &
  apval0,appval0,apppval0)
implicit double precision (a-h,o-z)
!
!  Perform the inverse of the operation performed by convert_from0.
!
!  Input parameters:
!    isubst - integer parameter specifying the change of variables (see above)
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



subroutine spheroidal_convert_to_w(apval,appval,apppval,wval,wpval,wppval)
implicit double precision (a-h,o-z)
!
!  Given the values alpha'(x), alpha''(x) and alpha''(x), compute the value
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
wppval = 2*appval**2/apval**3 - apppval / apval**2


end subroutine



subroutine spheroidal_convert_from_w(wval,wpval,wppval,apval,appval,apppval)
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

! subroutine spheroidal_phase_norm(phase,dnorm)
! implicit double precision (a-h,o-z)
! type(spheroidal_phase_data)       :: phase
! !
! !  Compute the integral
! !
! !        1   (    mu     )^2
! !    \int    (  ps  (x)  )    dx
! !        0   (    nu     )
! !
! !  given the phase function for (1).  The integral is computed via an adaptive procedure
! !  which DOES NOT RUN IN TIME INDEPENDENT OF THE PARAMETERS.
! !
! !  Input parameters:
! !    phase - the data structure describing the phase function
! !
! !  Output parameters:
! !    dnorm - the value of the integral
! !
! !
! type(chebexps_data)           :: chebdata
! double precision, allocatable :: stack(:,:), vals(:)

! eps0           = epsilon(0.0d0)
! eps            = eps0*100

! dnorm          = 0

! k = 30
! call chebexps(k,chebdata)

! nints = phase%nints

! maxstack = nints+100000
! allocate(stack(3,maxstack), vals(k) )

! nstack       = 0
! do int=1,nints
! a            = phase%ab(1,int)
! b            = phase%ab(2,int)
! print *,a,b

! val0 = 0
! do i=1,k
! x   = (b-a)/2 * chebdata%xs(i) + (b+a)/2
! wht = (b-a)/2 * chebdata%whts(i)
! call spheroidal_phase_eval_psi(phase,x,aval,apval,appval)
! valps = sin(aval)/sqrt(apval)
! val0  = val0 + valps**2 / (2*exp(x)-1) * wht
! end do

! if (val0 .ne. 0) then
! nstack          = nstack+1
! stack(1,nstack) = a
! stack(2,nstack) = b
! stack(3,nstack) = val0
! endif

! end do

! do while (nstack > 0) 

! a      = stack(1,nstack)
! b      = stack(2,nstack)
! val0   = stack(3,nstack)

! nstack = nstack-1
! c      = (a+b)/2

! ! val0 = 0
! ! do i=1,k
! ! x   = (b-a)/2 * chebdata%xs(i) + (b+a)/2
! ! wht = (b-a)/2 * chebdata%whts(i)
! ! call spheroidal_phase_eval_psi(phase,x,aval,apval,appval)
! ! valps = sin(aval)/sqrt(apval)
! ! val0  = val0 + valps**2  *wht
! ! end do

! val1 = 0
! do i=1,k
! x   = (c-a)/2 * chebdata%xs(i) + (c+a)/2
! wht = (c-a)/2 * chebdata%whts(i)
! call spheroidal_phase_eval_psi(phase,x,aval,apval,appval)
! valps = sin(aval)/sqrt(apval)
! val1  = val1 + valps**2 /(2*exp(x)-1) *wht
! end do

! val2 = 0
! do i=1,k
! x   = (b-c)/2 * chebdata%xs(i) + (b+c)/2
! wht = (b-c)/2 * chebdata%whts(i)
! call spheroidal_phase_eval_psi(phase,x,aval,apval,appval)
! valps = sin(aval)/sqrt(apval)
! val2  = val2 + valps**2 / (2*exp(x)-1) * wht
! end do

! diff = abs(val0 - val1 - val2)
! if (diff .lt. eps) then
! dnorm = dnorm+val0
! else

! if (nstack+2 .gt. maxstack) then
! call prina("dnorm computation stack overflow")
! stop
! endif

! nstack=nstack+1
! stack(1,nstack) = a
! stack(2,nstack) = c
! stack(3,nstack) = val1

! nstack=nstack+1
! stack(1,nstack) = c
! stack(2,nstack) = b
! stack(3,nstack) = val2

! endif

! end do

! dnorm = dnorm*2

! end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  The following code solves the Riccati equation corresponding to the spheroidal
!  wave equation going down the imaginary axis
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine spheroidal_riccati_tvp(ier,eps,chebdata,gamma,chi,dmu,nints,ab,rs,rders,rints)
implicit double precision (a-h,o-z)
type(chebexps_data)                           :: chebdata
integer                                       :: nints 
double precision, allocatable, intent(out)    :: ab(:,:)
double precision, allocatable, intent(out)    :: rs(:,:),rders(:,:),rints(:,:)
!
!  Construct the solution of the Riccati equation
!
!    r'(t) + (r(t))^2 + q(t) = 0,                                                         (10)
!
!  where
!
!               1-dmu^2       chi + gamma^2 t^2 
!    q(t) =  - ---------   -  ----------------- ,                                         (11)
!               (1+t)^2              (1+t^2)
!
!  such that
!
!    r(t) ~ -gamma    as t -> +- infty.
!
!  The solution of (10) and its derivative are represented as piecewise
!  Chebyshev expansions and constructed on a interval of the form
!  (-b,b) with b large.  The set of discretization intervals is determined
!  adaptively.  
!
!  The function r(t) is the logarithmic derivative of the radial spheroidal wave 
!  function of the third kind.
!
!  Input parameters:
!     eps - precision for the adaptive discretization procedure
!     (gamma,chi,dmu) - the parameters in (1)
!     chebdata - the structure returned by chebexps which, among other things,
!       specifies the order of the piecewise Chebyshev expansions to use
!
!  Output parameters:
!    ier - an error return code
!      ier = 0       indicates successful execution
!      ier = 4       means that the adaptive discretization of q failed
!      ier = 8       means that the maximum number of intervals was exceeded
!      ier = 16      means that an interval of length 0 was encountered
!
!    nints - the number of subintervals used to discretize the solution
!    ab - a (2,nints) array specifying the endpoints of said subintervals
!    rs - a (k,nints) array giving the values of the solution at the Chebyshev
!       nodes on each subinterval
!    rders - a (k,nints) array giving the values of the derivatives of the
!       solution at the Chebyshev nodes on each subinterval
!    rints - a (k,nints) array giving the antiderivative of the solutions
!       r(t) such that r(0) = 0
!
! 
double precision, allocatable         :: ab0(:,:),ts0(:)
double precision, allocatable         :: rs0(:),rders0(:),coefs0(:),qs0(:),amatr0(:,:)
double precision, allocatable         :: ps0(:),deltap(:),delta(:),fs0(:)
double precision                      :: rb0

double precision, allocatable         :: about(:,:)
double precision, allocatable         :: rsout(:,:),rdersout(:,:)
double precision                      :: rb

ier       = 0
k         = chebdata%k
ntail     = k * dtail
maxints   = 1000

eps0      = epsilon(0.0d0)

a         =  0.0d0
b         =  1.0d12
rb        = -gamma

allocate(rs0(k),rders0(k),coefs0(k),ts0(k),qs0(k),delta(k),deltap(k))
allocate(ps0(k),fs0(k),amatr0(k,k))

nints0   = 0
allocate(ab0(2,maxints))

nintsout = 0
allocate(about(2,maxints), rsout(k,maxints), rdersout(k,maxints))


!  Adaptively discretize the function q to form an initial set of intervals
!  WARNING: the initial list of intervals must appear in reverse order

nints0   = 1
ab0(1,1) = a
ab0(2,1) = b
nintsout = 0

do while (nints0 > 0)

a0 = ab0(1,nints0)
b0 = ab0(2,nints0)
nints0 = nints0-1

if ( (b0 - a0) .eq. 0) then
ier = 16
call prin2("in spheroidal_riccati_tvp,  zero length interval enountered while discretizing q, a0 = ",a0)
return
endif

ts0 = (b0-a0)/2 * chebdata%xs + (b0+a0)/2
call spheroidal_qaxis(gamma,chi,dmu,k,ts0,qs0)

coefs0 = matmul(chebdata%u,qs0)
dd1    = maxval(abs(coefs0))
dd2    = maxval(abs(coefs0(k-ntail+1:k)))

if (dd2 .lt. dd1 * eps) then

if(nintsout+1 .ge. maxints) then
ier = 4
return
endif

nintsout = nintsout+1
about(1,nintsout) = a0
about(2,nintsout) = b0
else

if (nints0 + 2 .ge. maxints) then
ier = 4
return
endif

nints0 = nints0+1
ab0(1,nints0) = (a0+b0)/2
ab0(2,nints0) = b0

nints0 = nints0+1
ab0(1,nints0) = a0
ab0(2,nints0) = (a0+b0)/2

endif

end do



!  Copy out the intervals
nints0 = nintsout

do int=1,nints0
int0 = int
ab0(:,int)    = about(:,int0)
end do
nintsout = 0


!
!  Now adaptively solve the initial value problem for (1)
!

do while (nints0 > 0 )
a0 = ab0(1,nints0)
b0 = ab0(2,nints0)
nints0 = nints0-1


if ( (b0 - a0) .eq. 0) then
ier = 16
call prin2("in spheroidal_riccati_tvp,  zero length interval enountered, a0 = ",a0)
return
endif

ts0 = (b0-a0)/2 * chebdata%xs + (b0+a0)/2
call  spheroidal_qaxis(gamma,chi,dmu,k,ts0,qs0)


if (nintsout .eq. 0) then
rb0 = rb
else
rb0 = rsout(1,nintsout)
endif


! use the trapezoidal rule to construct an initial guess for Newton iterations
call spheroidal_riccati_trap_tvp(ier,k,ts0,qs0,rb0,rs0,rders0)

if (ier .eq. 1024) then
call prin2("in spheroidal_riccati_tvp, NaN encounted, a0 = ",a0)
call prin2("in spheroidal_riccati_tvp, NaN encounted, b0 = ",b0)
return
endif

! integrate the derivative to obtain the values of rs --- this step is crucial
rs0 = rb0 + matmul(chebdata%aintr*(b0-a0)/2,rders0)

! perform Newton iterations 
do iter=1,niters

ps0 = 2 *rs0
fs0 = -rders0 - rs0**2 - qs0
rb0  = 0
call spheroidal_riccati_linear_tvp(k,a0,b0,ts0,chebdata%aintr,ps0,fs0,rb0,amatr0,delta,deltap)

rs0    = rs0 + delta
rders0 = rders0 + deltap

end do


ifsplit = 0

coefs0 = matmul(chebdata%u,rs0)
dd1    = maxval(abs(coefs0))
dd2    = maxval(abs(coefs0(k-ntail+1:k)))
if (dd2 .lt. eps * dd1) ifsplit = 1


coefs0 = matmul(chebdata%u,rders0)
dd1    = maxval(abs(coefs0))
dd2    = maxval(abs(coefs0(k-ntail+1:k)))
if (dd2 .lt. eps * dd1) ifsplit = 1

if(ifsplit .eq. 1) then

if (nintsout+1 .gt. maxints) then
ier = 8
return
endif

nintsout             = nintsout+1
about(1,nintsout)    = a0
about(2,nintsout)    = b0
rsout(:,nintsout)    = rs0
rdersout(:,nintsout) = rders0
else


if (nints0+2 .gt. maxints) then
ier = 8
return
endif

nints0               = nints0+1
ab0(1,nints0)        = a0
ab0(2,nints0)        = (a0+b0)/2

nints0               = nints0+1
ab0(1,nints0)        = (a0+b0)/2
ab0(2,nints0)        = b0


endif

end do

!
!  Copy out the intervals in reverse order
!


nints = nintsout
allocate(ab(2,nintsout))
allocate(rs(k,nintsout))
allocate(rders(k,nintsout))
allocate(rints(k,nintsout))

do int=1,nints
int0             = nints-int+1
ab(:,int)        = about(:,int0)
rs(:,int)        = rsout(:,int0)
rders(:,int)     = rdersout(:,int0)
end do

!
!  Compute the desired antiderivative
!

do int=1,nints
if ( ab(1,int) == 0.0d0) then
int0 = int
exit
endif
end do

rval = 0
do int=int0,nints
a0           = ab(1,int)
b0           = ab(2,int)
rints(:,int) = rval + (b0-a0)/2 * matmul(chebdata%aintl,rs(:,int))
rval         = rints(k,int)
end do

rval = 0
do int=int0-1,1,-1
a0           = ab(1,int)
b0           = ab(2,int)
rints(:,int) = rval + (b0-a0)/2 * matmul(chebdata%aintr,rs(:,int))
rval         = rints(1,int)
end do


end subroutine


subroutine spheroidal_riccati_trap_tvp(ier,k,ts,qs,rb,rs,rders)
implicit double precision (a-h,o-z)
integer                  :: k
double precision         :: ts(k)
double precision           :: rs(k),rders(k), qs(k), rb
!
!  Use the implicit trapezoidal method to solve the terminal value problem
! 
!     r'(t) + (r(t))^2 + q(t) = 0      a < t < b
!     r(b)                    = rb
!
!  The solution is tabulated at a collection of nodes in the interval (a,b)
!  specified by the user the right-most of which must be the right  endpoint
!  b.
!
!  Input parameters:
!    k - the number of nodes at which to approximate the solution
!    ts - an array containing the nodes
!    qs - the values of the function q at the nodes
!    rb - the terminal value for the solution
!
!  Output parameters:
!    ier - an error return code;
!      ier =    0    indicates successful execution
!      ier = 1024    means that NaN was encountered and the procedure aborted
!
!    rs - the values of the solution at the specified nodes
!    rders - the values of the derivative of the solution at the specified
!      nodes
!    
double precision   :: f0,r,rhs,delta,q0,q1

ier    = 0
!niters = 4       ! the number of Newton iterations for each step

rs(k)    = rb
rders(k) = -qs(k) - rb**2

do j=k-1,1,-1

t0 = ts(j)
t1 = ts(j+1)
h  = t1-t0

q0 = qs(j)
q1 = qs(j+1)

r   = rs(j+1)
rhs = rs(j+1) + h/2 * rs(j+1)**2 + h/2 * q1 + h/2*q0

do iter=1,niters
delta = (rhs - r +h/2*r**2 ) / (1 - h*r)
r     = r + delta
end do

rs(j)    = r
rders(j) = -r**2 - q0
end do


end subroutine


subroutine spheroidal_riccati_linear_tvp(k,a,b,ts,chebint,qs,fs,rb,amatr,rs,rders)
implicit double precision (a-h,o-z)
integer                  :: k
double precision         :: ts(k),chebint(k,k)
double precision           :: amatr(k,k),qs(k),fs(k),rs(k),rders(k),rb
!
!  Use a Chebyshev spectral method to solve the terminal value problem
!
!    r'(t) + qs(t) r(t) = f(t)
!    r(b)               = rb
!
!  over the interval (a,b).
!
!  Input parameters:
!    k - the number of Chebyshev nodes
!    ts - the k Chebyshev nodes on the interval (a,b)
!    chebint - the "right" Chebyshev spectral integration matrix
!    qs - the values of the coefficient at the Chebyshev nodes
!    ra - the initial value for the solution
!
!  Output parameters:
!    rs - the values of the solution at the  Chebyshevnodes
!    rders - the values of the derivative of the solution at the 
!      Chebyshev nodes
!
!  Work arrays:
!    amatr - a (k,k) matrix 
!

!
!  We solve the integral equation 
!
!     \sigma(t) + q(t) \int_b^t \sigma(s) ds = f(t) - rb q(t)
!
!  obtained by letting
!
!     r(t) = rb + \int_b^t \sigma(s) ds.
! 


rders = 0
amatr = 0

do i=1,k
amatr(i,i) = 1.0d0
end do

do i=1,k
rders(i)   = fs(i) - rb * qs(i)
amatr(i,:) = amatr(i,:) + qs(i) * chebint(i,:)*(b-a)/2
end do


call qrsolv(amatr,k,rders)
rs    = rb + (b-a)/2*matmul(chebint,rders)

end subroutine



subroutine spheroidal_ivp_int(a,b,k,xscheb,chebint,chebint2,chebint3,rs,ps,qs,&
  fs,ys,yders,yder2s,yder3s)
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
!    chebint? - the "left" Chebyshev spectral integration matrices as constructed
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


subroutine spheroidal_qaxis(gamma,chi,dmu,k,ts,qs)
implicit double precision (a-h,o-z)
integer                              :: k
double precision                     :: ts(k)
double precision                     :: qs(k)
qs = - ( 1.0d0 - dmu**2 ) / (1+ts**2)**2 - (chi + gamma**2 * ts**2) / ( 1 + ts**2)
end subroutine




end module
