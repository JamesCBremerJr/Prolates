!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for constructing a nonoscillatory phase function
!  for the differential equation 
!
!     z''(x) + q(x) z(x) = 0                                                              (1)
!
!  obtained from the reduced spheroidal wave equation 
!
!    (1-t^2) * y''(t) - 2x y'(t) + ( chi - gamma^2 t^2 ) y(t) = 0                         (2)
!
!  via the transformation
!
!    z(x) = y(1-exp(-x)) sqrt(2-exp(-x)).
!
!  The solutions of (1) are typically indexed by gamma and by an implicit parameter
!  dnu, called the characteristic exponent.  The relationship between the parameter
!  chi in (1) and dnu is fairly complicated (for instance, chi has branch  points at 
!  each half integer value of dnu) and the code here provide no means of computing the 
!  desired value of chi as a function of them.
!
!  See, however, prolates_exps.f90 which provides a mechanism for computing the values
!  of chi corresponding to the prolate spheroidal wave functions of order 0 and
!  integer characteristic exponents.
!
!  The nonoscillatory phase function Psi produced by this routine is determined
!  up to an additive constant by the fact that it corresponds to a solution of (2) 
!  which is in the Hardy space H^2 of functions analytic on the upper half of 
!  the complex plane whose boundary values are square integrable.  We fix the 
!  constant through the requirement that
!
!     \lim   Psi(x) = 0.
!    x \to 1-
!
!  This ensures that 
!
!                    sin (  Psi(x) )  
!    u(x) =       -------------------                                                     (3)
!                   sqrt(  Psi'(x) )
!
!  is related to angular prolate spheroidal wave function of the first kind via
!  the formula
!
!                                               u(-log(1-t))
!    ps    (t; gamma^2) =      C   (gamma^2)  ----------------                            (4)
!      dnu                       dnu              sqrt(1+t)
!
!  with C_dnu(\gamma^2) a normalization constant.  When dnu is an integer
!  it is the case that
!
!                                              v(-log(1-t))
!    qs    (t; gamma^2) =      D   (gamma^2)  -------------- ,
!      n                         n               sqrt(1+t)
!
!  where
!
!                    cos (  Psi(x) )  
!    v(x) =       -------------------                                                     (3)
!                   sqrt(  Psi'(x) )
!
!  and D_n(\gamma^2) is a normalization constant NOT equal to C_n(\gamma^2).
!
!  Psi and its first few derivatives are represented via piecewise Chebyshev
!  expansions on an interval of the form (0,b) with b chosen so that the
!  phase function is infintesmal at b.
!  
!  The following subroutines should be regarded as publicly callable:
!
!    prolates_phase1 - given the values of the parameters chi and gamma
!      in (2), find the values of the phase function Psi and its first few
!      derivatives at 0 and construct the phase function on an interval of the
!      form (0,b) with b chosen so that the derivatie of the phase function is
!      of extremely small magnitude at b
!
!    prolates_phase2 - given the values of the phase function and its first few
!      derivatives at 0, caclulate the phase function and its first few deriviatives 
!      on an interval of the form (0,b) where b is chosen so that either the 
!      phase function's derivative has decayed to magnitude less than machine 
!      precision by b or b = 7 if that does not happen before the point 7.
!
!    prolates_phase3 - given the values of the parameters chi and gamma
!      in (2), return only the values of the phase function and its first
!      few derivatives at 0
!
!    prolates_phase_eval_psi - evaluate the phase function Psi and its first 
!      derivative at a specified point x in the interval (0,b)
!
!    prolates_phase_eval_ps - evaluate the function (4) -- SANS NORMALIZATION CONSTANT -- 
!      at a point t in the interval (-1,1)
!
!    prolates_phase_eval_psder - evaluate the function (4) -- SANS NORMALIZATION CONSTANT -- 
!      and its derivative at a point t in the interval (-1,1)
!
!    prolates_phase_eval0 - evaluate the function (4) -- SANS NORMALIZATION CONSTANT --
!      and its derivative at the point 0; this routine is typically used to compute 
!      the normalization constant in (4)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module prolates_phase

use utils
use chebyshev
use prolates

!  The structure prolates_phase_data stores the phase function and all necessary
!  related data.
type         prolates_phase_data

! The parameters in (1)
double precision              :: chi
double precision              :: gamma
integer                       :: k                 ! order of the piecewise Chebyshev expansions
                                                   ! used to represent phase functions
double precision, allocatable :: xs(:)             ! the nodes of the kth order Chebyshev grid

! Representation of the phase function for (1) on the interval (0,b)
integer                       :: nints
double precision, allocatable :: ab(:,:)

! Piecewise Chebyshev expansions of the phase function and its first few derivatives
! for evaluation on the interval (-1,1)
double precision, allocatable :: acoefs(:,:)
double precision, allocatable :: apcoefs(:,:)
double precision, allocatable :: appcoefs(:,:)

! the values of the nonoscillatory phase function Psi and its derivatives at 0
double precision              :: aval
double precision              :: apval
double precision              :: appval
double precision              :: apppval

!
! Piecewise Chebyshev expansions of the phase function and its first few derivatives
! for evalution on the interval (1,\infty)
!
end type     prolates_phase_data


! parameters for the nonlinear solver:
!
!    maxiters     = maximum number of Newton iterates
!    ntrapiters   = number of iterations for trapezoidal method
!    dtail        = defines notion of "trailing coefficients"
!                                         
integer, private, parameter              :: maxiters     = 8
integer, private, parameter              :: ntrapiters   = 4
double precision, parameter, private     :: dtail        = 0.50d0
double precision, parameter, private     :: dtailappell  = 0.50d0
double precision, parameter, private     :: dtailappell2 = 0.25d0

contains

subroutine prolates_phase1(chebdata,gamma,chi,phase)
implicit double precision (a-h,o-z)
type(chebexps_data)                                      :: chebdata
type(prolates_phase_data), intent(out)                   :: phase
!
!  Compute the values of the phase function Psi and its first few derivatives
!  at 0 given the values of the parameters gamma and chi in (2).
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
!    (gamma,chi) - the parameters in (2)
!
!  Output parameters:
!    phase - a data structure which will hold the values of the phase
!      function and its first few derivatives upon the completion
!      of this routine
!
!type(chebexps_data)                      :: chebdata
double precision, allocatable            :: vals(:), vals2(:), vals3(:)

! for the solution of Appell's equation on (0,1)
double precision, allocatable            :: ws(:,:), wders(:,:), wder2s(:,:)
double precision, allocatable            :: wder3s(:,:), ab(:,:)

! for the solution of Riccati's equation on the imaginary axis
double precision, allocatable            :: ab_r(:,:)
double precision, allocatable            :: rs(:,:),rders(:,:),rints(:,:)
double complex                           :: ima

eps0           = epsilon(0.0d0)
eps            = 1.0d-12
if (eps0 .lt. 1.0d-16) eps = 1.0d-15
if (eps0 .lt. 1.0d-30) eps = 1.0d-20

if (gamma .ge. 10**5) eps=eps*10

ima            = (0.0d0,1.0d0)
pi             = acos(-1.0d0)
k              = chebdata%k

phase%gamma    = gamma
phase%k        = k
dmu            = 0


!  Solve Riccati's equation along the imaginary axis
call prolates_riccati_tvp(ier,eps,chebdata,gamma,chi,nints_r,ab_r,&
  rs,rders,rints)

if (ier .ne. 0) then
call prini("after riccati_tvp_adap, ier = ",ier)
call prind("gamma   = ",gamma)
call prin2("eps     = ",eps)
call prind("chi     = ",chi)
stop
endif

rval = rs(1,1)
rder = rders(1,1)

!  Compute the values of the phase function and its derivatives at 0 under
!  the change of variable isubst = 0

! apval0   = -rval
! appval0  = 0
! apppval0 = -rder * 2 * rval


! call prolates_convert_to_w(apval0,appval0,apppval0,wa,wpa,wppa)

! wval0   = -1/rval
! wpval0  = -appval0/apval0**2
! wppval0 = 2*appval0**2/apval0**3 - apppval0 / apval0**2



! !  Apply the transformation to the exponential change of variables and then
! !  convert the values of the first three derivatives of \alpha(x) to
! !  the values of w(0), w'(0) and w''(0)

! isubst         = 2
! call prolates_convert_from0(isubst,apval0,appval0,apppval0,apval,appval,apppval)
! call prolates_convert_to_w(apval,appval,apppval,wa,wpa,wppa)


apval   = -rval
appval  = rval
apppval = -2*rder*rval - rval

wa      = -1/rval
wpa     = -1/rval
wppa    = -1/rval + 2*rder/rval

!  Solve Appell's equation going forward using an adaptive Chebyshev spectral
!  method

call prolates_appell(ier,eps,chebdata,gamma,chi,wa,wpa,wppa, &
  nints,phase%ab,ws,wders,wder2s,wder3s)

if (ier .ne. 0) then
call prini("after prolates_appell, ier = ",ier)
call prind("gamma   = ",gamma)
call prin2("eps     = ",eps)
call prind("chi     = ",chi)
stop
endif


phase%nints  = nints
phase%k      = k

!  Construct the expansions of the derivatives of alpha


allocate(phase%acoefs(k,nints))
allocate(phase%apcoefs(k,nints))
allocate(phase%appcoefs(k,nints))
allocate(vals(k),vals2(k),vals3(k))


aval = 0
do int=nints,1,-1
a0       =  phase%ab(1,int)
b0       =  phase%ab(2,int)

vals     = 1.0d0/ws(:,int)
vals2    = aval + (b0-a0)/2 * matmul(chebdata%aintr,vals)
aval     = vals2(1)

phase%acoefs(:,int)   = matmul(chebdata%u,vals2)
phase%apcoefs(:,int)  = matmul(chebdata%u,vals)
vals     = 1.0d0/ws(:,int)**2
vals     = -wders(:,int)*vals
phase%appcoefs(:,int)  = matmul(chebdata%u,vals)

end do

! aval = 0.0d0

! aval     = vals2(1)

! phase%acoefs(:,int)   = matmul(chebdata%u,vals2)
! phase%apcoefs(:,int)  = matmul(chebdata%u,vals)


! end do

phase%aval    = aval
phase%apval   = apval
phase%appval  = appval
phase%apppval = apppval
phase%chi     = chi

end subroutine



subroutine prolates_phase2(chebdata,phase)
implicit double precision (a-h,o-z)
type(chebexps_data)                                      :: chebdata
type(prolates_phase_data), intent(inout)                 :: phase
!
!  Compute the values of the phase function Psi and its first few derivatives
!  at 0 given the values of the parameters gamma and chi in (2).
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
!    (gamma,chi) - the parameters in (2)
!
!  Output parameters:
!    phase - a data structure which will hold the values of the phase
!      function and its first few derivatives upon the completion
!      of this routine
!
!type(chebexps_data)                      :: chebdata
double precision, allocatable            :: vals(:), vals2(:), vals3(:)

! for the solution of Appell's equation on (0,1)
double precision, allocatable            :: ws(:,:), wders(:,:), wder2s(:,:)
double precision, allocatable            :: wder3s(:,:), ab(:,:)

! for the solution of Riccati's equation on the imaginary axis
double precision, allocatable            :: ab_r(:,:)
double precision, allocatable            :: rs(:,:),rders(:,:),rints(:,:)
double complex                           :: ima


eps0           = epsilon(0.0d0)
eps            = 1.0d-12
if (eps0 .lt. 1.0d-16) eps = 1.0d-14
if (eps0 .lt. 1.0d-30) eps = 1.0d-20


if (gamma .gt. 10**5) eps = eps*10

ima        = (0.0d0,1.0d0)
pi         = acos(-1.0d0)
k          = chebdata%k

gamma      = phase%gamma
chi        = phase%chi
aval       = phase%aval
apval      = phase%apval
appval     = phase%appval
apppval    = phase%apppval
call prolates_convert_to_w(apval,appval,apppval,wa,wpa,wppa)

!
!  Solve Appell's equation going forward using an adaptive Chebyshev spectral
!  method

call prolates_appell2(ier,eps,chebdata,gamma,chi,wa,wpa,wppa, &
  nints,phase%ab,ws,wders,wder2s,wder3s)

if (ier .ne. 0) then
call prini("after prolates_appell, ier = ",ier)
call prind("gamma   = ",gamma)
call prin2("eps     = ",eps)
call prind("chi     = ",chi)
stop
endif

phase%nints  = nints
phase%k      = k

!  Construct the expansions of the derivatives of alpha

!  Construct the expansions of alpha and alpha' going forward
if (allocated(phase%acoefs))  deallocate(phase%acoefs)
if (allocated(phase%apcoefs)) deallocate(phase%apcoefs)

allocate(phase%acoefs(k,nints))
allocate(phase%apcoefs(k,nints))
allocate(vals(k),vals2(k))

!
!  Integrate forward if alpha'(x) is not that small on the right
!

b   = phase%ab(2,nints)
val = ws(k,nints)


if (val .lt. 1.0d15) then

do int=1,nints
a0       =  phase%ab(1,int)
b0       =  phase%ab(2,int)

vals     = 1.0d0/ws(:,int)
vals2    = aval + (b0-a0)/2 * matmul(chebdata%aintl,vals)
aval     = vals2(k)

phase%acoefs(:,int)   = matmul(chebdata%u,vals2)
phase%apcoefs(:,int)  = matmul(chebdata%u,vals)

end do

else

!
!  Otherwise integrate backward
!

aval0 = aval
aval  = 0
do int=nints,1,-1
a0       =  phase%ab(1,int)
b0       =  phase%ab(2,int)

vals     = 1.0d0/ws(:,int)
vals2    = aval + (b0-a0)/2 * matmul(chebdata%aintr,vals)
aval     = vals2(1)

phase%acoefs(:,int)   = matmul(chebdata%u,vals2)
phase%apcoefs(:,int)  = matmul(chebdata%u,vals)

end do

endif



end subroutine



subroutine prolates_phase3(chebdata,gamma,chi,aval,apval,appval,apppval)
implicit double precision (a-h,o-z)
type(chebexps_data)                                      :: chebdata
!
!  Given the values of the paramters gamma and chi in (2), return the values of
!  the phase function and its first few derivatives at 0 .
!
!  Input parameters:
!    chebdata - a structure returned by the chebexps routine which, among other
!      things, specifies the order k of the Chebyshev expansions
!    (gamma,chi) - the parameters in (2)
!
!  Output parameters:
!    aval    - the value of Psi(0)
!    apval   - the value of Psi'(0)
!    appval  - the value of Psi''(0)
!    apppval - the value of Psi'''(0)
!

type(prolates_phase_data)           :: phase
call prolates_phase1(chebdata,gamma,chi,phase)
aval    = phase%aval
apval   = phase%apval
appval  = phase%appval
apppval = phase%apppval

end subroutine


subroutine prolates_phase_eval_psi(phase,x,aval,apval)
implicit double precision (a-h,o-z)
type(prolates_phase_data), intent(in)                   :: phase
!
!  Evaluate the phase function Psi and its first derivative at a point x.
!
!  Input parameters:
!    phase - the structure describing the phase function 
!    x - point at which to evaluate the phase function
!
!  Output parameters:
!    aval - the value of the phase function at the point x
!    apval - the derivative of the phase function at x
!    appval - the second derivative of the phase function at x 
!

call chebpw_eval22(phase%nints,phase%ab,phase%k,phase%acoefs,phase%apcoefs,&
  x,aval,apval)

end subroutine


subroutine prolates_phase_eval_ps(phase,t,valps)
implicit double precision (a-h,o-z)
type(prolates_phase_data), intent(in)                   :: phase
!
!  Evaluate the function (4) -- SANS NORMALIZATION CONSTANT -- at a point t in (-1,1).
!
!  Input parameters:
!    phase - the structure describing the phase function 
!    t - a point in the interval (-1,1)
!
!  Output parameters:
!    valps - the value of the function (4) 
!

if (t .lt. 0) then
t1 = -t

x = -log(1-t1)
call  chebpw_eval22(phase%nints,phase%ab,phase%k,phase%acoefs,phase%apcoefs,&
  x,aval,apval)

aval   = 2*phase%aval-aval
dd     = 1/sqrt(1+t1)
valps  = sin(aval)/sqrt(apval)*dd

else
x = -log(1-t)
call  chebpw_eval22(phase%nints,phase%ab,phase%k,phase%acoefs,phase%apcoefs, &
  x,aval,apval)
dd    = 1/sqrt(1+t)

valps  = sin(aval)/sqrt(apval)*dd
!derps  = 

endif


end subroutine

subroutine prolates_phase_eval_psqs(phase,t,valps,valqs)
implicit double precision (a-h,o-z)
type(prolates_phase_data), intent(in)                   :: phase
!
!  Evaluate the function (4) -- SANS NORMALIZATION CONSTANT -- at a point t in (-1,1).
!
!  Input parameters:
!    phase - the structure describing the phase function 
!    t - a point in the interval (-1,1)
!
!  Output parameters:
!    valps - the value of the function (4) 
!

if (t .lt. 0) then
t1 = -t

x = -log(1-t1)
call  chebpw_eval22(phase%nints,phase%ab,phase%k,phase%acoefs,phase%apcoefs,&
  x,aval,apval)

aval   = 2*phase%aval-aval
dd     = 1/sqrt(1+t1)
valps  = sin(aval)/sqrt(apval)*dd
valqs  = cos(aval)/sqrt(apval)*dd

else
x = -log(1-t)
call  chebpw_eval22(phase%nints,phase%ab,phase%k,phase%acoefs,phase%apcoefs, &
  x,aval,apval)
dd    = 1.0d0/sqrt(1+t)

valps  = sin(aval)/sqrt(apval)*dd
valqs  = cos(aval)/sqrt(apval)*dd


endif

end subroutine


subroutine prolates_phase_eval_psder(phase,t,valps,derps)
implicit double precision (a-h,o-z)
type(prolates_phase_data), intent(in)                   :: phase
!
!  Evaluate the function (4) -- SANS NORMALIZATION CONSTANT -- and its derivative
!  at a point t in (-1,1).
!
!  Input parameters:
!    phase - the structure describing the phase function 
!    t - a point in the interval (-1,1)
!
!  Output parameters:
!    valps - the value of the function (4) 
!    derps - the value of the derivative of the function (4)
!

if (t .lt. 0) then
t1 = -t
x = -log(1-t1)
call  chebpw_eval23(phase%nints,phase%ab,phase%k,phase%acoefs,phase%apcoefs,phase%appcoefs,&
  x,aval,apval,appval)

aval  = 2*phase%aval-aval
dd1   = 1.0d0/sqrt(1+t1)
uval  = sin(aval)/sqrt(apval)
uder  = cos(aval)*sqrt(apval) + 0.5d0*sin(aval)*appval/(apval**1.5d0)

valps = uval*dd1
derps = 0.5d0*valps/(1+t1) + uder/(1-t1)*dd1

else

x = -log(1-t)

call  chebpw_eval23(phase%nints,phase%ab,phase%k,phase%acoefs,phase%apcoefs,phase%appcoefs, &
  x,aval,apval,appval)

dd1   = 1.0d0/sqrt(1+t)
uval  = sin(aval)/sqrt(apval)
uder  = cos(aval)*sqrt(apval) - 0.5d0*sin(aval)*appval/(apval**1.5d0)

valps = uval*dd1
derps = -0.5d0*valps/(1+t) + uder/(1-t)*dd1

endif


end subroutine


subroutine prolates_phase_eval(phase,t,valps)
implicit double precision (a-h,o-z)
type(prolates_phase_data), intent(in)                   :: phase
!
!  Evaluate the function (3) at a point s.
!
!  Input parameters:
!    phase - the structure describing the phase function 
!    t - a point in the interval (0,1)
!
!  Output parameters:
!    valps - the value of the function u
!

nints = phase%nints
x     = -log(1-t)

if (x .gt. phase%ab(2,nints)) then
valps = 0
return
endif


call  chebpw_eval22(phase%nints,phase%ab,phase%k,phase%acoefs,phase%apcoefs,&
  x,aval,apval)
dd    = 1/sqrt(1+t)
valps  = sin(aval)/sqrt(apval)*dd

end subroutine

subroutine prolates_phase_eval2(phase,t,valps,aval,apval)
implicit double precision (a-h,o-z)
type(prolates_phase_data), intent(in)                   :: phase
!
!  Evaluate the function (3) at a point s.
!
!  Input parameters:
!    phase - the structure describing the phase function 
!    t - a point in the interval (0,1)
!
!  Output parameters:
!    valps - the value of the function u
!

x     = -log(1-t)
nints = phase%nints
b     = phase%ab(2,nints)

if (x .ge. b) then
valps = 0
aval  = 0
apval = 0
return
endif

call  chebpw_eval22(phase%nints,phase%ab,phase%k,phase%acoefs,phase%apcoefs,&
  x,aval,apval)
dd    = 1/sqrt(1+t)
valps  = sin(aval)/sqrt(abs(apval))*dd

end subroutine

subroutine prolates_phase_eval0(phase,valps,derps)
implicit double precision (a-h,o-z)
type(prolates_phase_data), intent(in)                   :: phase
!
!  Evaluate the function (3) and its derivative at 0.
!
!  Input parameters:
!    phase - the structure describing the phase function 
!    t - a point in the interval (0,1)
!
!  Output parameters:
!    valps - the value of u(0)
!    derps - the derivative u'(0)
!

aval   = phase%aval
apval  = phase%apval
appval = phase%appval

valps  = sin(aval)/sqrt(apval)
derps  = -0.5d0 * sin(aval) * (1.0d0/apval**0.5d0 + appval/apval**1.5d0) + cos(aval)*sqrt(apval)

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine prolates_appell(ier,eps,chebdata,gamma,chi,wa,wpa,wppa, &
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
!  on an interval of the form interval (0,b) using an adaptive Chebyshev spectral method.  
!  The function
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
!    (gamma,chi) - the parameters in (1)
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

dmu      = 0
maxints  = 1000
k        = chebdata%k
nn       = k * dtailappell
dstop    = 1d250

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

a = 0.0d0
b = 1.0d20

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
call qprolates(gamma,chi,k,ts,ps,qs)

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

call prolates_ivp_int(a0,b0,k,chebdata%xs,chebdata%aintl,chebdata%aintl2, &
  chebdata%aintl3,rs,ps,qs,fs,ys0,yders0,yder2s0,yder3s0)

ifstop  = 0
if (abs(ys0(k)) .gt. dstop) then
ifstop=1
endif
if (ifstop .eq. 1) exit

ifsplit = 0

coefs0 = matmul(chebdata%u,ys0)
dd1 = maxval(abs(coefs0))
dd2 = maxval(abs(coefs0(k-nn+1:k)))
if (dd2 .gt. eps*dd1) ifsplit = 1

!print *,a0,b0,dd2/dd1,eps
! coefs0 = matmul(chebdata%u,yders0)
! dd1 = maxval(abs(coefs0))
! dd2 = maxval(abs(coefs0(k-nn+1:k)))

! if (dd2 .gt. eps*dd1) ifsplit = 1
! coefs0 = matmul(chebdata%u,1/ys0)
! dd1 = maxval(abs(coefs0))
! dd2 = maxval(abs(coefs0(k-nn+1:k)))
! if (dd2 .gt. eps*dd1) ifsplit = 1

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



subroutine prolates_ivp_int(a,b,k,xscheb,chebint,chebint2,chebint3,rs,ps,qs,&
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




subroutine prolates_appell2(ier,eps,chebdata,gamma,chi,wa,wpa,wppa, &
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
!  on an interval of the form interval (0,b) using an adaptive Chebyshev spectral method.  
!  The function
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
!    (gamma,chi) - the parameters in (1)
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
eps0     = epsilon(0.0d0)
dstop    = 1d30
dmu      = 0
maxints  = 1000
k        = chebdata%k
nn       = k * dtailappell2

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

a =   0.0d0
b =  30.0d0

nints0 = 4
dd     = log(b-a)
nints0 = ceiling(dd)
dd     = dd/(nints0+0.0d0)

! nints0 = 8
! dd     = log(b-a)/(nints0+0.0d0)


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
call qprolates(gamma,chi,k,ts,ps,qs)

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

call prolates_ivp_int(a0,b0,k,chebdata%xs,chebdata%aintl,chebdata%aintl2, &
  chebdata%aintl3,rs,ps,qs,fs,ys0,yders0,yder2s0,yder3s0)

if (abs(ys0(k)) .gt. dstop) exit

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



subroutine qprolates(gamma,chi,k,xs,vals,ders)
implicit double precision (a-h,o-z)
double precision :: xs(k),vals(k),ders(k)
!
!  Evaluate the coefficient q(x) in the equation (1) and its derivative at
!  a collection of points.
!
!  Input parameters:
!    (c,chi) - the parameters in (1)
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


if (z .gt. 1.0d7)  then

val = (3+0.2d1*chi-0.2d1*gamma**2)/(0.64d2*z**5)+(5+0.4d1*chi-  &
      0.4d1*gamma**2)/(0.64d2*z**4)+(1+chi-  &
      0.1d1 *gamma**2)/(0.8d1*z**3)+(3+0.4d1*chi+  &
      0.12d2*gamma**2)/(0.16d2*z**2)+(1+0.2d1*chi-  &
      0.2d1*gamma**2)/(0.4d1*z)

der = (5*(-3-0.2d1*chi+0.2d1*gamma**2))/(0.64d2*z**5)+(-5-  &
      0.4d1*chi+0.4d1*gamma**2)/(0.16d2*z**4)+(3*(-1-0.1d1*chi+  &
      gamma**2))/(0.8d1*z**3)+(-3-0.4d1*chi-  &
      0.12d2*gamma**2)/(0.8d1*z**2)+(-1-0.2d1*chi+  &
      0.2d1*gamma**2)/(0.4d1*z)

else


val = ((1-0.2d1*z)**(-2)+(0.4d1*gamma**2)/z**2+(2+0.4d1*chi-  &
      0.4d1*gamma**2)/(-1+0.2d1*z))/0.4d1


der  = (0.2d1*(z**3*(chi-z-2*chi*z)+gamma**2*(-1+z)*(-1+2*z)*(1+(-3  &
       +z)*z)))/(z**2*(-1+2*z)**3)

endif

vals(i) = val
ders(i) = der

end do




end subroutine




subroutine prolates_convert_from0(isubst,apval0,appval0,apppval0, &
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


subroutine prolates_convert_to0(isubst,apval,appval,apppval, &
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



subroutine prolates_convert_to_w(apval,appval,apppval,wval,wpval,wppval)
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



subroutine prolates_convert_from_w(wval,wpval,wppval,apval,appval,apppval)
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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  The following code solves the Riccati equation corresponding to the spheroidal
!  wave equation going down the imaginary axis
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine prolates_riccati_tvp(ier,eps,chebdata,gamma,chi,nints,ab,rs,rders,rints)
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
!                 1           chi + gamma^2 t^2 
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
!     (gamma,chi) - the parameters in (1)
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
b         =  1.012
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
call prolates_qaxis(gamma,chi,k,ts0,qs0)

coefs0 = matmul(chebdata%u,qs0)
dd1    = maxval(abs(coefs0))
dd2    = maxval(abs(coefs0(k-ntail+1:k)))
!dd2    = maxval(abs(coefs0(k/2+1:k)))


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
call  prolates_qaxis(gamma,chi,k,ts0,qs0)


if (nintsout .eq. 0) then
rb0 = rb
else
rb0 = rsout(1,nintsout)
endif


! use the trapezoidal rule to construct an initial guess for Newton iterations
call prolates_riccati_trap_tvp(ier,k,ts0,qs0,rb0,rs0,rders0)

if (ier .eq. 1024) then
call prin2("in spheroidal_riccati_tvp, NaN encounted, a0 = ",a0)
call prin2("in spheroidal_riccati_tvp, NaN encounted, b0 = ",b0)
return
endif

! integrate the derivative to obtain the values of rs --- this step is crucial
rs0 = rb0 + matmul(chebdata%aintr*(b0-a0)/2,rders0)

! perform Newton iterations 
do iter=1,maxiters

ps0 = 2 *rs0
fs0 = -rders0 - rs0**2 - qs0
rb0  = 0
call prolates_riccati_linear_tvp(k,a0,b0,ts0,chebdata%aintr,ps0,fs0,rb0,amatr0,delta,deltap)

dd2 = norm2(delta)
dd1 = norm2(rs0)
dd  = dd2/dd1
if (dd .lt. eps0*10) exit

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


subroutine prolates_riccati_trap_tvp(ier,k,ts,qs,rb,rs,rders)
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
eps    = epsilon(0.0d0)
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

do iter=1,ntrapiters
delta = (rhs - r +h/2*r**2 ) / (1 - h*r)
r     = r + delta

! dd1  = abs(r)
! dd2  = abs(delta)
! dd   = dd2/dd1
! if (dd .lt. eps0*10) exit

end do

rs(j)    = r
rders(j) = -r**2 - q0
end do


end subroutine


subroutine prolates_riccati_linear_tvp(k,a,b,ts,chebint,qs,fs,rb,amatr,rs,rders)
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


subroutine prolates_qaxis(gamma,chi,k,ts,qs)
implicit double precision (a-h,o-z)
integer                              :: k
double precision                     :: ts(k)
double precision                     :: qs(k)

qs = - ( 1.0d0) / (1+ts**2)**2 - (chi + gamma**2 * ts**2) / ( 1 + ts**2)

end subroutine

end module
