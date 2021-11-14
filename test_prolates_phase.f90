program test_prolates_phase
use utils
use prolates
use chebyshev
use prolates_phase

implicit double precision (a-h,o-z)

integer                             :: nints, k
type(chebexps_data)                 :: chebdata
type(prolates_phase_data)           :: phase
double precision, allocatable       :: coefsps(:)

eps0    = epsilon(0.0d0)
k       = 30
pi      = acos(-1.0d0)
inormal = 1

k = 30
call chebexps(k,chebdata)

gamma = 100d0
n     = gamma/2
dnu   = n

!  Use the Osipov-Rokhlin algorithm to construct coefficient expansions
call elapsed(t1)
do ii=1,10
call prolates_integer(inormal,gamma,n,chi,n1,n2,coefsps)
end do
call elapsed(t2)
dtime = (t2-t1)/10

call prolates_integer_evalder(n,n1,n2,coefsps,0.0d0,valps0,derps0)
call prini("n1 = ",n1)
call prini("n2 = ",n2)
call prin2("prolates_integer time = ",dtime)
!call prind("after prolates_integer, chi = ",chi)
!call prin2("coefps = ",coefsps)


!  Calculate the values of the phase function and its first few derivatives at 0
call elapsed(t1)
do ii=1,10
call prolates_phase1(chebdata,gamma,chi,phase)
end do
call elapsed(t2)
dtime = (t2-t1)/10
call prin2("prolates_phase1 time = ",dtime)

! Calculate the normalization constant
call prolates_phase_eval0(phase,valps,derps)
call prolates_lege0(dnu,valps0,derps0)

if ( mod(n,2) .eq. 1 ) then
dnorm = derps0 / derps
else
dnorm = valps0 / valps
endif


call elapsed(t1)
do ii=1,10
call prolates_phase2(chebdata,phase)
end do
call elapsed(t2)
dtime = (t2-t1)/10
call prin2("prolates_phase2 time  = ",dtime)

dd   = -phase%aval*2/pi
derr = abs(dd-(n+1))
if (n .gt. 0) derr = derr / abs(n+1)
call prin2("relative error in the number of zeros = ",derr)



! check the error
t = 0.10d0
call elapsed(t1)
call prolates_phase_eval_ps(phase,t,valps)
valps = valps*dnorm
call elapsed(t2)
call prin2("phase eval time = ",t2-t1)


call elapsed(t1)
call prolates_integer_eval(n,n1,n2,coefsps,t,valps0)
call elapsed(t2)
call prin2("Legendre expansion eval time = ",t2-t1)


!derr = abs((valps-valps0)/valps0)
!call prin2("relative error in evaluation = ",derr)
derr = abs((valps-valps0))
call prin2("absolute error in evaluation = ",derr)


end program
