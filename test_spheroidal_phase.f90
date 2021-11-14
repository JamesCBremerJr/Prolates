program test_spheroidal_phase
use utils
use spheroidal
use spheroidal_phase

implicit double precision (a-h,o-z)

integer                             :: nints, k
type(spheroidal_phase_data)         :: phase
double precision, allocatable       :: coefsps(:)

k      = 30
pi     = acos(-1.0d0)

! gamma  = 1000000.0d0
! gamma  =  983292.0d0
gamma  = 10000
n      = 0.50d0*gamma
m      = 0

dnu    = n
dmu    = m


!  Use the Osipov-Rokhlin algorithm to construct coefficient expansions
call elapsed(t1)
call spheroidal_integer(gamma,n,m,chi,n1,n2,coefsps)
call elapsed(t2)

call prini("n1 = ",n1)
call prini("n2 = ",n2)
call prin2("spheroidal_integer time = ",t2-t1)
call prind("after spheriodal_integer, chi = ",chi)
!call prin2("coefps = ",coefsps)

!  Calculate the values of the phase function and its first few derivatives at 0
call elapsed(t1)
call spheroidal_phase1(gamma,dmu,chi,phase)
call elapsed(t2)
call prin2("spheroidal_phase1 time = ",t2-t1)

call elapsed(t1)
call spheroidal_phase2(phase)
call elapsed(t2)
call prin2("spheroidal_phase2 time = ",t2-t1)


dd =-phase%aval*2/pi-1 
derr = abs(dd-(n-m))
call prin2("relative error in the number of zeros = ",derr/abs(n-m))

! Compute the normalization constant
t = 0.001d0
call spheroidal_phase_eval(phase,t,valps,valqs)
call spheroidal_integer_eval(n,m,n1,n2,coefsps,t,valps0)
dnorm = valps0/valps


! check the error
t = 0.1d0
call elapsed(t1)
call spheroidal_phase_eval(phase,t,valps,valqs)
call elapsed(t2)
call prin2("phase eval time = ",t2-t1)


valps = valps*dnorm
call elapsed(t1)
call spheroidal_integer_eval(n,m,n1,n2,coefsps,t,valps0)
call elapsed(t2)
call prin2("Legendre expansion eval time = ",t2-t1)


valps  = abs(valps)
valps0 = abs(valps0)
derr   = abs((valps-valps0)/valps0)
!derr   = abs(valps-valps0)
call prin2("relative error in evaluation = ",derr)

! ! Construct the phase function again and test again
! call elapsed(t1)
! call spheroidal_phase2(phase)
! call elapsed(t2)
! call prin2("spheroidal_phase2 time = ",t2-t1)

t = 0.3d0
call spheroidal_phase_eval(phase,t,valps,valqs)
valps = valps*dnorm
call spheroidal_integer_eval(n,m,n1,n2,coefsps,t,valps0)
valps  = abs(valps)
valps0 = abs(valps0)
derr   = abs((valps-valps0)/valps0)
call prin2("relative error in evaluation = ",derr)


end program
