program test_prolates_window
use utils
use chebyshev
use riccati
use linearode
use prolates
use prolates_window

implicit double precision (a-h,o-z)

type(chebexps_data)                 :: chebdata
double precision, allocatable       :: coefsps(:),coefsqs(:)
double complex                      :: clambda




k   = 30
call chebexps(k,chebdata)

eps    = 1.0d-12
pi     = acos(-1.0d0)

!
!  Test chi_n > c^2 regime
!

c      = 500.0d0
n      = c
dnu    = n
isubst = 0

call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)
call prolates_integer_coefficients(c,n,ncoefs,coefsps,coefsqs,c1,c2)
call prolates_integer_evalder(n,ncoefs,coefsps,coefsqs,0.0d0,valps,derps,valqs,derqs)

wa      = ((c1*valps)**2   + (c2*valqs)**2)
wpa     = 2*c1*valps*derps + 2*c2*valqs*derqs
apval0  = 1/wa

call elapsed(t1)
call prolates_window_avals(c,chi,isubst,aval,apval,appval,apppval)
call elapsed(t2)


time   = t2-t1
dnu0   = -aval/(pi/2)-1
errrel = abs((dnu-dnu0)/dnu0)
errrel2  = abs((apval-apval0)/apval0)

write (*,"(A,D12.6,A,D12.6,A,D12.6,A,D12.6,A,D12.6)") "prolates_windows_avals:  c= ",&
  c, ", chi = ",chi,", errrel = ",errrel,", errrel2 = ",errrel2, ",  time = ",time

c      = 5000.0d0
n      = c
dnu    = n
isubst = 0

call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)


call prolates_integer_coefficients(c,n,ncoefs,coefsps,coefsqs,c1,c2)
call prolates_integer_evalder(n,ncoefs,coefsps,coefsqs,0.0d0,valps,derps,valqs,derqs)

wa      = ((c1*valps)**2   + (c2*valqs)**2)
wpa     = 2*c1*valps*derps + 2*c2*valqs*derqs
apval0  = 1/wa

call elapsed(t1)
call prolates_window_avals(c,chi,isubst,aval,apval,appval,apppval)
call elapsed(t2)


time   = t2-t1
dnu0   = -aval/(pi/2)-1
errrel = abs((dnu-dnu0)/dnu0)
errrel2  = abs((apval-apval0)/apval0)

write (*,"(A,D12.6,A,D12.6,A,D12.6,A,D12.6,A,D12.6)") "prolates_windows_avals:  c= ",&
  c, ", chi = ",chi,", errrel = ",errrel,", errrel2 = ",errrel2, ",  time = ",time


c      = 50000.0d0
n      = c
dnu    = n
isubst = 0

call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)


call prolates_integer_coefficients(c,n,ncoefs,coefsps,coefsqs,c1,c2)
call prolates_integer_evalder(n,ncoefs,coefsps,coefsqs,0.0d0,valps,derps,valqs,derqs)

wa      = ((c1*valps)**2   + (c2*valqs)**2)
wpa     = 2*c1*valps*derps + 2*c2*valqs*derqs
apval0  = 1/wa

call elapsed(t1)
call prolates_window_avals(c,chi,isubst,aval,apval,appval,apppval)
call elapsed(t2)


time   = t2-t1
dnu0   = -aval/(pi/2)-1
errrel = abs((dnu-dnu0)/dnu0)
errrel2  = abs((apval-apval0)/apval0)

write (*,"(A,D12.6,A,D12.6,A,D12.6,A,D12.6,A,D12.6)") "prolates_windows_avals:  c= ",&
  c, ", chi = ",chi,", errrel = ",errrel,", errrel2 = ",errrel2, ",  time = ",time

c      = 500000.0d0
n      = c
dnu    = n
isubst = 0

call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)


call prolates_integer_coefficients(c,n,ncoefs,coefsps,coefsqs,c1,c2)
call prolates_integer_evalder(n,ncoefs,coefsps,coefsqs,0.0d0,valps,derps,valqs,derqs)

wa      = ((c1*valps)**2   + (c2*valqs)**2)
wpa     = 2*c1*valps*derps + 2*c2*valqs*derqs
apval0  = 1/wa

call elapsed(t1)
call prolates_window_avals(c,chi,isubst,aval,apval,appval,apppval)
call elapsed(t2)


time   = t2-t1
dnu0   = -aval/(pi/2)-1
errrel = abs((dnu-dnu0)/dnu0)
errrel2  = abs((apval-apval0)/apval0)

write (*,"(A,D12.6,A,D12.6,A,D12.6,A,D12.6,A,D12.6)") "prolates_windows_avals:  c= ",&
  c, ", chi = ",chi,", errrel = ",errrel,", errrel2 = ",errrel2, ",  time = ",time


c      = 5000000.0d0
n      = c
dnu    = n
isubst = 0

call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)


call prolates_integer_coefficients(c,n,ncoefs,coefsps,coefsqs,c1,c2)
call prolates_integer_evalder(n,ncoefs,coefsps,coefsqs,0.0d0,valps,derps,valqs,derqs)

wa      = ((c1*valps)**2   + (c2*valqs)**2)
wpa     = 2*c1*valps*derps + 2*c2*valqs*derqs
apval0  = 1/wa

call elapsed(t1)
call prolates_window_avals(c,chi,isubst,aval,apval,appval,apppval)
call elapsed(t2)


time   = t2-t1
dnu0   = -aval/(pi/2)-1
errrel = abs((dnu-dnu0)/dnu0)
errrel2  = abs((apval-apval0)/apval0)

write (*,"(A,D12.6,A,D12.6,A,D12.6,A,D12.6,A,D12.6)") "prolates_windows_avals:  c= ",&
  c, ", chi = ",chi,", errrel = ",errrel,", errrel2 = ",errrel2, ",  time = ",time

!
!   Test small n and medium sized c 
!

c = 250.0d0
do i=1,100

n      = 100
dnu    = n
isubst = 0

call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)

call elapsed(t1)
call prolates_window_avals(c,chi,isubst,aval,apval,appval,apppval)
call elapsed(t2)

time   = t2-t1
dnu0   = -aval/(pi/2)-1
errrel = abs((dnu-dnu0)/dnu0)

write (*,"(A,D12.6,A,D12.6,A,D12.6,A,D12.6)") "prolates_windows_avals:  c= ",&
  c, ", chi = ",chi,", errrel = ",errrel, ",  time = ",time

c      = c + 100.0d0
end do

c      = 50000000.0d0
n      = c
dnu    = n
isubst = 0

call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)


call prolates_integer_coefficients(c,n,ncoefs,coefsps,coefsqs,c1,c2)
call prolates_integer_evalder(n,ncoefs,coefsps,coefsqs,0.0d0,valps,derps,valqs,derqs)

wa      = ((c1*valps)**2   + (c2*valqs)**2)
wpa     = 2*c1*valps*derps + 2*c2*valqs*derqs
apval0  = 1/wa

call elapsed(t1)
call prolates_window_avals(c,chi,isubst,aval,apval,appval,apppval)
call elapsed(t2)


time   = t2-t1
dnu0   = -aval/(pi/2)-1
errrel = abs((dnu-dnu0)/dnu0)
errrel2  = abs((apval-apval0)/apval0)

write (*,"(A,D12.6,A,D12.6,A,D12.6,A,D12.6,A,D12.6)") "prolates_windows_avals:  c= ",&
  c, ", chi = ",chi,", errrel = ",errrel,", errrel2 = ",errrel2, ",  time = ",time

!
!   Test small n and medium sized c 
!

c = 250.0d0
do i=1,100

n      = 100
dnu    = n
isubst = 0

call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)

call elapsed(t1)
call prolates_window_avals(c,chi,isubst,aval,apval,appval,apppval)
call elapsed(t2)

time   = t2-t1
dnu0   = -aval/(pi/2)-1
errrel = abs((dnu-dnu0)/dnu0)

write (*,"(A,D12.6,A,D12.6,A,D12.6,A,D12.6)") "prolates_windows_avals:  c= ",&
  c, ", chi = ",chi,", errrel = ",errrel, ",  time = ",time

c      = c + 100.0d0
end do

!
!  Large c and small n
!

c = 1000000.0d0

do i=1,100

n      = 100+i
dnu    = n
isubst = 0

call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)

call elapsed(t1)
call prolates_window_avals(c,chi,isubst,aval,apval,appval,apppval)
call elapsed(t2)

time   = t2-t1
dnu0   = -aval/(pi/2)-1
errrel = abs((dnu-dnu0)/dnu0)

write (*,"(A,D12.6,A,D12.6,A,D12.6,A,D12.6)") "prolates_windows_avals:  c= ",&
  c, ", chi = ",chi,", errrel = ",errrel, ",  time = ",time

c      = c + 50000

end do


!
!  Large c and large n
!



c      = 1000000.0d0
n      = c

dnu0   = n
isubst = 0


ifroots = 0


call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)
call prolates_integer_coefficients(c,n,ncoefs,coefsps,coefsqs,c1,c2)
call prolates_integer_evalder(n,ncoefs,coefsps,coefsqs,0.0d0,valps,derps,valqs,derqs)

wa      = ((c1*valps)**2   + (c2*valqs)**2)
wpa     = 2*c1*valps*derps + 2*c2*valqs*derqs
apval0  = 1/wa


call elapsed(t1)
call prolates_window_avals(c,chi,isubst,aval,apval,appval,apppval)
call elapsed(t2)

dnu   = -2/pi*aval-1

! print *,dnu0,dnu,(dnu-dnu0)/dnu0


time   = t2-t1
dnu0   = -2/pi*aval-1
errrel = abs((dnu-dnu0)/dnu0)
errrel2  = abs((apval-apval0)/apval0)

write (*,"(A,D12.6,A,D12.6,A,D12.6,A,D12.6,A,D12.6)") "prolates_windows_avals:  c= ",&
  c, ", chi = ",chi,", errrel = ",errrel,", errrel2 = ",errrel2, ",  time = ",time


end program
