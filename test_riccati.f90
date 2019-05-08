module test_riccati_functions
use utils
use chebyshev
use riccati

contains


subroutine test_riccati_trap_ivp()
implicit double precision (a-h,o-z)

type(chebexps_data)           :: chebdata
double precision, allocatable :: ts(:)
double complex, allocatable   :: rs(:),rders(:),qs(:)
double complex                :: ra, val0, ima

ima = (0.0d0,1.0d0)
k   = 2

print *,"test_riccati_trap_ivp"
print *,"-------------------------------------------------------"

do l=1,8

call chebexps(k,chebdata)

allocate(ts(k),rs(k),rders(k),qs(k))

a = 0.0d0
b = 0.5d0

ts = (b-a)/2*chebdata%xs + (b+a)/2
qs = ts**2*ima
fs = 0
ra = (0.0d0,1.0d0)

call elapsed(t1)
call riccati_trap_ivp(ier,k,ts,qs,ra,rs,rders)
call elapsed(t2)

val0 = (0.3924233525020934d0,0.762954025722340d0)
dd2  = abs(rs(k)-val0)

errrel = dd2 / abs(val0)

if (l .gt. 1) then
print "('   k = ',I3.3,'  errrel = ',D12.5,'  rate = ',F12.5)", k,errrel,dd1/dd2
else
! print "('   k = ',I3.3,'  errrel = ',D12.5)", k,errel
endif

!else
!print *,k,abs(rs(k)-val0)
!endif

dd1  = dd2
deallocate(ts,rs,rders,qs)

k = k *2
end do

print *,"-------------------------------------------------------"
print *,""

end subroutine



subroutine test_riccati_trap_tvp()
implicit double precision (a-h,o-z)

type(chebexps_data)           :: chebdata
double precision, allocatable :: ts(:)
double complex, allocatable   :: rs(:),rders(:),qs(:)
double complex                :: rb, val0, ima

ima = (0.0d0,1.0d0)
k   = 2

print *,"test_riccati_trap_tvp"
print *,"-------------------------------------------------------"


val0 = (-0.495424202731495715803372685453516624d0, 0.261968071845953165029871942343462586d0)

do l=1,8

call chebexps(k,chebdata)

allocate(ts(k),rs(k),rders(k),qs(k))

a = 0.0d0
b = 1.0d0

ts = (b-a)/2*chebdata%xs + (b+a)/2
qs = ts**2
rb =-1.0d0+ima

call riccati_trap_tvp(ier,k,ts,qs,rb,rs,rders)

dd2  = abs(rs(1)-val0)
errrel = dd2 / abs(val0)

if (l .gt. 1) then
print "('   k = ',I3.3,'  errrel = ',D12.5,'  rate = ',F12.5)", k,errrel,dd1/dd2
else
! print "('   k = ',I3.3,'  errrel = ',D12.5)", k,errel
endif


dd1  = dd2
deallocate(ts,rs,rders,qs)

k = k *2
 end do

print *,"-------------------------------------------------------"
print *,""

end subroutine

subroutine test_riccati_linear_ivp()
implicit double precision (a-h,o-z)

type(chebexps_data)           :: chebdata
double precision, allocatable :: ts(:)
double complex, allocatable   :: qs(:),fs(:),amatr(:,:),rs(:),rders(:)
double complex                :: ra,ima,val0

ima = (0.0d0,1.0d0)

print *,"test_riccati_linear_ivp"
print *,"-------------------------------------------------------"

do k=3,16
call chebexps(k,chebdata)

allocate(qs(k),fs(k),ts(k),amatr(k,k),rs(k),rders(k))

a  = 1
b  = 2

ra = (0.0d0,1.0d0)
ts = (b-a)/2 * chebdata%xs + (b+a)/2
fs = cos(ts)
qs = ts**2 + ima * ts

call riccati_linear_ivp(k,a,b,ts,chebdata%aintl,qs,fs,ra,amatr,rs,rders)

val0 = (0.04806701969293754358422759737806224d0,-0.000059058819605814648115588368041328d0)
errrel = abs(val0 - rs(k)) / abs(val0)

print "('   k = ',I3.3,'  errrel = ',D12.5)", k,errrel


deallocate(qs,fs,ts,amatr,rs,rders)
end do

print *,"-------------------------------------------------------"
print *,""

end subroutine



subroutine test_riccati_linear_tvp()
implicit double precision (a-h,o-z)

type(chebexps_data)           :: chebdata
double precision, allocatable :: ts(:)
double complex, allocatable   :: qs(:),fs(:),amatr(:,:),rs(:),rders(:)
double complex                :: rb,ima,val0

ima = (0.0d0,1.0d0)

print *,"test_riccati_linear_tvp"
print *,"-------------------------------------------------------"

val0 = (-10.3220966428745140173050004969428548d0,1.2350655417296033364529544977684766d0)

do k=3,16
call chebexps(k,chebdata)
allocate(qs(k),fs(k),ts(k),amatr(k,k),rs(k),rders(k))

a  = 1
b  = 2

rb = (0.0d0,1.0d0)
ts = (b-a)/2 * chebdata%xs + (b+a)/2
fs = cos(ts)
qs = ts**2 + ima * ts

call riccati_linear_tvp(k,a,b,ts,chebdata%aintr,qs,fs,rb,amatr,rs,rders)
errrel = abs(val0 - rs(1)) / abs(val0)
print "('   k = ',I3.3,'  errrel = ',D12.5)", k,errrel

deallocate(qs,fs,ts,amatr,rs,rders)
end do

print *,"-------------------------------------------------------"
print *,""

end subroutine



subroutine odefun(k,ts,qs,userptr)
implicit double precision (a-h,o-z)
integer          :: k
double precision :: ts(k)
double complex   :: qs(k)
type(c_ptr)      :: userptr
qs = 1.0d0/(.1d0+ts**2)
end subroutine

end module



program test_riccati
use utils
use chebyshev
use riccati
use iso_c_binding
use test_riccati_functions
implicit double precision (a-h,o-z)

type(chebexps_data)           :: chebdata
double precision, allocatable :: ts(:)
double precision, allocatable :: ab(:,:)
double complex, allocatable   :: rs(:,:),rders(:,:)
double complex                :: ra,ima,val0
type(c_ptr)                   :: userptr

! call test_riccati_linear_ivp()
! call test_riccati_linear_tvp()
! call test_riccati_trap_ivp()
! call test_riccati_trap_tvp()


ima = (0.0d0,1.0d0)
k   = 30
call chebexps(k,chebdata)


eps = epsilon(0.0d0)*100

val0 = (0.389289274652975199864315422821480130d0,0.279315229612124128900676289546634409d0)


ra  = 1.0d0 + ima

a   = 1.0d0
b   = 2.0d0

call riccati_ivp_adap(ier,eps,chebdata,a,b,odefun,ra,nints,ab,rs,rders,userptr)

if (ier .ne. 0) then
call prini("after riccati_ivp_adap, ier = ",ier)
stop
endif

errrel = abs(rs(k,nints)-val0)/abs(val0)
call prin2("riccati_ivp_adap relative error =",errrel)

!!! now solve the same system backwards and check the initial value
call riccati_tvp_adap(ier,eps,chebdata,a,b,odefun,val0,nints,ab,rs,rders,userptr)

if (ier .ne. 0) then
call prini("after riccati_tvp_adap, ier = ",ier)
stop
endif

errrel = abs(rs(1,1)-ra)/abs(ra)
call prin2("riccati_tvp_adap relative error =",errrel)

end program
