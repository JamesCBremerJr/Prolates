module test_linearode_functions
use utils
use chebyshev
use linearode

contains

subroutine fun3(k,ts,rs,ps,qs,fs,userptr)
implicit double precision (a-h,o-z)
integer                       :: k
double precision              :: ts(k),rs(k),ps(k),qs(k),fs(k)
type(c_ptr)                   :: userptr

rs = 1.0d0/(.1d0+ts)
ps = ts**2
qs = ts**3
fs = cos(ts)

end subroutine

subroutine fun2(k,ts,ps,qs,fs,userptr)
implicit double precision (a-h,o-z)
integer                       :: k
double precision              :: ts(k),rs(k),ps(k),qs(k),fs(k)
type(c_ptr)                   :: userptr

ps = cos(7*ts)
qs = sin(13*ts)
fs = ts

! ps = 0
! qs = 1000
! fs = ts

end subroutine


end module




program test_linearode
use test_linearode_functions
implicit double precision (a-h,o-z)
type(c_ptr)                   :: userptr
type(chebexps_data)           :: chebdata
double precision, allocatable :: ab(:,:), ys(:,:),yders(:,:),yder2s(:,:),yder3s(:,:)


k    = 30
call chebexps(k,chebdata)
eps  = 1.0d-12

!

!  Test the second order solve routines
!

a    = 1
b    = 2

ya   = 1
ypa  = 2

! ya   = 0.001d0
! ypa  = 0.001d0

val0 = 3.8096154401135364162314338163416987d0

call solve2_ivp_adap(ier,eps,a,b,chebdata,nints,ab,ys,yders,yder2s,&
  fun2,ya,ypa,userptr)

call prini("after solve2_ivp_adap, ier = ",ier)
call prin2("ab = ",ab)

errrel = abs(ys(k,nints)-val0)/abs(val0)
call prin2("errrel = ",errrel)
call prina("")


yb  = ys(k,nints)
ypb = yders(k,nints)

call solve2_tvp_adap(ier,eps,a,b,chebdata,nints,ab,ys,yders,yder2s, &
  fun2,yb,ypb,userptr)

call prini("after solve2_ivp_adap, ier = ",ier)
call prin2("ab = ",ab)


errrel = abs(ya - ys(1,1))/abs(ya)
call prin2("errrel = ",errrel)
call prina("")


!
!  Test the third order solvers
!

a    = 0
b    = 1
ya   = 1
ypa  = 2
yppa = 2


val0 = 3.37346574656496314475897879523967605d0

call solve3_ivp_adap(ier,eps,a,b,chebdata,nints,ab,ys,yders,yder2s,yder3s,&
  fun3,ya,ypa,yppa,userptr)

call prini("after solve3_ivp_adap, ier = ",ier)
call prin2("ab = ",ab)

errrel = abs(ys(k,nints)-val0)/abs(val0)
call prin2("errrel = ",errrel)
call prina("")


yb   = ys(k,nints)
ypb  = yders(k,nints)
yppb = yder2s(k,nints)

call solve3_tvp_adap(ier,eps,a,b,chebdata,nints,ab,ys,yders,yder2s,yder3s,&
  fun3,yb,ypb,yppb,userptr)

call prini("after solve3_tvp_adap, ier = ",ier)
call prin2("ab = ",ab)

errrel = abs(ys(1,1)-ya)/abs(ya)
call prin2("errrel = ",errrel)

end program


