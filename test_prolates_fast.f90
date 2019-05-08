program test_prolates_fast
use utils
use chebyshev
use tensor
use prolates
use prolates_fast
implicit double precision (a-h,o-z)

type(chebexps_data)                :: chebdata
double precision, allocatable      :: coefsps(:),coefsqs(:)
double precision, allocatable      :: ab(:,:),acoefs(:,:),apcoefs(:,:),appcoefs(:,:)
double complex                     :: clambda,ima
double precision, allocatable      :: ts(:),roots(:),vals(:),vals0(:)
!double complex, allocatable        :: vals(:),vals0(:)

ima = (0.0d0,1.0d0)
pi  = acos(-1.0d0)
k   = 30
call chebexps(k,chebdata)

! c  = 1000000.0d0
! c  = 525333.0d0
c  = 1000000

call prolexps_eta_range(c,eta1,eta2,eta3,eta4)
n  = (eta1+eta2)/2
n  = (eta1+eta2)/2

call prind("c = ",c)
call prini("n = ",n)

call elapsed(t1)
call prolates_integer(c,n,chi0,clambda,ncoefs,coefsps,coefsqs)
call elapsed(t2)
dtime1 = t2-t1

! find the turning point
a      =  0.0d0
b      =  1.0d0
call qprolates_find_roots(0,c,chi0,a,b,nroots,roots)
call prin2("roots = ",roots)

if (nroots .gt. 0) then
b      = roots(1)
else
b      = 0.99d0
endif

b = min(.99d0,b)


call elapsed(t1)
call prolates_fast_phase(chebdata,c,n,chi,nints,ab,acoefs,apcoefs,appcoefs)
call elapsed(t2)
dtime2 = t2-t1

call prini("nints = ",nints)
call prin2("ab = ",ab)

call prini("nints = ",nints)
call prin2("prolates_integer time = ",dtime1)
call prin2("prolates_phase   time = ",dtime2)
ratio = dtime1/dtime2
call prind("speedup ratio                  = ",ratio)

!
!  Normalize PS_n(x) so that PS_n(0) = P_n(0) if n is even and PS_n'(0) = P_n'(0) if
!  n is odd
!

call prolates_integer_evalder(n,ncoefs,coefsps,coefsqs,0.0d0,valps0,derps0,valqs0,derqs0)
call prolates_fast_legendre0(n,valp,derp,valq,derq)

if ( mod(n,2) .eq. 0) then
coefsps = coefsps * valp/valps0
cc2     = valp/valps0
else
coefsps = coefsps * derp/derps0
cc2     = derp/derps0
endif

!
!  Normalize the fast_prolate functions as well
!
t = 0.0d0
call prolates_integer_evalder(n,ncoefs,coefsps,coefsqs,t,valps0,derps0,valqs0,derqs0)
call prolates_fast_evalder(k,nints,ab,acoefs,apcoefs,appcoefs, &
  0.0d0,valps00,derps00,valqs00,derqs00)

if ( mod(n,2) .eq. 0) then
cc1 = valp/valps00
else
cc1 = derp/derps00
endif
call prina("")


!
!  Do the speed trial / accuracy check
!

call prina("")

nn = 100
allocate(ts(nn),vals(nn),vals0(nn))

do i=1,nn
ts(i) =  a + (b-a) * (i-1.0d0)/(nn-1.0d0)
end do

call elapsed(t1)
do i=1,nn
call prolates_integer_eval0(n,ncoefs,coefsps,ts(i),valps)
vals0(i) = valps
end do
call elapsed(t2)
dtime1 = (t2-t1)/nn
call prin2("prolates_integer_eval average time = ",dtime1)

call elapsed(t1)
do i=1,nn
call prolates_fast_eval0(k,nints,ab,acoefs,apcoefs,ts(i),valps)
vals(i) = cc1*valps
end do
call elapsed(t2)
dtime2 = (t2-t1)/nn
call prin2("prolates_fast_eval average time = ",dtime2)
ratio = dtime1/dtime2
call prind("speedup ratio                  = ",ratio)

errabs = maxval(abs(vals-vals0))
call prin2("maximum absolute error = ",errabs)


! call prina("")
! t = 1-exp(-ab(2,nints))
! call prolates_integer_eval0(n,ncoefs,coefsps,t,valps0)
! call prolates_fast_eval0(c,chi,k,nints,ab,acoefs,apcoefs,t,valps)
! valps = valps * cc1

! ! val0 = 0.0340618082915039463880757889561099598d0
! ! print *,val0
! ! stop

! ! print *,c,n,t
! print *,t
! print *,valps0
! print *,valps

end program
