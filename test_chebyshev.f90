module test_chebyshev_subroutines

use utils
use chebyshev

contains


subroutine testfun(ts,vals)
implicit double precision (a-h,o-z)
double precision, intent(in)     :: ts(:)
double precision, intent(out)    :: vals(:)
vals = cos(ts) + sin(ts)**2 + ts**2 - 1
end subroutine


end module



program test_chebyshev

use utils
use chebyshev
use test_chebyshev_subroutines
use iso_c_binding

implicit double precision (a-h,o-z)

type(chebexps_data)           :: chebdata
double precision, allocatable :: ab(:,:),vals0(:,:),ts(:),vals(:),vals1(:)
double precision, allocatable :: z(:),w(:,:),coefs(:)


k = 20
call chebexps(k,chebdata)

!
!  Build a piecewise discretization scheme
!

a = 0.0d0
b = 1.0d0

nints = 30
allocate(ab(2,nints))
do i=1,nints
ab(1,i) = (i-1.0d0)/(nints+0.0d0)
ab(2,i) = (i+0.0d0)/(nints+.0d0)
end do


!
!  Sample a test function on this grid
!

allocate(vals0(k,nints))

do int=1,nints
a = ab(1,int)
b = ab(2,int)
do i=1,k
x = chebdata%xs(i) * (b-a)/2 + (b+a)/2
vals0(i,int) = cos(133*x)
end do
end do

m = 100
allocate(ts(m),vals(m))
do i=1,m
ts(i) = (i-1.0d0)/(m-1.0d0)
end do

!
!  Interpolate it to the equispaced grid
!

call elapsed(t1)
call chebpw_aint(chebdata,nints,ab,m,ts,vals0,vals)
call elapsed(t2)
call prin2("aint average time = ",(t2-t1)/m)

errmax = 0
do i=1,m
t = ts(i)
errabs = abs(vals(i)-cos(133*t))
errmax = max(errabs,errmaX)
end do
call prin2("aint errmax = ",errmax)

!
!  Test the transpose interpolation routine
!

allocate(z(m),w(k,nints))

do i=1,m
call random_number(z(i))
end do

dsum1 = 0.0d0
do i=1,m
t = ts(i)
dsum1 = dsum1 + z(i)*cos(133*t)
end do

call elapsed(t1)
call chebpw_aintt(chebdata,nints,ab,m,ts,z,w)
call elapsed(t2)
call prin2("aintl average time = ",(t2-t1)/m)

dsum2 = 0.0d0
do int=1,nints
do i=1,k
t = ts(i)
dsum2 = dsum2 + w(i,int)*vals0(i,int)
end do
end do

call prin2("aintl error = ",dsum1-dsum2)

!
!  Test integration
!



end program
