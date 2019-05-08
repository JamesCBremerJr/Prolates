module test_tensor_functions

implicit double precision (a-h,o-z)

contains

subroutine testfun(x,y,val)
implicit double precision (a-h,o-z)
double precision, intent(in)  :: x,y
double precision, intent(out) :: val

val = (1+y**2-y)*cos(x**2+y**2+1)

end subroutine

subroutine testfun2(x,y,val,derx,dery)
implicit double precision (a-h,o-z)
double precision, intent(in)  :: x,y
double precision, intent(out) :: val
val = y**2*x
derx = y**2
dery = 2*y*x
end subroutine

end module

program test_tensor

use utils
use chebyshev
use tensor
use test_tensor_functions

implicit double precision (a-h,o-z)

double precision, allocatable :: xs(:),ys(:),u(:,:)
double precision, allocatable :: vals(:,:),coefs(:)

double precision, allocatable  :: ab(:,:),cd(:,:),coefs2(:,:,:)

!double precision, allocatable :: ab(:,:),vals2(:,:,:),coefs2(:,:)
integer*1, allocatable        :: compis(:),compjs(:)
double precision, allocatable :: compcoefs(:)

!
!  Construct the Chebyshev nodes and u matrix
!

n = 21
call chebpts(n,xs)
call elapsed(t1)
call tensor_umatrix(n,ncoefs,u)
call elapsed(t2)

call prin2("tensor_umatrix time = ",t2-t1)

!
!  Test an expansion
!

allocate(vals(n,n),coefs(ncoefs))

a =  0.11d0
b =  1.21d0

c = -0.23d0
d =  1.00d0

do i=1,n
x = xs(i) * (b-a)/2 + (b+a)/2
do j=1,n
y = xs(j) * (d-c)/2 + (c+d)/2
call testfun(x,y,val)
vals(i,j) = val
end do
end do

call tensor_coefs(n,u,vals,coefs)

!call prin2("vals = ",vals)
!call prin2("coefs = ",coefs)


x = 0.1d0
y = 0.5d0
call tensor_eval(n,a,b,c,d,coefs,x,y,val)
call testfun(x,y,val0)

call prin2("evaluation error = ",val0-val)



x = 0.412112d0
y = 0.53232d0
call tensor_eval(n,a,b,c,d,coefs,x,y,val)
call testfun(x,y,val0)
call prin2("evaluation error = ",val0-val)


a = 1
b = 2

c = 1
d = 2

nintsab = 12
nintscd = 12
allocate(ab(2,nintsab))
allocate(cd(2,nintscd))
allocate(coefs2(ncoefs,nintsab,nintscd))

do int=1,nintsab
ab(1,int) = a + (b-a)/(nintsab+0.0d0) * (int-1)
ab(2,int) = a + (b-a)/(nintsab+0.0d0) * (int)
end do

do int=1,nintscd
cd(1,int) = c + (d-c)/(nintscd+0.0d0) * (int-1)
cd(2,int) = c + (d-c)/(nintscd+0.0d0) * (int)
end do

call prin2("ab = ",ab)
call prin2("cd = ",cd)

do intab=1,nintsab
a0 = ab(1,intab)
b0 = ab(2,intab)
do intcd=1,nintscd

c0 = cd(1,intcd)
d0 = cd(2,intcd)

do i=1,n
x = xs(i) * (b0-a0)/2 + (b0+a0)/2
do j=1,n

y = xs(j) * (d0-c0)/2 + (d0+c0)/2

val = cos(24*x**2+y)
vals(i,j) = val

end do
end do

call tensor_coefs(n,u,vals,coefs2(:,intab,intcd))

end do
end do


x = 1.1d0
y = 1.4d0

call tensorpw_eval(n,nintsab,ab,nintscd,cd,coefs2,x,y,val)

val0   = cos(24*x**2+y)
errabs = val -val0
call prin2("evaluation error = ",errabs)

call tensorpw_evalder(n,nintsab,ab,nintscd,cd,coefs2,x,y,val,derx,dery)
val0   = cos(24*x**2+y)
derx0  = -sin(24*x**2+y)*24*2*x
dery0  = -sin(24*x**2+y)

errabs = val -val0
call prin2("evaluation error = ",errabs)

errabs = derx-derx0
call prin2("evaluation error = ",errabs)

errabs = dery-dery0
call prin2("evaluation error = ",errabs)

end program
