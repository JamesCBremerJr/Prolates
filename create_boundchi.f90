module boundexp_funs

use utils
use adapquad
use chebyshev
use iso_c_binding

type fun_data
integer                                 :: nints, k 
double precision, allocatable           :: ab(:,:), xs(:), vals(:,:)
end type

contains

subroutine fun1(n,ts,vals,userptr)
implicit double precision (a-h,o-z)
double precision, pointer      :: dk
type(c_ptr)                    :: userptr
double precision               :: ts(n),vals(n)

call c_f_pointer(userptr,dk)

do i=1,n
vals(i) = sqrt( (1 - dk**2 * ts(i)**2 ) / ( 1- ts(i)**2) )
end do

end subroutine


subroutine fun2(t,val,userptr)
implicit double precision (a-h,o-z)
double precision, target       :: dk
type(c_ptr)                    :: userptr
double precision, intent(in)   :: t 
double precision, intent(out)  :: val

dk      = t 
eps     = 1.0d-15
dd      = 1.0d-30

a       = 0.0d0
b       = min(1.0d0-dd,1.0d0/dk)
userptr = c_loc(dk)

call adapint(ier,eps,a,b,fun1,userptr,val,ntotal)
! val    = dk/val

end subroutine

subroutine fun3(t,val,userptr)
implicit double precision (a-h,o-z)
double precision, target       :: dk
type(c_ptr)                    :: userptr
double precision, intent(in)   :: t 
double precision, intent(out)  :: val

dk      = t 
eps     = 1.0d-15
dd      = 1.0d-30

a       = 0.0d0
b       = min(1.0d0-dd,1.0d0/dk)
userptr = c_loc(dk)

call adapint(ier,eps,a,b,fun1,userptr,val,ntotal)
val    = dk/val

end subroutine



subroutine funeval(t,val,userptr)
implicit double precision (a-h,o-z)
double precision, pointer      :: dk
type(c_ptr)                    :: userptr
type(fun_data), pointer        :: fundata
double precision, intent(in)   :: t 
double precision, intent(out)  :: val
call c_f_pointer(userptr,fundata)
call chebpw_eval(fundata%nints,fundata%ab,fundata%k,fundata%xs,fundata%vals,t,val)
end subroutine


subroutine inverse(nints,ab,k,chebdata,vals,ders,nintsinv,abinv,valsinv)
implicit double precision (a-h,o-z)
type(chebexps_data)                        :: chebdata
integer, intent(in)                        :: nints,k
double precision, intent(in)               :: vals(k,nints),ders(k,nints),ab(2,nints)
double precision, intent(out), allocatable :: abinv(:,:),valsinv(:,:)
integer, intent(out)                       :: nintsinv



nintsinv = nints
allocate(valsinv(k,nints),abinv(2,nints))

eps0     = 1.0d-15
maxiters = 20
nextra   = 2

do int=1,nints
a = ab(1,int)
b = ab(2,int)

call chebpw_eval(nints,ab,k,chebdata%xs,vals,a,a0)
call chebpw_eval(nints,ab,k,chebdata%xs,vals,b,b0)

abinv(1,int) = a0
abinv(2,int) = b0

end do

call prin2("abinv = ",abinv)

!
!  Use Newton's method to evaluate the inverse at each of the Chebyhev nodes 
!  on the grid defined by abinv
!

do int = nints,1,-1

a0 = abinv(1,int)
b0 = abinv(2,int)
a  = ab(1,int)
b  = ab(2,int)

do i = k,1,-1

x  = (b0-a0)/2*chebdata%xs(i) + (b0+a0)/2
t  = (b-a)/2*chebdata%xs(i) + (b+a)/2

do iter=1,maxiters+1

if (iter .eq. maxiters+1) then
call prina("in kummer_phase_invert: maximum number of Newton iterations exceeded")
stop
endif

call chebpw_eval(nints,ab,k,chebdata%xs,vals,t,val)
call chebpw_eval(nints,ab,k,chebdata%xs,ders,t,der)

delta = (val-x)/der
if (abs(delta) .lt. eps0*(1+abs(t))) exit
t     = t - delta

end do


do iextra=1,nextra

call chebpw_eval(nints,ab,k,chebdata%xs,vals,t,val)
call chebpw_eval(nints,ab,k,chebdata%xs,ders,t,der)

delta = (val-x)/der
t     = t - delta

end do

valsinv(i,int) = t 

end do
end do

end subroutine


end module


program boundexp
use utils
use adapquad
use chebyshev
use utils
use boundexp_funs
implicit double precision (a-h,o-z)

type(fun_data), target        :: psi
type(fun_data), target        :: psider
type(fun_data), target        :: phi
type(c_ptr)                   :: userptr
type(chebexps_data)           :: chebdata
double precision, allocatable :: ab(:,:),ab2(:,:),vals(:,:),ders(:,:),abin(:,:)
double precision, allocatable :: abinv(:,:),valsinv(:,:)

k = 30
call chebexps(k,chebdata)

eps = 1.0d-15
a   = 0.0d0
b   = 1.0d5

nintsin = 2
allocate(abin(2,nintsin))
abin(1,1) = a
abin(2,1) = 1.0d0

abin(1,2) = 1.0d0
abin(2,2) = b


!
!  Discretize the function  
!
!    Psi(x) =   x / E(x),
!
!  where
!
!                  min(1,1/x)  (   1 - x^2 t^2  ) ^(1/2) 
!    E(x) =   \int             (  ------------- )          dt
!                  0           (    1 - t^2     )
!  
!  and over the interval [a,b].
!

!call chebadap(ier,eps,a,b,fun2,k,chebdata,nints,ab,userptr)
call chebadap2(ier,eps,nintsin,abin,fun2,k,chebdata,nints,ab,userptr)

call prin2("ab = " ,ab)

!
!  Store the data necessary to evaluate Psi and Psi'
! 

allocate(psi%xs(k), psi%ab(2,nints), psi%vals(k,nints) )

psi%nints = nints
psi%k     = k
psi%xs    = chebdata%xs
psi%ab    = ab

do int=1,nints
a0 = ab(1,int)
b0 = ab(2,int)
do i=1,k
t = (b0-a0)/2 * chebdata%xs(i) + (b0+a0)/2
call fun3(t,val,userptr)
psi%vals(i,int) = val
end do
end do

allocate(psider%xs(k), psider%ab(2,nints), psider%vals(k,nints) )
psider%nints = nints
psider%k     = k
psider%xs    = chebdata%xs
psider%ab    = ab

do int=1,nints
a0 = ab(1,int)
b0 = ab(2,int)
psider%vals(:,int) = matmul(chebdata%adiff, psi%vals(:,int))
psider%vals(:,int) = 2/(b0-a0) * psider%vals(:,int)
end do

userptr   = c_loc(psi)
call funeval(0.1d0,val,userptr)

userptr   = c_loc(psider)
call funeval(0.1d0,der,userptr)

! print *,""
! print *,val-0.0638218322355659653598858419352843566d0
! print *,der-0.641429498958836457413871743830184255d0

!
!  Invert the function Psi(x) using Newton's method to form the function Phi(x)
!

call inverse(nints,ab,k,chebdata,psi%vals,psider%vals,nintsinv,abinv,valsinv)


allocate(phi%xs(k), phi%ab(2,nints), phi%vals(k,nints) )
phi%nints = nintsinv
phi%k     = k
phi%xs    = chebdata%xs
phi%ab    = abinv
phi%vals  = valsinv

userptr   = c_loc(phi)
call funeval(0.1d0,val,userptr)

val0 = 0.15611809376344279514d0
print *,val
print *,val0
print *,val-val0


print *,""
print *,""

call funeval(0.1d0,val,userptr)
val0 = 2.82544588302492423566949143279634023854624752267588456758314d0
print *,val
!
!  Output the fortran code for evaluating Phi(x)
!


iw = 20
open(iw,FILE="boundchi.f90")


0900 format(A)
1000 format(A,I8.8,A)
1100 format(A,I8.8,A,I8.8,A)
1200 format(A,D24.15,A)

write(iw,"(A)") "subroutine boundchi(c,n,chi1,chi2)"
write(iw,"(A)") "implicit double precision (a-h,o-z)"
write(iw,"(A)") "data pi / 3.14159265358979323846264338327950288d0 /"
write(iw,"(A)") "x1 = 2*c / ( pi * n )"
write(iw,"(A)") "x2 = 2*c / ( pi * (n+1) )"
write(iw,"(A)") "call evalphi(x1,val1)"
write(iw,"(A)") "call evalphi(x2,val2)"
write(iw,"(A)") "chi1 = ( c / val1 )**2"
write(iw,"(A)") "chi2 = ( c / val2 )**2"
write(iw,"(A)") "end subroutine"
write(iw,"(A)") ""

write(iw,"(A)") "subroutine evalphi(x,val)"
write(iw,"(A)") "implicit double precision (a-h,o-z)"
write(iw,"(A)") "integer             :: nints, k"
write(iw,1000)  "double precision    :: xscheb(",phi%k,")"
write(iw,1100)  "double precision    :: ab(",2,",",phi%nints,")"
write(iw,1100)  "double precision    :: vals(",phi%k,",",phi%nints,")"

write(iw,1000)  "data nints      / ",phi%nints," /"
write(iw,1000)  "data k          / ",phi%k," /"

write(iw,1200)  "data xscheb     / ",phi%xs(1),", &"
do i=2,phi%k-1
write(iw,1200)  "                  ",phi%xs(i),", &"
end do
write(iw,1200)  "                  ",phi%xs(phi%k),"  /"


write(iw,1200)  "data ab         / ",phi%ab(1,1),", &"
write(iw,1200)  "                  ",phi%ab(2,1),", &"

do int=2,phi%nints-1
do i=1,2
write(iw,1200)  "                  ",phi%ab(i,int),", &"
end do
end do

write(iw,1200)  "                  ",phi%ab(1,phi%nints),", &"
write(iw,1200)  "                  ",phi%ab(2,phi%nints),"  /"



write(iw,1200)  "data vals       / ",phi%vals(1,1),", &"
do i=2,k
write(iw,1200)  "                  ",phi%vals(i,1),", &"
end do

do int=2,phi%nints-1
do i=1,k
write(iw,1200)  "                  ",phi%vals(i,int),", &"
end do
end do

do i=1,phi%k-1
write(iw,1200)  "                  ",phi%vals(i,phi%nints),", &"
end do

write(iw,1200)  "                  ",phi%vals(k,phi%nints),"  /"

write(iw,"(A)") ""
write(iw,"(A)") "do int=1,nints-1"
write(iw,"(A)") "if(x .lt. ab(2,int)) exit"
write(iw,"(A)") "end do"
write(iw,"(A)") "a = ab(1,int)"
write(iw,"(A)") "b = ab(2,int)"
write(iw,"(A)") "xx   = (2*x - (b+a) ) /(b-a)"
write(iw,"(A)") "sum1=0"
write(iw,"(A)") "sum2=0"
write(iw,"(A)") "dd1 = 1.0d0"
write(iw,"(A)") "do i=1,k"
write(iw,"(A)") " dd=1.0d0"
write(iw,"(A)")" if (i .eq. 1 .OR. i .eq. k) dd = 0.5d0"
write(iw,"(A)") "diff = xx-xscheb(i)"
write(iw,"(A)") "if(abs(diff) .le. eps0 ) then"
write(iw,"(A)") "val = vals(i,int)"
write(iw,"(A)") "return"
write(iw,"(A)") " endif"
write(iw,"(A)") "dd   = (dd1*dd)/diff"
write(iw,"(A)") "dd1  = - dd1"
write(iw,"(A)") "sum1 = sum1+dd*vals(i,int)"
write(iw,"(A)") "sum2 = sum2+dd"
write(iw,"(A)") "dd   = - dd"
write(iw,"(A)") "end do"
write(iw,"(A)") "val = sum1/sum2"
write(iw,"(A)") "end subroutine"
write(iw,"(A)") ""


close(iw)

end program
