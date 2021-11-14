!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for constructing and evaluating expansions of the angular
!  prolate spheroidal wave functions of the first kind of integer orders and 
!  characteristic values in associated Legendre polynomials.
!
!  These prolate functions are the solutions of the spheroidal wave equation
!
!    (1-x^2) * y''(x) - 2x y'(x) + ( \chi - m^2 /(1-x^2) - \gamma^2 x^2 ) y(x) = 0          (1)
!
!  on the interval (-1,1) which satisfy the boundary conditions
!
!    \lim   y'(x) (1-x^2) = 0.                                                              (2)
!    x \to \pm 1
!
!  The equation (1) together with (2) is a singular self-adjoint Sturm-Liouville
!  problem.  Consequently, for fixed values of \gamma and m, there exist a countable
!  number of values of chi for which solutions of (1) satisfying (2) exist.  These 
!  values of chi are known as Sturm-Liouville eigenvalues, and we denote them by
!
!    \chi_m^m(\gamma^2) < \chi_{m_1}^m (\gamma^2) < \chi_{m+2}^m (\gamma^2) < ...
!
!  For each such eigenvalue, there is a corresponding one-dimensional eigenspace.
!  Here, we compute the eigenfunction
!
!      m
!    ps   (x; gamma^2)                                                                    (3)
!      n 
!
!  whose L^2(-1,1) norm is 1 and such that either (3) has the same sign as
!  the associated Legendre function
!
!      m
!    P  (x)                                                                               (4)
!      n
!
!  at 0 or (if the associated Legendre function is 0 at 0), the derivative of
!  (3) has the same sign as that of (4) at 0.
!
!  The function (3) admits an expansion of the form
!
!      m                 \infty     ~ m
!    ps (x,\gamma^2) =   \sum    a  P           (x),                                      (5)
!      n                 j = N    j   n + 2*j
!
!  where N = -Floor((n-m)/2) and
!
!    ~ m
!    P  (x)  
!      n 
!
!  denotes the L^2(-1,1) normalized associated Legendre polynomial of degree
!  n and order m.    This code operates by constructing a truncated version
!  (5) via the Osipov-Rokhlin method.
!
!  The following subroutines should be regarded as publicly callable:
!
!    spheroidal_integer - compute the the Sturm-Liouville eigenvalue chi_n^m(\gamma^2) 
!      and the coefficients in the associated Legendre expansions (5) of the function
!      ps_n^m(\gamma^2) using the Rokhlin-Osipov method
!
!    spheroidal_integer_eval - evaluate an expansion of the form (5) at a specified
!      point in (-1,1); the ALFs are calculated using the three-term recurrence relation
!      and some care is taken to mitigate problems associated with numerical underflow
!
!    spheroidal_alfs - evaluate the L^2 normalized associated Legendre polynomials
!
!        ~ m           ~ m
!        P  (x) , ..., P   (x)
!          m             m+n
!
!      at a specified point on the interval (-1,1)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module     spheroidal
use utils
contains

subroutine spheroidal_integer(gamma,n,m,chi,n1,n2,coefsps)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out)  :: coefsps(:)
!
!  Calculate the Sturm-Liouville eigenvalue chi_n^m(\gamma^2) and construct an
!  expansion of the form 
!
!      m                   n2    ~ m
!    ps (x,\gamma^2) ≈   \sum    a  P           (x)                                         (6)
!      n                 j = n1    j   n + 2*j
!
!  of the corresponding eigenfunction, which is the L^2(-1,1) normalized
!  angular spheroidal wave function of the first kind of bandlimit gamma,
!  order m and characteristic exponent n.
!                                                                                         
!  Input parameters:
!    gamma - the "bandlimit" of the PSWF function to construct, which must be
!      a nonnegative real number
!    n - the characteristic exponent, which must be a nonnegative integer
!    m - the order of the PSWF to construct, which must be a nonnegative integer
!
!  Output parameters:
!    chi - the Sturm-Liouville eigenvalue
!    (n1, n2)  - the lower and upper indices in the sum (6)
!    coefsps - an array containing the coefficients in the expansion (6)
!
double precision, allocatable :: as(:),bs(:),coefsps0(:)

pi     = acos(-1.0d0)
eps    = epsilon(0.0d0)

!  Estimate the number of terms in the expansion
!mm = 10 + (n-m) + floor(sqrt(gamma*(n-m)))
mm = 10 + 2/pi*(n-m) + floor(sqrt(gamma*(n-m)))

!mm = 20

1000 continue
!  Determine the limits of the index in the expansion
nn1 = max(-floor((n-m)/2.0d0),-mm) 
nn2 = mm


!  Construct the symmetric tridiagonal matrix; the diagonal elements
!  are stored in the array as and the off-diagonal entrys in the array
!  bs
allocate(coefsps0(nn1:nn2),as(nn1:nn2),bs(nn1+1:nn2))

!  Form the appropriate symmetric tridiagonal matrix 

dn = n
dm = m

do k=nn1,nn2
dk    = k
as(k) = (0.2d1*dk+dn)*(1+0.2d1*dk+dn)+gamma**2-(0.2d1*(-1+dm**2+  &
        (2*dk+dn)*(1+2*dk+dn))*gamma**2)/((-1+4*dk+2*dn)*(3+4*dk+  &
        2*dn))
end do

do k=nn1+1,nn2
dk    = k
bs(k) = (sqrt(((-0.1d1+0.2d1*dk-0.1d1*dm+dn)*(0.2d1*dk-0.1d1*dm+  &
        dn)*(-0.1d1+0.2d1*dk+dm+dn)*(0.2d1*dk+dm+dn))/((-0.3d1+  &
        0.4d1*dk+0.2d1*dn)*(0.1d1+0.4d1*dk+0.2d1*dn)))*gamma**2)/(-1  &
        +0.4d1*dk+0.2d1*dn)
end do


! Compute the specified eigenvalue using the Rokhlin-Osipov technique
nn  = nn2-nn1+1
idx = -nn1+1

call spheroidal_trieigen(nn,as,bs,idx,chi,coefsps0)

! Truncate the coefficient expansion so as to eliminate coefficients of 
! extremely small magnitude


dd  = abs(coefsps0(0))
eps = epsilon(0.0d0)**2


do i=nn1,0
if ( abs(coefsps0(i)) .ge. eps*dd) then
mm1 = i
exit
endif
end do

do i=nn2,0,-1
if ( abs(coefsps0(i)) .gt. eps*dd) then
mm2 = i
exit
endif
end do


!  If the number of coefficients is insufficient, repeat the process


if ( mm2 .eq. nn2) then
mm = mm*2
deallocate(coefsps0,as,bs)
goto 1000
endif

n1 = mm1
n2 = mm2

allocate(coefsps(n1:n2))
coefsps = coefsps0(n1:n2)

! this fixed the sign
dmax = 0
idx  = 0
do i=n1,n2
if (abs(coefsps(i)) .gt. dmax) then
dmax = abs(coefsps(i))
idx  = i
endif
end do

coefsps = coefsps * sign(1.0d0,coefsps(idx))

end subroutine


subroutine spheroidal_integer_eval(n,m,n1,n2,coefsps,x,valps)
implicit double precision (a-h,o-z)
double precision :: coefsps(n1:n2)
!
!  Evaluate the expansion (6) at a specified point x in order to
!  calculate the value of the angular prolate spheroidal wave function
!  ps_n^m(x;gamma^2) there.
!
!  Input parameters:
!    n - the characteristic exponent of the PSWF
!    m - the degree of the PSWF
!    n1 - the lower index in the sum in (4) as returned by spheroidal_integer
!    n2 - the upper index in the sum in (4) as returned by spheroidal_integer
!    coefsps - the coefficients in the expans (4) as returned by spheroidal_integer
!    x - the point in (-1,1) at which to evaluate the prolate function
!     ps_n(x,c) 
!
!  Output parameters:
!    valps - the value of the angular PSWF at the specified point x
!

double precision, allocatable :: alfs(:)


! Determine how many ALFs we need to evaluate
nn = n+2*n2
nn = nn-m+1


! Evaluate 'em
allocate(alfs(nn))
call spheroidal_alfs(m,nn,x,alfs)

! Sum up the expansion
valps = 0
idx   = n+2*n1-m+1

do i=n1,n2
valps = valps + alfs(idx)*coefsps(i)
idx   = idx + 2
end do


end subroutine


subroutine spheroidal_alfs(m,n,x,vals)
implicit double precision (a-h,o-z)
double precision :: vals(n)
! 
!  Evaluate the L^2 normalized associated Legendre functions
!  
!    ~ m       ~ m             ~ m  
!    P  (x) ,  P    (x), ...,  P   (x)
!      m         m+1             m+n-1
!
!  at a specified point x using the well-known three-term recurrence
!  relations.  Any values which are smaller than the smallest
!  representable double precision number are truncated to 0.
!  
!  Input parameters:
!    m - specifies the order of the ALFs to calculate
!    n - specifies the number of ALFs to calculate
!    x - the point in the interval (-1,1) at which to calculate them
!
!  Output parameters:
!    vals - an array whose jth entry contains the value of P_{m+j-1}^m(x)
!
double precision, allocatable :: z1(:),z2(:)
data pi / 3.14159265358979323846264338327950288d0 /

allocate(z1(n),z2(n))

! Evaluate the first two normalized ALFs
call spheroidal_pmm(m,x,z1(1),z2(1),z1(2),z2(2))

vals(1) = z1(1) * 10.0d0**(z2(1))
vals(2) = z1(2) * 10.0d0**(z2(2))

! Use the recurrence relation to evaluate the rest ... conduct
! the operations using extended precision arithemtic to avoid
! underflow

dm = m
dn = dm

do i=3,n

bb = -0.1d1/Sqrt(((2-dm+dn)*(2+dm+dn))/(15+16*dn+4*dn**2))
cc = ((1+dm+dn)*sqrt(((0.1d1-0.1d1*dm+dn)*(0.2d1-0.1d1*dm+  &
     dn)*(0.5d1+0.2d1*dn))/((0.1d1+dm+dn)*(0.2d1+dm+dn)*(0.1d1+  &
     0.2d1*dn))))/(2-0.1d1*dm+dn)
     
w2   = z2(i-1)
w1   = -bb*x*z1(i-1) -cc*z1(i-2)

do while ( abs(w1) .gt. 1.0d200)
w1=w1/1.0d200
w2=w2+200

z1(i-1) = z1(i-1)/1.0d200
z2(i-1) = z2(i-1)+200
end do

z1(i) = w1
z2(i) = w2

vals(i) = w1*10.0d0**(w2)
dn      = dn+1

end do

end subroutine


subroutine spheroidal_pmm(m,x,z1,z2,w1,w2)
implicit double precision (a-h,o-z)
! 
!  Evaluate the L^2 normalized associated Legendre functions
!  
!    ~ m                 ~ m
!    P  (x)    and       P   (x)
!      m                   m+1
!
!  at a specified point x.  Because of the potential for underflow,
!  the values of these functions are returned in the form
!
!     z1 * 10^z2.
!  
!  Input parameters:
!    m - specifies the order of the ALFs to calculate
!    x - the point in the interval (-1,1) at which to calculate them
!
!  Output parameters:
!    (z1,z2) - the value of the first ALF is z1 * exp(z2)
!    (w1,w2) - the value of the second ALF is z1 * exp(z2)
!

data pi / 3.14159265358979323846264338327950288d0 /

! ! Handle the special case m = 0
! if (m == 0) then
! val1 = 1.0d0/sqrt(2.0d0)
! val2 = sqrt(3.0d0/2.0d0) * x
! return
! endif

! Compute the value of dd = Gamma(m+3/2) / Gamma(m+1) using an asymptotic
! expansion when m is large and via a direct calculation when m is small.

dm = m
if (dm .gt. 100) then
dd = 0.1d1+0.220407270229753263690239406535908984d-1/dm**16+  &
     0.598095648752670211445447989717649762d-1/dm**15-  &
     0.476979015394840555330802089883945882d-2/dm**14-  &
     0.13033081777562571801354351919144392d-1/dm**13+  &
     0.142075389042872757272562012076377869d-2/dm**12+  &
     0.392863565457446384243667125701904297d-2/dm**11-  &
     0.61860934147262014448642730712890625d-3/dm**10-  &
     0.17494493513368070125579833984375d-2/dm**9+  &
     0.4302351735532283782958984375d-3/dm**8+  &
     0.12755692005157470703125d-2/dm**7-  &
     0.5538463592529296875d-3/dm**6-0.1842498779296875d-2/dm**5+  &
     0.1800537109375d-2/dm**4+0.87890625d-2/dm**3-  &
     0.546875d-1/dm**2+0.375d0/dm
dd = dd * sqrt(dm)
else
dd = gamma(dm+1.5d0)/ gamma(dm+1.0d0)
endif

! The normalization constant is z1 = sqrt(dd/sqrt(pi))
! dd = sqrt(dd/sqrt(pi))
! z1 = dd * (1-x**2)**(dm/2)
! z2 = 0
z1   = sqrt(dd/sqrt(pi))
z2   = dm/2*log10(1-x**2)

if (mod(m,2) == 1) z1=-z1

do while ( abs(z1) .gt. 1.0d6)
z1=z1/1.0d6
z2=z2+6
end do


! Now compute the second ALF, which is sqrt(3+2*dm)*x time the first one
w1 = z1 * sqrt(3+2*dm)*x
w2 = z2

end subroutine


subroutine spheroidal_trieigen(n,as,bs,j,chi,y)
implicit double precision (a-h,o-z)
double precision             :: as(n),bs(n),y(n)
!
!  Calculate a single eigenvalue and a corresponding eigenvector of a
!  symmetric tridiagonal matrix via the Rokhlin-Osipov approach.
!
!  Input parameters:
!    n - the dimensions of the matrix 
!    as - the n diagonal entries of the matrix
!    bs - the n-1 subdiagonal entries of the matrix
!    j - the index of the eigenvalue to compute
!
!  Output parameters:
!    chi - the eigenvalue
!    y - the eigenvector
!

double precision, allocatable   :: as0(:),bs0(:),x(:)
!
!  Use sturm bisection to estimate  the eigenvalue
!

eps = epsilon(0.0d0)*10
call spheroidals_sturm_bisection(eps,n,as,bs,j,chi)
!
!  Apply the inverse power method, starting with a random right hand side
!
allocate(as0(n),bs0(n),x(n))

call random_seed()
do i=1,n
call random_number(y(i))
end do
dd = 0
do i=1,n
dd = dd + y(i)**2
end do 
dd = sqrt(dd)
y  = y / dd

chi0 = chi
do iter  = 1,6

! solve the system ( A - chi * I ) x = y 
as0 = as - chi0
bs0 = bs

! eliminate the subdiagonal entries
do i=2,n
dd     = bs0(i-1) / as0(i-1)
as0(i) = as0(i) - dd * bs0(i-1)
y(i)   = y(i)   - dd * y(i-1)
end do

x(n) = y(n) / as0(n)

do i=n-1,1,-1
x(i) = y(i) - bs0(i) * x(i+1)
x(i) = x(i) / as0(i)
end do

! normalize the eigenvector
dd   = 0
imax = 0
dmax = 0

do i=1,n
dd = dd + x(i)**2

if (abs(x(i)) .gt. dmax) then
dmax = abs(x(i))
imax = i
endif

end do
dd = sqrt(dd)
y   = x / dd

end do


! ! update the eigenvalue
! dsum1=0

! x(1) =  as(1)*y(1)  + bs(1)*y(2)

! do i=2,n-1
! x(i) = bs(i-1)*y(i-1) + as(i) * y(i) + bs(i) * y(i+1)
! end do

! x(n) = bs(n-1)*y(n-1)+as(n)*y(n)

! call prin2("y = ",y)
! call prin2("x = ",x)

! chi0 = 0
! do i=1,n
! chi0 = chi0 + y(i) * x(i)
! end do

!chi = chi0
end subroutine




subroutine spheroidals_sturm_bisection(eps,n,as,bs,j,dlambda)
implicit double precision (a-h,o-z)

integer, intent(in)            :: n,j
double precision, intent(in)   :: eps,as(n),bs(n)
double precision, intent(out)  :: dlambda
!
!  Use the classic Sturm bisection algorithm in order to approximate the
!  jth eigenvalue of a symmetric tridiagonal matrix.
!
!  Input parameters:
!
!    eps - the desired precision for the calculations
!    n - the dimensions of the input matrix
!    as - an array of length n speciying the diagonal entries of the input
!      matrix
!    bs - an array of length n-1 specifying the off-diagonal entries of the
!      input matrix
!    j - the index (between 1 and n) of the eigenvalue to compute
!
!
!  Output parameters:
!    dlambda - an estimate of the jth eigenvalue of the input matrix
!

double precision, allocatable :: pols(:),pols0(:)


allocate(pols(0:n))

!
!  Bracket the eigenvalue using the Gershgorin theorem
!

x1 = as(1) - abs(bs(1))
x2 = as(1) + abs(bs(1))

do i=2,n-1
x1  = min(x1, as(i) - abs(bs(i-1)) - abs(bs(i)))
x2  = max(x2, as(i) + abs(bs(i-1)) + abs(bs(i)))
end do

x1  = min(x1, as(i) - abs(bs(i-1)))
x2  = max(x2, as(i) + abs(bs(i-1)))


!
!  Perform repeated bisection
!

do while ( abs(x2-x1) .gt. eps * (abs(x1)+abs(x2)))
x       = (x1+x2)/2

if (abs(x1) + abs(x2) .lt. 1d-200) exit

!
!  Evaluate the characteristic polynomials of the principal submatrices
!  of the input matrix via the standard recurrence relation at the point x.
!

pols(0) = 1.0d0
pols(1) = as(1)-x

do l=2,n

pols(l)  = (as(l)-x)*pols(l-1) - bs(l-1)**2 * pols(l-2)
dd       = pols(l)**2+pols(l-1)**2

if (dd .lt. 1.0d-13) then
pols(l)   = pols(l)   * 1.0d+24
pols(l-1) = pols(l-1) * 1.0d+24
endif

if (dd .gt. 1.0d13) then
pols(l)   = pols(l)   * 1.0d-24
pols(l-1) = pols(l-1) * 1.0d-24
endif

end do

!
!  Count the number of sign changes
!

nn = 0
do l=1,n 
if (pols(l)*pols(l-1) .lt. 0) nn = nn + 1
if (pols(l) .eq. 0)           nn = nn - 1
end do

!
!  There are nn eigenvalues less than x
!

if (j .le. nn) then
x2 = x
else
x1 = x
endif

end do

dlambda = x

end subroutine


end module 
