!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for evaluating the angular prolate spheroidal wave functions
!  of order 0 and integer characteristic exponents -- that is, the functions
!
!  It does contain code for computing the Sturm-Liouville eigenvalues of
!  prolate functions of noninteger characterisitc exponents, but not mechanism for
!  evaluating such functions is provided.
!
!  This code uses the Xiao-Rokhlin approach and is meant to be used in the regime
!  in which the bandlimit and characteristic exponent are relatively small.
!
!  The functions ps_n(x,c) and qs_n(x,c) are normalized as in Meixner and Schafke's
!  "Mathieucsche Funktionen und Spharoidfunktionen" so that, for instance,
!
!         1
!    \int    ps_n(x,c)^2 dx = 2 / (2n+1).
!        -1
!
!  The user is cautioned that this definition of qs_n(x,c) is problematic.  It
!  is of extremely small magnitude when c is large, as is the Wronskian of the
!  pair ps_n(x,c) and qs_n(x,c).
!
!  The following subroutines are publicly callable:
!
!  prolates_integer - construct the coefficients in the Legendre expansions of the
!    spheroidal wave function of the first and second kinds of a nonnegative integer
!    characteristic exponent n
!
!  prolates_integer_eigenvalue - compute the Sturm-Liouville eigenvalue for the
!    prolate spheroidal wave function ps_n(x)
!
!  prolates_integer_eval - evaluate the angular prolate functions whose Legendre
!    expansions are constructed by prolates_integer at a specified point
!
!  prolates_integer_eval0 - evaluate only the function ps_n(x,c) using the
!    expansion constructed by prolates_integer
!
!  prolates_integer_evalder - evaluate the angular prolate functions whose
!    Legendre expansions are constructed by prolates_integer and their
!    derivatives at a specified point
!
!  prolates_integer_bound - use an asymptotic estimate of Bonami and Karoui to
!    bound the Sturm-Liouville eigenvalue of one of the prolate sperhoidal
!    wave functions of integer characteristic value
!
!  prolates_integer_coefficients - compute coefficients c1 and c2 such that
!     the square modulus of the function
!
!      c1 ps_n(x) + i c2 qs_n(x)
!
!    is completely monotone
!
!    IMPORTANT WARNING: this calculation is numerically unstable for large
!    c and certain values of n -- use it at your own risk.  There are
!    better ways of computing the nonoscillatory phase but this routine
!    is useful for testing in certain circumstances
!
!  prolates_noninteger_all - compute the Legendre expansions of the angular
!    spheroidal wave functions of the first and second kinds of characteristic
!    exponents 
!
!        dnu0 + k     for    k=0,1,2,3,...,n
!
!    where -1/2 < dnu0 < 1/2 and n is specified by the user
!
!  prolates_noninteger_all_eval - evaluate one of the angular prolate spheroidal
!    wave functions whose Legendre expansions are constructed by 
!    prolates_noninteger_all
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module     prolates

use utils

contains



subroutine prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out)  :: coefsps(:), coefsqs(:)
double complex                              :: clambda
!
!  Calculate approximate expansions
!
!                    M
!    ps_n(x,c) =   \sum  a_j P_{2j + mod(n,2) }(x)
!                   j=0
!                                                                                         (1)
!                    M                                    M
!    qs_n(x,c) =   \sum  b_j P_{2j + 1 - mod(n,2)}(x) + \sum  a_j Q_{2j+mod(n,2)} (x)
!                   j=0                                  j=0
!
!  of the angular prolate spheroidal wave functions of the first and second
!  kinds of a nonnegative integer characteristic exponent  This routine determines an 
!  appropriate value for M.
!
!  Input parameters:
!    c - the bandlimit of the prolate fnuctions to construct
!    n - the characteristic exponent, which must be a nonnegative integer
!
!  Output parameters:
!    chi - the Sturm-Liouville eigenvalue of ps_n
!    clambda - the integral eigenvalue of ps_n
!    ncoefs - the number of coefficients in each term of the expansions
!      in (1); i.e., the integer M
!    coefsps - an array containing the first set of coefficients {a_n} in (1)
!    coefsqs - an array containing the second set of coefficinets {b_n} in (1)
!

k = n/2

!
!  Determine the size of the matrices and allocate the output 
!  arrays
!

m      = max(ceiling(1.0d0*(n+sqrt(n*c))),50)
ncoefs = m

allocate(coefsps(m),coefsqs(m))
!
!  Call an auxillary routine depending on whether n is even or odd
!  to compute the coefficients
!

if (2*k .eq. n) then
call prolates_integer_even(c,n,m,k,chi,coefsps,coefsqs)
else
call prolates_integer_odd(c,n,m,k,chi,coefsps,coefsqs)
endif

!
!  Compute the coefficients and the integral eigenvalue
!

call prolates_integral_eigenvalue(c,n,chi,ncoefs,coefsps,coefsqs,clambda)

end subroutine


subroutine prolates_integer_even(c,n,m,k,chi,coefsps,coefsqs)
implicit double precision (a-h,o-z)
double precision :: coefsps(m), coefsqs(m)
!
!  An auxillary rotuine called by prolates_integer which computes the
!  coefficient expansions in the case of an even characteristic exponent.
!

double precision, allocatable         :: as(:), bs(:), cs(:), work(:), w(:), z(:,:)
integer, allocatable                  :: iwork(:)
logical                               :: tryrac


!
!  Form the symmetric tridiagonal matrix to compute the first set of
!  coefficients -- the diagonal is bs, the off-diagonal is as
!

allocate(as(m),bs(m),cs(m))

do j=1,m-1
as(j) = (2.0d0*c**2*j*(-1.0d0 + 2*j))/((-1.0d0 + 4*j)*Sqrt((-3.0d0 + 4*j)*(1.0d0 + 4*j)))
end do

do j=0,m-1
bs(j+1) = 2d0*j*(1d0 + 2*j) + (c**2*(-1d0 + 4*j + 8*j**2))/(-3d0 + 8*j + 16*j**2)
end do


!
!  Compute the eigenvalue
!

call prolate_trieigen(m,bs,as,k+1,chi,coefsps)

! il     = k+1
! iu     = k+1
! tryrac = .TRUE.
! liwork = 10*m
! lwork  = 20*m
! nzc    = 1


! allocate(iwork(liwork),work(lwork))
! allocate(w(m),z(m,1))


!  call dstemr('V','I',m,bs,as,vl,vu,il,iu,mm,w,z,m,nzc,isuppz,tryrac, &
!    work,lwork,iwork,liwork,info)


! chi = w(1)
! do j=0,m-1
! coefsps(j+1) = z(j+1,1) 
! end do

! deallocate(iwork,work,w,z)

!
!  Normalize the eigenvectors 
!

do j=0,m-1
coefsps(j+1) = coefsps(j+1) * sqrt( (4*j+1.0d0)/2.0d0)
end do

rn  = 0
do j=0,m-1
rn = rn + coefsps(j+1)**2 * (2*n+1.0d0)/(4*j+1.0d0)
end do
coefsps = coefsps / sqrt(rn)
coefsps = coefsps * sign(1.0d0,coefsps(1))


!
!  Compute the second set of coefficients by solving the appropriate linear system
!  of equations
!

allocate(w(m))


do j=1,m-1
as(j) = (2*c**2*j*(1.0d0 + 2*j))/(-1.0d0 + 16*j**2)
end do

do j=0,m-2
cs(j+1) = (2*c**2*(1.0d0 + j)*(3.0d0 + 2*j))/((5.0d0 + 4*j)*(7.0d0 + 4*j))
end do

do j=0,m-1
bs(j+1) = 2.0d0*(1.0d0 + j)*(1.0d0 + 2*j) + (c**2*(3.0d0 + 4*j*(3.0d0 + 2*j))) &
   /((1.0d0 + 4*j)*(5.0d0 + 4*j)) - chi
end do



w    = 0 
w(1) = -c**2 * coefsps(1)

!call dgtsv(m,1,as,bs,cs,w,m,info)
call prolates_trisolve(m,as,bs,cs,w)

do i=1,m
coefsqs(i) = w(i)
end do


deallocate(w)



end subroutine


subroutine prolates_integer_odd(c,n,m,k,chi,coefsps,coefsqs)
implicit double precision (a-h,o-z)
double precision :: coefsps(m), coefsqs(m)
!
!  An auxillary rotuine called by prolates_integer which computes the
!  coefficient expansions in the case of an odd characteristic exponent.
!

double precision, allocatable         :: as(:), bs(:), cs(:), work(:), w(:), z(:,:)
integer, allocatable                  :: iwork(:)
logical                               :: tryrac

!
!  Form the symmetric tridiagonal matrix to compute the first set of
!  coefficients -- the diagonal is bs, the off-diagonal is as
!


allocate(as(m),bs(m),cs(m))

do j=1,m-1
as(j) = (2*c**2*j*(1.0d0 + 2*j))/((1.0d0 + 4*j)*Sqrt((-1.0d0 + 4*j)*(3.0d0 + 4*j)))
end do

do j=0,m-1
bs(j+1) = 2.0d0 + 6*j + 4*j**2 + (c**2*(3.0d0 + 4*j*(3.0d0 + 2*j)))/((1.0d0 + 4*j)*(5.0d0 + 4*j))
end do


!
!  Compute the appropriate eigenvalue & eigenvector
!
call prolate_trieigen(m,bs,as,k+1,chi,coefsps)

! il     = k+1
! iu     = k+1
! tryrac = .TRUE.
! liwork = 10*m
! lwork  = 20*m
! nzc    = 1

! allocate(iwork(liwork),work(lwork))
! allocate(w(m),z(m,1))


! call dstemr('V','I',m,bs,as,vl,vu,il,iu,mm,w,z,m,nzc,isuppz,tryrac, &
!   work,lwork,iwork,liwork,info)


! !
! ! Copy out the relevant eigenvalue; scale and normalize the coefficients
! !


! chi = w(1)
! do j=0,m-1
! coefsps(j+1) = z(j+1,1) 
! end do

! deallocate(w,z,iwork,work)


!
!  Normalize the coefficients
!

do j=0,m-1
coefsps(j+1) = coefsps(j+1) * sqrt( (4*j+3.0d0)/2.0d0)
end do

rn  = 0
do j=0,m-1
rn = rn + coefsps(j+1)**2 * (4*k+3.0d0)/(4*j+3.0d0)
end do
coefsps = coefsps / sqrt(rn)
coefsps = coefsps * sign(1.0d0,coefsps(1))


!
!  Compute the second set of coefficients by solving the appropriate linear system
!  of equations
!


do j=1,m-1
as(j) = (2*c**2*j*(-1.0d0 + 2*j))/((-3.0d0 + 4*j)*(-1.0d0 + 4*j))
end do

do j=0,m-2
cs(j+1)=(2*c**2*(1.0d0 + j)*(1.0d0 + 2*j))/((3.0d0 + 4*j)*(5.0d0 + 4*j))
end do

do j=0,m-1
bs(j+1) = 2*j*(1.0d0 + 2*j) + (c**2*(-1.0d0 + 4*j + 8*j**2))/((-1.0d0 + 4*j)*(3.0d0 + 4*j)) - chi
end do

allocate(w(m))
w    = 0 
w(1) = -c**2/3 * coefsps(1)

!call dgtsv(m,1,as,bs,cs,w,m,info)
call prolates_trisolve(m,as,bs,cs,w)

do i=1,m
coefsqs(i) = w(i)
end do

end subroutine



subroutine prolates_integral_eigenvalue(c,n,chi,ncoefs,coefsps,coefsqs,clambda)
implicit double precision (a-h,o-z)

double precision                    :: c,chi
integer                             :: ncoefs
double precision                    :: coefsps(ncoefs),coefsqs(ncoefs)
double complex                      :: clambda

!
!  Return the eigenvalue of the restricted Fourier operator corresponding
!  to PS_n(x).  That is, return the value of lambda such that
!
!                              1
!    T_c [PS_n](x) = \lambda \int   exp(ixt) PS_n(t) dt.
!                             -1
!
!  Input parameters:
!    c - the bandlimit
!    n - the characteristic exponent of the 
!    chi - the Sturm-Liouville eigenvalue corresponding of PS_n
!    ncoefs - the number of coefficients in the Legendre expansion of PS_n
!    coefsp - the coefficients in the Legendre expansion of PS_n
!
!  Output parameters:
!    clambda - the eigenvalue of T_c corresponding to PS_n
!

double complex, parameter      :: ima = (0.0d0,1.0d0)

ifeven = 0
if (mod(n,2) .eq. 0) ifeven = 1

call prolates_integer_evalder(n,ncoefs,coefsps,coefsqs,0.0d0,valps,derps,valqs,derqs)

if (ifeven .eq. 1) then
rint    = 2*coefsps(1)
clambda = rint/valps
else
rint    = 2*coefsps(1)
clambda = rint/derps * c/3.0d0 * ima
endif

end subroutine


subroutine prolates_integer_eval(n,ncoefs,coefsps,coefsqs,x,valps,valqs)
implicit double precision (a-h,o-z)

double precision :: coefsps(:), coefsqs(:)
!
!  Evaluate the expansions
!
!                    M
!    ps_n(x,c) =   \sum  a_j P_{2j + mod(n,2) }(x)
!                   j=0
!                                                                                         (3)
!
!                    M                                    M
!    qs_n(x,c) =   \sum  b_j P_{2j + 1 - mod(n,2)}(x) + \sum  a_j Q_{2j+mod(n,2)} (x)
!                   j=0                                  j=0
!
!  given the sequences {a_n} and {b_n} of coefficients returned by prolates_integer.
!
!  Input parameters:
!    n - the characteristic exponent of the prolate functions being evaluated
!    ncoefs - the number of coefficients in each  sum in the expansions (3)
!    coefsps - an array containing the coefficients {a_j}
!    coefsqs - an array containing the coefficients {b_j}
!    x - the point in (-1,1) at which to evaluate the prolate functions
!     ps_n(x,c) and qs_n(x,c)
!
!  Output parameters:
!    valps - the value of ps_n(x,c) 
!    valqs - the value of qs_n(c,x)
!

!double precision :: pvals(0:2*ncoefs-1),qvals(0:2*ncoefs-1)

double precision, allocatable :: pvals(:),qvals(:)

allocate(pvals(0:2*ncoefs),qvals(0:2*ncoefs))


!
!  Evaluate P_j and Q_j for j=0,1,2,... via the three-term recurrence relations
!


pvals(0) = 1
pvals(1) = x

qvals(0) = 0.5d0 * log( (1+x) / (1-x))
qvals(1) = -1.0d0 + x*qvals(0)

do j=1,2*ncoefs-2
dd1=  (2*j+1.0d0)/ (j+1.0d0)*x
dd2 = (j+0.0d0)/ (j+1.0d0)

pvals(j+1) = dd1* pvals(j) - dd2*pvals(j-1) 
qvals(j+1) = dd1* qvals(j) - dd2*qvals(j-1)
end do



!
!  Now sum the two series
!

ii = mod(n,2)

valps = 0
valqs = 0

do i=0,ncoefs-1
valps = valps + coefsps(i+1) * pvals(2*i+ii)
valqs = valqs + coefsps(i+1) * qvals(2*i+ii)
end do

do i=0,ncoefs-1
valqs = valqs + coefsqs(i+1) * pvals(2*i+1-ii)
end do

end subroutine


subroutine prolates_integer_eval0(n,ncoefs,coefsps,x,valps)
implicit double precision (a-h,o-z)

double precision :: coefsps(:)
!
!  Evaluate the expansion
!
!                    M
!    ps_n(x,c) =   \sum  a_j P_{2j + mod(n,2) }(x)
!                   j=0
!
!  given the sequences {a_n} of coefficients returned by prolates_integer.
!
!  Input parameters:
!    n - the characteristic exponent of the prolate functions being evaluated
!    ncoefs - the number of coefficients in each  sum in the expansions (3)
!    coefsps - an array containing the coefficients {a_j}
!    x - the point in (-1,1) at which to evaluate the prolate functions
!     ps_n(x,c) and qs_n(x,c)
!
!  Output parameters:
!    valps - the value of ps_n(x,c) 
!

double precision, allocatable :: pvals(:)

allocate(pvals(0:2*ncoefs))


!
!  Evaluate P_j and Q_j for j=0,1,2,... via the three-term recurrence relations
!


pvals(0) = 1
pvals(1) = x


do j=1,2*ncoefs-2
dd1=  (2*j+1.0d0)/ (j+1.0d0)*x
dd2 = (j+0.0d0)/ (j+1.0d0)

pvals(j+1) = dd1* pvals(j) - dd2*pvals(j-1) 
end do


!
!  and sum the two series
!

ii = mod(n,2)

valps = 0

do i=0,ncoefs-1
valps = valps + coefsps(i+1) * pvals(2*i+ii)
end do


end subroutine


subroutine prolates_integer_evalder(n,ncoefs,coefsps,coefsqs,x,valps,derps,valqs,derqs)
implicit double precision (a-h,o-z)

double precision :: coefsps(:), coefsqs(:)
!
!  Evaluate the expansions
!
!                    M
!    ps_n(x,c) =   \sum  a_j P_{2j + mod(n,2) }(x)
!                   j=0
!                                                                                         (3)
!
!                    M                                    M
!    qs_n(x,c) =   \sum  b_j P_{2j + 1 - mod(n,2)}(x) + \sum  a_j Q_{2j+mod(n,2)} (x)
!                   j=0                                  j=0
!
!  and their derivatives given the sequences {a_n} and {b_n} of coefficients returned by 
!  prolates_integer.
!
!  Input parameters:
!    n - the characteristic exponent of the prolate functions being evaluated
!    ncoefs - the number of coefficients in each  sum in the expansions (3)
!    coefsps - an array containing the coefficients {a_j}
!    coefsqs - an array containing the coefficients {b_j}
!    x - the point in (-1,1) at which to evaluate the prolate functions
!     ps_n(x,c) and qs_n(x,c)
!
!  Output parameters:
!    valp - the value of ps_n(x,c) 
!    derp - the derivative of ps_n(x,c) 
!    valq - the value of qs_n(c,x)
!    derq - the derivative of qs_n(x,c)
!

! double precision :: pvals(0:2*ncoefs),qvals(0:2*ncoefs)
! double precision:: pders(0:2*ncoefs),qders(0:2*ncoefs)

double precision, allocatable :: pvals(:),qvals(:)
double precision, allocatable :: pders(:),qders(:)

allocate(pvals(0:2*ncoefs),qvals(0:2*ncoefs))
allocate(pders(0:2*ncoefs),qders(0:2*ncoefs))

!
!  Evaluate P_j and Q_j for j=0,1,2,..., 2*ncoefs via the three-term recurrence relations
!


 if (abs(1-x) .lt. 1.0d-14) then


pvals = 1
do i=1,2*ncoefs
pders(i) = (i-1)*i/2.0d0
end do


 else
pvals(0) = 1
pvals(1) = x

qvals(0) = 0.5d0 * log( (1+x) / (1-x))
qvals(1) = -1.0d0 + x*qvals(0)

do j=1,2*ncoefs-1
dd1=  (2*j+1.0d0)/ (j+1.0d0)*x
dd2 = (j+0.0d0)/ (j+1.0d0)

pvals(j+1) = dd1* pvals(j) - dd2*pvals(j-1) 
qvals(j+1) = dd1* qvals(j) - dd2*qvals(j-1)
end do


!
!  Evaluate the derivatives of P_j and Q_j for j=0,1,...,2ncoefs-1 via the
!  well-known formulas
!


dd1 = 1/(1-x**2)
dd2 = x*dd1

do j=0,2*ncoefs-1
pders(j) = (j+1.0d0) * (dd2 *pvals(j) - dd1 * pvals(j+1))
qders(j) = (j+1.0d0) * (dd2 *qvals(j) - dd1 * qvals(j+1))
end do
endif


!
!  Now sum the relevant series
!

ii = mod(n,2)

valps = 0
derps = 0
valqs = 0
derqs = 0

do i=0,ncoefs-1
valps = valps + coefsps(i+1) * pvals(2*i+ii)
derps = derps + coefsps(i+1) * pders(2*i+ii)
valqs = valqs + coefsps(i+1) * qvals(2*i+ii)
derqs = derqs + coefsps(i+1) * qders(2*i+ii)
end do

do i=0,ncoefs-1
valqs = valqs + coefsqs(i+1) * pvals(2*i+1-ii)
derqs = derqs + coefsqs(i+1) * pders(2*i+1-ii)
end do

end subroutine






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine prolates_integer_eignvalue(c,n,chi)
implicit double precision (a-h,o-z)
!
!  Compute the Sturm-Liouville eigenvalue chi for the prolate spheroidal
!  wave function ps_n(x,c).
!
!  Input parameters:
!    c - the bandlimit of the prolate fnuctions to construct
!    n - the characteristic exponent, which must be a nonnegative integer
!
!  Output parameters:
!    chi - the Sturm-Liouville eigenvalue of ps_n
!

k = n/2

!
!  Determine the size of the matrices and allocate the output 
!  arrays
!

m      = max(ceiling(1.0d0*(n+sqrt(n*c))),50)

!
!  Call an auxillary routine depending on whether n is even or odd
!  to compute the coefficients
!
if (2*k .eq. n) then
call prolates_integer_even0(c,n,m,k,chi)
else
call prolates_integer_odd0(c,n,m,k,chi)
endif


end subroutine


subroutine prolates_integer_even0(c,n,m,k,chi)
implicit double precision (a-h,o-z)
!
!  An auxillary rotuine called by prolates_integer which computes the
!  coefficient expansions in the case of an even characteristic exponent.
!

double precision, allocatable         :: as(:), bs(:), coefsps(:),work(:),w(:),z(:,:)
integer, allocatable                  :: iwork(:)
logical                               :: tryrac

!
!  Form the symmetric tridiagonal matrix to compute the first set of
!  coefficients -- the diagonal is bs, the off-diagonal is as
!
allocate(as(m),bs(m),coefsps(m))

do j=1,m-1
as(j) = (2.0d0*c**2*j*(-1.0d0 + 2*j))/((-1.0d0 + 4*j)*Sqrt((-3.0d0 + 4*j)*(1.0d0 + 4*j)))
end do

do j=0,m-1
bs(j+1) = 2d0*j*(1d0 + 2*j) + (c**2*(-1d0 + 4*j + 8*j**2))/(-3d0 + 8*j + 16*j**2)
end do

!
!  Compute the eigenvalue
!
!call prolate_trieigen(m,bs,as,k+1,chi,coefsps)


il     = k+1
iu     = k+1
tryrac = .TRUE.
liwork = 10*m
lwork  = 20*m
nzc    = 1


allocate(iwork(liwork),work(lwork))
allocate(w(m),z(m,1))


 call dstemr('V','I',m,bs,as,vl,vu,il,iu,mm,w,z,m,nzc,isuppz,tryrac, &
   work,lwork,iwork,liwork,info)




chi = w(1)
! do j=0,m-1
! coefsps(j+1) = z(j+1,1) 
! end do

deallocate(iwork,work,w,z)

end subroutine


subroutine prolates_integer_odd0(c,n,m,k,chi)
implicit double precision (a-h,o-z)
!
!  An auxillary rotuine called by prolates_integer which computes the
!  coefficient expansions in the case of an odd characteristic exponent.
!

double precision, allocatable         :: as(:), bs(:), coefsps(:),work(:),w(:),z(:,:)
integer, allocatable                  :: iwork(:)
logical                               :: tryrac

!
!  Form the symmetric tridiagonal matrix to compute the first set of
!  coefficients -- the diagonal is bs, the off-diagonal is as
!


allocate(as(m),bs(m),coefsps(m))

do j=1,m-1
as(j) = (2*c**2*j*(1.0d0 + 2*j))/((1.0d0 + 4*j)*Sqrt((-1.0d0 + 4*j)*(3.0d0 + 4*j)))
end do

do j=0,m-1
bs(j+1) = 2.0d0 + 6*j + 4*j**2 + (c**2*(3.0d0 + 4*j*(3.0d0 + 2*j)))/((1.0d0 + 4*j)*(5.0d0 + 4*j))
end do


!
!  Compute the appropriate eigenvalue & eigenvector
!

!call prolate_trieigen(m,bs,as,k+1,chi,coefsps)

il     = k+1
iu     = k+1
tryrac = .TRUE.
liwork = 10*m
lwork  = 20*m
nzc    = 1


allocate(iwork(liwork),work(lwork))
allocate(w(m),z(m,1))


 call dstemr('V','I',m,bs,as,vl,vu,il,iu,mm,w,z,m,nzc,isuppz,tryrac, &
   work,lwork,iwork,liwork,info)




chi = w(1)
! do j=0,m-1
! coefsps(j+1) = z(j+1,1) 
! end do

deallocate(iwork,work,w,z)


end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine prolates_noninteger_all(c,dnu0,n,m,chis,coefs)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out)  :: coefs(:,:),chis(:)
!
!  Compute the coefficients in the expansions
! 
!                          m
!    ps_{dnu0+k}(x,c) =  \sum    a_{j,k} P_{dnu0+k+2j} (x)
!                         j=-m
!
!                          m
!    qs_{dnu0+k}(x,c) =  \sum    a_{j,k} Q_{dnu0+k+2j} (x)
!                         j=-m
!
!  of the angular prolate spheroidal wave functions of the first and
!  second kinds of characteristic exponent dnu0+k for k=0,1,...,n.
!  An appropriate value for m is determined by this routine.
!
!  WARNING: this routine is quite slow for large values of n and c
!  as it uses involves two nonsymmetric dense eigensolves.
!
!  Input parameters:
!    c - the bandlimit for the prolate functions
!    dnu0 - the value of the parameter dnu0, which must be in the
!      range -1/2 < dnu0 < 1/2
!    n - the parameter controlling the number of expansion to construct
!
!  Output parameters:
!    m - a parameter controlling the number of terms in the expansions
!    chis - an zero-indexed array whose jth entry contains the value of
!        chi_{dnu0+j}(c)
!    coefs - a matrix dimensioned (-m:m,0:n) whose jth column
!      coefficients in the expansions of ps_{dnu0+k}
!

double precision, allocatable :: amatr(:,:), wi(:), wr(:), vr(:,:), vl(:,:),work(:)
integer, allocatable          :: idxs(:)

!
!  Determine the dimension of the matrix and allocate the output matrices
!

nn     = ceiling((n+0.0d0)/2.0d0)
m      = 1.2*max(ceiling((n+sqrt(n*c))),20)
l      = 2*m+1

allocate(chis(0:2*nn),coefs(-m:m,0:2*nn))


!
!  Form the first nonsymmetric tridiagonal matrix
!


allocate(amatr(-m:m,-m:m))
amatr  = 0 


dnu = dnu0
do i=-m,m

if (i .ne. -m) then
amatr(i,i-1)   = c**2 * ( (dnu + 2*i)*(dnu+2*i-1)) / ( (2*dnu+4*i-1)* (2*dnu+4*i-3) )
endif


amatr(i,i)   = (dnu+2*i)*(dnu+2*i+1) + c**2 * ( 2 * (dnu + 2*i)*(dnu+2*i+1)-1) / &
  ( (2*dnu+4*i+3)* (2*dnu+4*i-1) )

if (i .ne. m) then
amatr(i,i+1)   = c**2 * ((dnu + 2*i+2)*(dnu+2*i+1)) / ( (2*dnu+4*i+3)* (2*dnu+4*i+5) )
endif

end do


!
!  Compute its eigenvalues and eigenvectors via LAPACK
!
allocate(wr(l),wi(l),vr(-m:m,l))
wr = 0
wi = 0

lwork=-1
call dgeev('N','V',l,amatr,l,wr,wi,vl,l,vr,l,ww,lwork,info)
lwork = ww+1000
allocate(work(lwork))
call dgeev('N','V',l,amatr,l,wr,wi,vl,l,vr,l,work,lwork,info)



!
! Sort the list of indices and copy them out to the appropriate place
!


allocate(idxs(l))
do i=1,l
idxs(i) = i
end do

call quicksort2(l,wr,idxs)

do i=0,nn
chis(2*i) = wr(2*i+1)

do j=-m,m
coefs(j,2*i) = vr(j,idxs(2*i+1))
end do

end do


!
!  Now form the second matrix 
!

amatr  = 0 
dnu = dnu0+1.0d0
do i=-m,m

if (i .ne. -m) then
amatr(i,i-1)   = c**2 * ( (dnu + 2*i)*(dnu+2*i-1)) / ( (2*dnu+4*i-1)* (2*dnu+4*i-3) )
endif


amatr(i,i)   = (dnu+2*i)*(dnu+2*i+1) + c**2 * ( 2 * (dnu + 2*i)*(dnu+2*i+1)-1) / &
  ( (2*dnu+4*i+3)* (2*dnu+4*i-1) )

if (i .ne. m) then
amatr(i,i+1)   = c**2 * ((dnu + 2*i+2)*(dnu+2*i+1)) / ( (2*dnu+4*i+3)* (2*dnu+4*i+5) )
endif

end do


!
!  Compute its eigenvalues and eigenvectors via LAPACK and copy them
!  into the correct place in the output arrays
!

call dgeev('N','V',l,amatr,l,wr,wi,vl,l,vr,l,work,lwork,info)

!
! Sort the list of indices and copy them out to the appropriate place
!


do i=1,l
idxs(i) = i
end do

call quicksort2(l,wr,idxs)

do i=0,nn-1
chis(2*i+1)  = wr(2*i+2)

do j=-m,m
coefs(j,2*i+1) = vr(j,idxs(2*i+2))
end do

end do



!
!  Normalize the coefficient expansions
!

do j=0,n
dnu  = dnu0 + j
rn   = 0

do i=-m,m
rn = rn + coefs(i,j)**2 * (2*dnu+1.0d0) / (2*dnu0+4*i+1.0d0)
end do


coefs(:,j) = coefs(:,j) / sqrt(rn)

end do


do j=0,nn
coefs(:,2*j) = coefs(:,2*j) * sign(1.0d0,coefs(0,2*j))
end do

do j=0,nn-1
coefs(:,2*j+1) = coefs(:,2*j+1) *  sign(1.0d0,coefs(0,2*j+1))
end do

end subroutine



! subroutine prolates_noninteger_all_eval(dnu0,n,m,coefs,j,x,valps,valqs)
! implicit double precision (a-h,o-z)
! double precision  :: coefs(-m:m,0:n)
! !
! !  Evaluate 
! !
! !     ps_{dnu0+j} (x,c)    
! !
! !  and     
! !
! !     qs_{dnu0+j}(x,c)
! !
! !  given the coefficient expansions constructed by prolates_noninteger_all.
! !
! !
! !  Input parameters:
! !    dnu0 - the
! !    n - the number of 
! !    m - parameter returned by prolates_noninteger_all which indicates the
! !      size of the expansions
! !    j - index of the angular prolate function to evaluate
! !    x - the point in the interval (-1,1) at which to evaluate the functions
! !
! !  Output parameters:
! !    valps - the value of ps_{dnu0+j}(x,c)
! !    valqs - the value of qs_{dnu0+j}(x,c)
! !



! valps = 0
! valqs = 0

! do i=-m,m
! call prolates_lege(dnu0+2*i,x,valp,valq)

! valps = valps  + coefs(i,j) * valp
! valqs = valqs  + coefs(i,j) * valq

! end do

! end subroutine


! subroutine prolates_noninteger_all_evalder(dnu0,n,m,coefs,j,x,valps,derps,valqs,derqs)
! implicit double precision (a-h,o-z)
! double precision  :: coefs(-m:m,0:n)
! !
! !  Evaluate 
! !
! !     ps_{dnu0+j} (x,c)    
! !
! !  and     
! !
! !     qs_{dnu0+j}(x,c)
! !
! !  and their derivatives given the coefficient expansions constructed by
! !   prolates_noninteger_all.
! !
! !
! !  Input parameters:
! !    dnu0 - the
! !    n - the number of 
! !    m - parameter returned by prolates_noninteger_all which indicates the
! !      size of the expansions
! !    j - index of the angular prolate function to evaluate
! !    x - the point in the interval (-1,1) at which to evaluate the functions
! !
! !  Output parameters:
! !    valps - the value of ps_{dnu0+j}(x,c)
! !    valqs - the value of qs_{dnu0+j}(x,c)
! !



! valps = 0
! valqs = 0

! derps = 0
! derqs = 0

! do i=-m,m
! call prolates_legeder(dnu0+2*i,x,valp,derp,valq,derq)

! valps = valps  + coefs(i,j) * valp
! derps = derps  + coefs(i,j) * derp
! valqs = valqs  + coefs(i,j) * valq
! derqs = derqs  + coefs(i,j) * derq

! end do

! end subroutine



! subroutine prolates_lege(dnu,x,valp,valq)
! implicit double precision (a-h,o-z)
! !
! !  Evaluate the Legendre functions of the first and second kinds of real order
! !  dnu and argument 0 <= x < 1.  When dnu is a negative integer, the function
! !  of the second kind is undefined and the value returned by this routine is
! !  undefined.
! !
! !  Input parameters:
! !    dnu - the 
! !    x - the point in the interval [0,1) at which to evaluate the Legendre
! !      functions
! !
! !  Output parameters:
! !    valp - the value of P_dnu(x)
! !    valq - the value of Q_dnu(x)
! !

! data pi / 3.14159265358979323846264338327950288d0 /

! t = acos(x)

! call alegendre_eval_init(dd)

! ! Handle the case of positive dnu
! if (dnu .ge. 0) then

! call alegendre_eval(dnu,0.0d0,t,alpha,alphader,vallogp,vallogq,valp0,valq0)

! valp = valp0 / sqrt(sin(t)) 
! valp = valp  / sqrt(dnu+0.5d0)

! valq = valq0 / sqrt(sin(t)) 
! valq = valq  / sqrt(dnu+0.5d0) * pi/2

! return
! endif

! ! Handle the case in which dnu < -1
! if (dnu .le. -1) then

! dnu0 = -dnu-1
! call alegendre_eval(dnu0,0.0d0,t,alpha,alphader,vallogp,vallogq,valp0,valq0)
! valp0 = valp0 / sqrt(sin(t)) 
! valp0 = valp0  / sqrt(dnu0+0.5d0)

! valq0 = valq0 / sqrt(sin(t)) 
! valq0 = valq0  / sqrt(dnu0+0.5d0)*pi/2

! valp = valp0
! valq = valq0 - pi * 1/tan(pi*dnu0)*valp0

! return
! return

! endif

! ! Handle the case in which -1 < dnu < 0

! dnu1 = dnu+1
! call alegendre_eval(dnu1,0.0d0,t,alpha,alphader,vallogp,vallogq,valp1,valq1)

! valp1 = valp1 / sqrt(sin(t)) 
! valp1 = valp1  / sqrt(dnu1+0.5d0)

! valq1 = valq1 / sqrt(sin(t)) 
! valq1 = valq1  / sqrt(dnu1+0.5d0)*pi/2

! dnu2 = dnu+2
! call alegendre_eval(dnu2,0.0d0,t,alpha,alphader,vallogp,vallogq,valp2,valq2)
! valp2 = valp2 / sqrt(sin(t)) 
! valp2 = valp2  / sqrt(dnu2+0.5d0)

! valq2 = valq2 / sqrt(sin(t)) 
! valq2 = valq2  / sqrt(dnu2+0.5d0)*pi/2


! valp = (2*dnu+3)/(dnu+1) * x*valp1 - (dnu+2)/(dnu+1) * valp2
! valq = (2*dnu+3)/(dnu+1) * x*valq1 - (dnu+2)/(dnu+1) * valq2

! end subroutine



! subroutine prolates_legeder(dnu,x,valp,derp,valq,derq)
! implicit double precision (a-h,o-z)
! !
! !  Evaluate the Legendre functions of the first and second kinds of real order
! !  dnu and argument 0 <= x < 1.  When dnu is a negative integer, the function
! !  of the second kind is undefined and the values returned in valq and derq
! !  are undefined.
! !
! !  Input parameters:
! !    dnu - the degree of the Legendre functions to evaluation
! !    x - the point in the interval [0,1) at which to evaluate them
! !
! !  Output parameters:
! !    valp - the value of P_dnu(x)
! !    derp - the value of P_dnu'(x)
! !    valq - the value of Q_dnu(x)
! !    derq - the value of Q_dnu'(x)
! !

! call prolates_lege(dnu,x,valp,valq)
! call prolates_lege(dnu+1.0d0,x,valp1,valq1)

! derp = (dnu+1)*valp1 - (dnu+1)*x*valp
! derp = derp / (x**2-1)

! derq = (dnu+1)*valq1 - (dnu+1)*x*valq
! derq = derq / (x**2-1)

! end subroutine


subroutine prolate_trieigen(n,as,bs,j,chi,y)
implicit double precision (a-h,o-z)
double precision             :: as(n),bs(n),y(n)
!
!  Calculate a single eigenvalue and a corresponding eigenvector of a
!  symmetric tridiagonal matrix.
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
call prolates_sturm_bisection(eps,n,as,bs,j,chi)


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

!
! Solve the system ( A - chi * I ) x = y 
!

as0 = as - chi0
bs0 = bs


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


!
!  Normalize the eigenvector
!  

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




end subroutine



subroutine prolates_sturm_bisection(eps,n,as,bs,j,dlambda)
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

do while ( abs(x2-x1) .gt. eps * (abs(x1)+abs(x2)) )
x       = (x1+x2)/2

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

subroutine prolates_trisolve(n,as,bs,cs,y)
implicit double precision (a-h,o-z)
double precision :: as(n-1),bs(n),cs(n-1),y(n)
!
!  Solve the linear system of equations
!
!     A x = y
!
!  where A is an (n,n) tridiagonal matrix.  
!
!  Input parameters:
!    n - the dimension of the matrix A
!    as - the subdiagonal entries
!    bs - the diagonal entries
!    cs - the superdiagonal entries
!    y - upon input, this vector should be the right-hand
!       side
!
!  Output parameters:
!    y - upon return, this vector will contain the solution
!
!

!
!  Eliminate the subdiagonal to produce an upper triangular matrix
!


do i=1,n-1
bs(i+1) = bs(i+1) - as(i)/bs(i) * cs(i)
y(i+1)  = y(i+1)  - as(i)/bs(i) * y(i)
end do

!
!  Do back substitution 
!

y(n) = y(n) / bs(n)


do i=n-1,1,-1
y(i) = ( y(i) - cs(i)*y(i+1) ) / bs(i)
end do

end subroutine



subroutine prolates_integer_bound(c,n,chi1,chi2)
implicit double precision (a-h,o-z)
data pi / 3.14159265358979323846264338327950288d0 /
!
!  Use the bounds
!
!    A. Bonami and A. Kaouri, "Uniform Approximation and Explicit Estimates 
!       for the Prolate Spheroidal Wave Functions,"  Constructive Approximation
!       43 (2016), pages 15-45
!
!  to find an interval containing a Sturm-Liouville eigenvalue.  The midpoint 
!  of the interval approximates the eigenvalue.
!
!
!  Input parameters:
!    c - the bandlimit
!    n - the characteristic exponent of the prolate function
!
!  Output parameters:
!    (chi1, chi2) - an interval containing chi
!
x1 = 2*c / ( pi * n )
x2 = 2*c / ( pi * (n+1) )
call evalphi(x1,val1)
call evalphi(x2,val2)
chi1 = ( c / val1 )**2
chi2 = ( c / val2 )**2
end subroutine


subroutine prolates_integer_coefficients(c,n,ncoefs,coefsps,coefsqs,c1,c2)
implicit double precision (a-h,o-z)
double precision       :: coefsps(ncoefs),coefsqs(ncoefs)

!
!  Compute coefficients c1 and c2 such that the square of the modulus of 
!
!    c1 ps_n(x) + i c2 qs_n(x)                                                            (3)
!
!  is completely monotone and the pair {c1 ps_n, c2 qs_n} has Wronskian 1.
!  
!  Input parameters:
!    c - the bandlimti
!    n -  the characteristic exponent 
!    ncoefs - the number of coefficients 
!    coefsps - the coefficients for ps constructed by prolates_integer
!    coefsqs - the coefficients for qs constructed by prolates_integer
!
!  Output parameters:
!    c1 - the first coefficient in (3)
!    c2 - the second coefficient in (3)
!
!

data pi / 3.14159265358979323846264338327950288d0 /

k = n/2

call prolates_integer_evalder(n,ncoefs,coefsps,coefsqs,0.0d0,valps,derps,valqs,derqs)
aa = sum(coefsps)


if (k*2 .eq. n) then

dk1 = (-1)**k     * coefsps(1) * 1.0d0 * 1/valps
dk2 = (-1)**(k+1) * valps * 1.0d0/c * 1/coefsps(1)



c1  = dk1 + pi/2*dk2
c2  = dk2


W   = aa**2/c*(1-pi/2*1/c*valps**2/coefsps(1)**2)

c1  = c1 / sqrt(abs(W))
c2  = c2 / sqrt(abs(W))

else

dk1 = (-1)**k * c/3 * coefsps(1) * 1/aa * 1/derps
dk2 = (-1)**(k+1) * 3/c**2 * 1/coefsps(1) * derps/aa

c1  = dk1 + pi/2*dk2
c2  = dk2

W   = c1*c2*aa**2

c1  = c1/sqrt(abs(W))
c2  = c2/sqrt(abs(W))


endif

end subroutine



subroutine evalphi(x,val)
implicit double precision (a-h,o-z)
integer             :: nints, k
double precision    :: xscheb(00000030)
double precision    :: ab(00000002,00000089)
double precision    :: vals(00000030,00000089)
data nints      / 00000089 /
data k          / 00000030 /
data xscheb     /   -0.100000000000000D+01, &
                    -0.994137957154360D+00, &
                    -0.976620555710087D+00, &
                    -0.947653171182802D+00, &
                    -0.907575419670957D+00, &
                    -0.856857176167589D+00, &
                    -0.796093065705644D+00, &
                    -0.725995491923131D+00, &
                    -0.647386284781828D+00, &
                    -0.561187065362382D+00, &
                    -0.468408440699790D+00, &
                    -0.370138155339914D+00, &
                    -0.267528338529221D+00, &
                    -0.161781996552765D+00, &
                    -0.541389085854175D-01, &
                     0.541389085854175D-01, &
                     0.161781996552765D+00, &
                     0.267528338529221D+00, &
                     0.370138155339914D+00, &
                     0.468408440699790D+00, &
                     0.561187065362382D+00, &
                     0.647386284781828D+00, &
                     0.725995491923131D+00, &
                     0.796093065705644D+00, &
                     0.856857176167589D+00, &
                     0.907575419670957D+00, &
                     0.947653171182802D+00, &
                     0.976620555710087D+00, &
                     0.994137957154360D+00, &
                     0.100000000000000D+01  /
data ab         /    0.000000000000000D+00, &
                     0.340724276794204D+00, &
                     0.340724276794204D+00, &
                     0.568840247322897D+00, &
                     0.568840247322897D+00, &
                     0.728492428264635D+00, &
                     0.728492428264635D+00, &
                     0.835864576284616D+00, &
                     0.835864576284616D+00, &
                     0.904190331858824D+00, &
                     0.904190331858824D+00, &
                     0.945568826231068D+00, &
                     0.945568826231068D+00, &
                     0.969687034264429D+00, &
                     0.969687034264429D+00, &
                     0.983360689439900D+00, &
                     0.983360689439900D+00, &
                     0.990962450186940D+00, &
                     0.990962450186940D+00, &
                     0.995130134588190D+00, &
                     0.995130134588190D+00, &
                     0.997392031722123D+00, &
                     0.997392031722123D+00, &
                     0.998610274706748D+00, &
                     0.998610274706748D+00, &
                     0.999262501510237D+00, &
                     0.999262501510237D+00, &
                     0.999610002498256D+00, &
                     0.999610002498256D+00, &
                     0.999794397495666D+00, &
                     0.999794397495666D+00, &
                     0.999891902758736D+00, &
                     0.999891902758736D+00, &
                     0.999943305068899D+00, &
                     0.999943305068899D+00, &
                     0.999970329856727D+00, &
                     0.999970329856727D+00, &
                     0.999984503723871D+00, &
                     0.999984503723871D+00, &
                     0.999991921297247D+00, &
                     0.999991921297247D+00, &
                     0.999995795376715D+00, &
                     0.999995795376715D+00, &
                     0.999997815055288D+00, &
                     0.999997815055288D+00, &
                     0.999998866211903D+00, &
                     0.999998866211903D+00, &
                     0.999999412448299D+00, &
                     0.999999412448299D+00, &
                     0.999999695895383D+00, &
                     0.999999695895383D+00, &
                     0.999999842783325D+00, &
                     0.999999842783325D+00, &
                     0.999999918809485D+00, &
                     0.999999918809485D+00, &
                     0.999999958113655D+00, &
                     0.999999958113655D+00, &
                     0.999999978411285D+00, &
                     0.999999978411285D+00, &
                     0.999999988882872D+00, &
                     0.999999988882872D+00, &
                     0.999999994280052D+00, &
                     0.999999994280052D+00, &
                     0.999999997059335D+00, &
                     0.999999997059335D+00, &
                     0.999999998489322D+00, &
                     0.999999998489322D+00, &
                     0.999999999224490D+00, &
                     0.999999999224490D+00, &
                     0.999999999602160D+00, &
                     0.999999999602160D+00, &
                     0.100000000000000D+01, &
                     0.100000000000000D+01, &
                     0.100000000030653D+01, &
                     0.100000000030653D+01, &
                     0.100000000059767D+01, &
                     0.100000000059767D+01, &
                     0.100000000116456D+01, &
                     0.100000000116456D+01, &
                     0.100000000226756D+01, &
                     0.100000000226756D+01, &
                     0.100000000441200D+01, &
                     0.100000000441200D+01, &
                     0.100000000857774D+01, &
                     0.100000000857774D+01, &
                     0.100000001666298D+01, &
                     0.100000001666298D+01, &
                     0.100000003234095D+01, &
                     0.100000003234095D+01, &
                     0.100000006271188D+01, &
                     0.100000006271188D+01, &
                     0.100000012148373D+01, &
                     0.100000012148373D+01, &
                     0.100000023508739D+01, &
                     0.100000023508739D+01, &
                     0.100000045441467D+01, &
                     0.100000045441467D+01, &
                     0.100000087730928D+01, &
                     0.100000087730928D+01, &
                     0.100000169157897D+01, &
                     0.100000169157897D+01, &
                     0.100000325708067D+01, &
                     0.100000325708067D+01, &
                     0.100000626201377D+01, &
                     0.100000626201377D+01, &
                     0.100001201975768D+01, &
                     0.100001201975768D+01, &
                     0.100002303106681D+01, &
                     0.100002303106681D+01, &
                     0.100004404556361D+01, &
                     0.100004404556361D+01, &
                     0.100008405915291D+01, &
                     0.100008405915291D+01, &
                     0.100016005848248D+01, &
                     0.100016005848248D+01, &
                     0.100030401180082D+01, &
                     0.100030401180082D+01, &
                     0.100057586365868D+01, &
                     0.100057586365868D+01, &
                     0.100108758093291D+01, &
                     0.100108758093291D+01, &
                     0.100204745961269D+01, &
                     0.100204745961269D+01, &
                     0.100384149810881D+01, &
                     0.100384149810881D+01, &
                     0.100718271714454D+01, &
                     0.100718271714454D+01, &
                     0.101338623917495D+01, &
                     0.101338623917495D+01, &
                     0.102488241059831D+01, &
                     0.102488241059831D+01, &
                     0.104619936395725D+01, &
                     0.104619936395725D+01, &
                     0.108593191491656D+01, &
                     0.108593191491656D+01, &
                     0.116096438943856D+01, &
                     0.116096438943856D+01, &
                     0.130636851576710D+01, &
                     0.130636851576710D+01, &
                     0.160110846571051D+01, &
                     0.160110846571051D+01, &
                     0.224232104640424D+01, &
                     0.224232104640424D+01, &
                     0.378278302185347D+01, &
                     0.378278302185347D+01, &
                     0.795735518311802D+01, &
                     0.795735518311802D+01, &
                     0.207404790043730D+02, &
                     0.207404790043730D+02, &
                     0.640865287451800D+02, &
                     0.640865287451800D+02, &
                     0.221922177256830D+03, &
                     0.221922177256830D+03, &
                     0.822178030554441D+03, &
                     0.822178030554441D+03, &
                     0.316103154309294D+04, &
                     0.316103154309294D+04, &
                     0.123921068668076D+05, &
                     0.123921068668076D+05, &
                     0.490677310051761D+05, &
                     0.490677310051761D+05, &
                     0.195272873323581D+06, &
                     0.195272873323581D+06, &
                     0.779098734146655D+06, &
                     0.779098734146655D+06, &
                     0.311241276054450D+07, &
                     0.311241276054450D+07, &
                     0.124416900323221D+08, &
                     0.124416900323221D+08, &
                     0.497508414515967D+08, &
                     0.497508414515967D+08, &
                     0.198971531793825D+09, &
                     0.198971531793825D+09, &
                     0.795822462479599D+09, &
                     0.795822462479599D+09, &
                     0.318316252374976D+10, &
                     0.318316252374976D+10, &
                     0.127323954462943D+11  /
data vals       /    0.000000000000000D+00, &
                     0.156870644537564D-02, &
                     0.625637675852292D-02, &
                     0.140075451901850D-01, &
                     0.247294260423736D-01, &
                     0.382915786747984D-01, &
                     0.545257764494490D-01, &
                     0.732263702473541D-01, &
                     0.941514542880642D-01, &
                     0.117025114164934D+00, &
                     0.141540961954919D+00, &
                     0.167367045377864D+00, &
                     0.194152069504357D+00, &
                     0.221532709952876D+00, &
                     0.249141650027058D+00, &
                     0.276615865165335D+00, &
                     0.303604625692188D+00, &
                     0.329776702856362D+00, &
                     0.354826341232667D+00, &
                     0.378477689124385D+00, &
                     0.400487536093320D+00, &
                     0.420646368820110D+00, &
                     0.438777901210697D+00, &
                     0.454737346325793D+00, &
                     0.468408768617484D+00, &
                     0.479701884991848D+00, &
                     0.488548677906897D+00, &
                     0.494900151641341D+00, &
                     0.498723513200566D+00, &
                     0.500000000000000D+00, &
                     0.500000000000000D+00, &
                     0.500853778538474D+00, &
                     0.503401114694551D+00, &
                     0.507600325417093D+00, &
                     0.513382985351205D+00, &
                     0.520655479003414D+00, &
                     0.529301135558517D+00, &
                     0.539182901237797D+00, &
                     0.550146483556812D+00, &
                     0.562023878967490D+00, &
                     0.574637172528953D+00, &
                     0.587802478507256D+00, &
                     0.601333877345595D+00, &
                     0.615047199962687D+00, &
                     0.628763516461008D+00, &
                     0.642312203277161D+00, &
                     0.655533489310994D+00, &
                     0.668280415080065D+00, &
                     0.680420176090877D+00, &
                     0.691834858778758D+00, &
                     0.702421611241751D+00, &
                     0.712092319031848D+00, &
                     0.720772876915114D+00, &
                     0.728402160238455D+00, &
                     0.734930804685426D+00, &
                     0.740319901720521D+00, &
                     0.744539710185425D+00, &
                     0.747568473655174D+00, &
                     0.749391419489967D+00, &
                     0.750000000000000D+00, &
                     0.750000000000000D+00, &
                     0.750425495981697D+00, &
                     0.751694871246134D+00, &
                     0.753786973885788D+00, &
                     0.756667164070157D+00, &
                     0.760288209486762D+00, &
                     0.764591490759774D+00, &
                     0.769508470923087D+00, &
                     0.774962373191870D+00, &
                     0.780870003921178D+00, &
                     0.787143653235050D+00, &
                     0.793693004690136D+00, &
                     0.800426987639649D+00, &
                     0.807255511577988D+00, &
                     0.814091030307099D+00, &
                     0.820849894672226D+00, &
                     0.827453465090735D+00, &
                     0.833828968276528D+00, &
                     0.839910095581180D+00, &
                     0.845637352461984D+00, &
                     0.850958179137678D+00, &
                     0.855826871093597D+00, &
                     0.860204334538572D+00, &
                     0.864057716159240D+00, &
                     0.867359948653389D+00, &
                     0.870089253715875D+00, &
                     0.872228642587268D+00, &
                     0.873765451140817D+00, &
                     0.874690941948345D+00, &
                     0.875000000000000D+00, &
                     0.875000000000000D+00, &
                     0.875207670246404D+00, &
                     0.875827369502558D+00, &
                     0.876849244959802D+00, &
                     0.878257144242098D+00, &
                     0.880029007953397D+00, &
                     0.882137397994942D+00, &
                     0.884550141703579D+00, &
                     0.887231067935316D+00, &
                     0.890140808477038D+00, &
                     0.893237636725231D+00, &
                     0.896478315456880D+00, &
                     0.899818926687687D+00, &
                     0.903215658944623D+00, &
                     0.906625530588336D+00, &
                     0.910007031874398D+00, &
                     0.913320672983651D+00, &
                     0.916529430020801D+00, &
                     0.919599085734005D+00, &
                     0.922498466237889D+00, &
                     0.925199579163028D+00, &
                     0.927677662287408D+00, &
                     0.929911154750530D+00, &
                     0.931881605357568D+00, &
                     0.933573534210990D+00, &
                     0.934974264921198D+00, &
                     0.936073744898382D+00, &
                     0.936864370662675D+00, &
                     0.937340833686904D+00, &
                     0.937500000000000D+00, &
                     0.937500000000000D+00, &
                     0.937601212200564D+00, &
                     0.937903326065887D+00, &
                     0.938401807057172D+00, &
                     0.939089208677886D+00, &
                     0.939955335049811D+00, &
                     0.940987460664256D+00, &
                     0.942170599900474D+00, &
                     0.943487817394713D+00, &
                     0.944920569238999D+00, &
                     0.946449064325779D+00, &
                     0.948052634945016D+00, &
                     0.949710105971714D+00, &
                     0.951400152619434D+00, &
                     0.953101637726706D+00, &
                     0.954793920822924D+00, &
                     0.956457132716584D+00, &
                     0.958072410988815D+00, &
                     0.959622093490942D+00, &
                     0.961089868676272D+00, &
                     0.962460883293025D+00, &
                     0.963721809585391D+00, &
                     0.964860875655812D+00, &
                     0.965867863994988D+00, &
                     0.966734084339970D+00, &
                     0.967452327914788D+00, &
                     0.968016810668274D+00, &
                     0.968423113268805D+00, &
                     0.968668125271905D+00, &
                     0.968750000000000D+00, &
                     0.968750000000000D+00, &
                     0.968799555451128D+00, &
                     0.968947513173483D+00, &
                     0.969191761307715D+00, &
                     0.969528826447254D+00, &
                     0.969953941723727D+00, &
                     0.970461139281424D+00, &
                     0.971043364434913D+00, &
                     0.971692608220813D+00, &
                     0.972400054599792D+00, &
                     0.973156238251600D+00, &
                     0.973951208741065D+00, &
                     0.974774696816146D+00, &
                     0.975616278723896D+00, &
                     0.976465534684710D+00, &
                     0.977312198034375D+00, &
                     0.978146292009619D+00, &
                     0.978958251698371D+00, &
                     0.979739029283878D+00, &
                     0.980480181367008D+00, &
                     0.981173937840043D+00, &
                     0.981813252494995D+00, &
                     0.982391836264968D+00, &
                     0.982904174698123D+00, &
                     0.983345531920735D+00, &
                     0.983711943917499D+00, &
                     0.984000204390282D+00, &
                     0.984207846690560D+00, &
                     0.984333125298210D+00, &
                     0.984375000000000D+00, &
                     0.984375000000000D+00, &
                     0.984399396160528D+00, &
                     0.984472249116780D+00, &
                     0.984592558344053D+00, &
                     0.984758676315446D+00, &
                     0.984968337980173D+00, &
                     0.985218700966069D+00, &
                     0.985506395482713D+00, &
                     0.985827582668972D+00, &
                     0.986178019937187D+00, &
                     0.986553131720677D+00, &
                     0.986948083935006D+00, &
                     0.987357860417799D+00, &
                     0.987777339616333D+00, &
                     0.988201369844987D+00, &
                     0.988624841532887D+00, &
                     0.989042755022662D+00, &
                     0.989450282660549D+00, &
                     0.989842824133333D+00, &
                     0.990216054256078D+00, &
                     0.990565962693996D+00, &
                     0.990888885409118D+00, &
                     0.991181527952667D+00, &
                     0.991440981068131D+00, &
                     0.991664729412865D+00, &
                     0.991850654524464D+00, &
                     0.991997033421252D+00, &
                     0.992102534397457D+00, &
                     0.992166211617174D+00, &
                     0.992187500000000D+00, &
                     0.992187500000000D+00, &
                     0.992199564701738D+00, &
                     0.992235597487660D+00, &
                     0.992295117105718D+00, &
                     0.992377330426896D+00, &
                     0.992481145667613D+00, &
                     0.992605190478951D+00, &
                     0.992747834494821D+00, &
                     0.992907215834407D+00, &
                     0.993081270971183D+00, &
                     0.993267767313374D+00, &
                     0.993464337790494D+00, &
                     0.993668516708253D+00, &
                     0.993877776120314D+00, &
                     0.994089561969982D+00, &
                     0.994301329278025D+00, &
                     0.994510575693996D+00, &
                     0.994714872787854D+00, &
                     0.994911894536189D+00, &
                     0.995099442553492D+00, &
                     0.995275467733865D+00, &
                     0.995438088102193D+00, &
                     0.995585602824836D+00, &
                     0.995716502494529D+00, &
                     0.995829475975216D+00, &
                     0.995923414258245D+00, &
                     0.995997411925202D+00, &
                     0.996050766914894D+00, &
                     0.996082979332283D+00, &
                     0.996093750000000D+00, &
                     0.996093750000000D+00, &
                     0.996099735926716D+00, &
                     0.996117615249964D+00, &
                     0.996147153846162D+00, &
                     0.996187965640763D+00, &
                     0.996239518716225D+00, &
                     0.996301143687010D+00, &
                     0.996372044170439D+00, &
                     0.996451309140259D+00, &
                     0.996537926912693D+00, &
                     0.996630800483305D+00, &
                     0.996728763907884D+00, &
                     0.996830599402120D+00, &
                     0.996935054823561D+00, &
                     0.997040861195422D+00, &
                     0.997146749935564D+00, &
                     0.997251469465604D+00, &
                     0.997353800895114D+00, &
                     0.997452572504696D+00, &
                     0.997546672790015D+00, &
                     0.997635061877289D+00, &
                     0.997716781179687D+00, &
                     0.997790961233293D+00, &
                     0.997856827729311D+00, &
                     0.997913705842708D+00, &
                     0.997961023040713D+00, &
                     0.997998310629159D+00, &
                     0.998025204350432D+00, &
                     0.998041444372923D+00, &
                     0.998046875000000D+00, &
                     0.998046875000000D+00, &
                     0.998049851511796D+00, &
                     0.998058742586211D+00, &
                     0.998073433428167D+00, &
                     0.998093734660832D+00, &
                     0.998119385209223D+00, &
                     0.998150056260339D+00, &
                     0.998185356224799D+00, &
                     0.998224836606101D+00, &
                     0.998267998666652D+00, &
                     0.998314300764893D+00, &
                     0.998363166225478D+00, &
                     0.998413991594788D+00, &
                     0.998466155127211D+00, &
                     0.998519025343897D+00, &
                     0.998571969505233D+00, &
                     0.998624361841312D+00, &
                     0.998675591391560D+00, &
                     0.998725069315822D+00, &
                     0.998772235554909D+00, &
                     0.998816564739610D+00, &
                     0.998857571273503D+00, &
                     0.998894813546929D+00, &
                     0.998927897276333D+00, &
                     0.998956478003605D+00, &
                     0.998980262831026D+00, &
                     0.998999011505083D+00, &
                     0.999012536991378D+00, &
                     0.999020705697848D+00, &
                     0.999023437500000D+00, &
                     0.999023437500000D+00, &
                     0.999024919735907D+00, &
                     0.999029347489525D+00, &
                     0.999036664179983D+00, &
                     0.999046776438923D+00, &
                     0.999059555492600D+00, &
                     0.999074839062091D+00, &
                     0.999092433747580D+00, &
                     0.999112117854068D+00, &
                     0.999133644607866D+00, &
                     0.999156745706236D+00, &
                     0.999181135136454D+00, &
                     0.999206513195672D+00, &
                     0.999232570639218D+00, &
                     0.999258992882586D+00, &
                     0.999285464181443D+00, &
                     0.999311671714603D+00, &
                     0.999337309497363D+00, &
                     0.999362082057005D+00, &
                     0.999385707808976D+00, &
                     0.999407922081494D+00, &
                     0.999428479748342D+00, &
                     0.999447157444465D+00, &
                     0.999463755356508D+00, &
                     0.999478098599819D+00, &
                     0.999490038213513D+00, &
                     0.999499451823843D+00, &
                     0.999506244040842D+00, &
                     0.999510346661266D+00, &
                     0.999511718750000D+00, &
                     0.999511718750000D+00, &
                     0.999512457577829D+00, &
                     0.999514664687875D+00, &
                     0.999518312095978D+00, &
                     0.999523353613450D+00, &
                     0.999529725516314D+00, &
                     0.999537347466163D+00, &
                     0.999546123666840D+00, &
                     0.999555944237090D+00, &
                     0.999566686775578D+00, &
                     0.999578218091286D+00, &
                     0.999590396069337D+00, &
                     0.999603071639849D+00, &
                     0.999616090815476D+00, &
                     0.999629296761943D+00, &
                     0.999642531865202D+00, &
                     0.999655639758873D+00, &
                     0.999668467276486D+00, &
                     0.999680866294891D+00, &
                     0.999692695438083D+00, &
                     0.999703821614870D+00, &
                     0.999714121369334D+00, &
                     0.999723482029994D+00, &
                     0.999731802651900D+00, &
                     0.999738994755223D+00, &
                     0.999744982873646D+00, &
                     0.999749704935037D+00, &
                     0.999753112504267D+00, &
                     0.999755170922228D+00, &
                     0.999755859375000D+00, &
                     0.999755859375000D+00, &
                     0.999756227882789D+00, &
                     0.999757328759903D+00, &
                     0.999759148133935D+00, &
                     0.999761663104288D+00, &
                     0.999764842068477D+00, &
                     0.999768645171374D+00, &
                     0.999773024869945D+00, &
                     0.999777926604094D+00, &
                     0.999783289562453D+00, &
                     0.999789047530282D+00, &
                     0.999795129805246D+00, &
                     0.999801462165593D+00, &
                     0.999807967874264D+00, &
                     0.999814568701766D+00, &
                     0.999821185950214D+00, &
                     0.999827741460881D+00, &
                     0.999834158587890D+00, &
                     0.999840363121479D+00, &
                     0.999846284145518D+00, &
                     0.999851854815919D+00, &
                     0.999857013049070D+00, &
                     0.999861702112784D+00, &
                     0.999865871116153D+00, &
                     0.999869475399231D+00, &
                     0.999872476828157D+00, &
                     0.999874844005868D+00, &
                     0.999876552412162D+00, &
                     0.999877584489063D+00, &
                     0.999877929687500D+00, &
                     0.999877929687500D+00, &
                     0.999878113570001D+00, &
                     0.999878662910783D+00, &
                     0.999879570822176D+00, &
                     0.999880825932507D+00, &
                     0.999882412545945D+00, &
                     0.999884310862670D+00, &
                     0.999886497255810D+00, &
                     0.999888944600670D+00, &
                     0.999891622650896D+00, &
                     0.999894498455438D+00, &
                     0.999897536809476D+00, &
                     0.999900700731845D+00, &
                     0.999903951961024D+00, &
                     0.999907251461386D+00, &
                     0.999910559931136D+00, &
                     0.999913838303351D+00, &
                     0.999917048231596D+00, &
                     0.999920152551953D+00, &
                     0.999923115713874D+00, &
                     0.999925904173129D+00, &
                     0.999928486741348D+00, &
                     0.999930834888187D+00, &
                     0.999932922994067D+00, &
                     0.999934728553552D+00, &
                     0.999936232331784D+00, &
                     0.999937418478523D+00, &
                     0.999938274606203D+00, &
                     0.999938791839466D+00, &
                     0.999938964843750D+00, &
                     0.999938964843750D+00, &
                     0.999939056628282D+00, &
                     0.999939330835380D+00, &
                     0.999939784041513D+00, &
                     0.999940410593754D+00, &
                     0.999941202688343D+00, &
                     0.999942150478936D+00, &
                     0.999943242212833D+00, &
                     0.999944464393028D+00, &
                     0.999945801963492D+00, &
                     0.999947238514747D+00, &
                     0.999948756506405D+00, &
                     0.999950337503071D+00, &
                     0.999951962419768D+00, &
                     0.999953611772821D+00, &
                     0.999955265932060D+00, &
                     0.999956905370104D+00, &
                     0.999958510904578D+00, &
                     0.999960063929208D+00, &
                     0.999961546630056D+00, &
                     0.999962942183517D+00, &
                     0.999964234933295D+00, &
                     0.999965410544314D+00, &
                     0.999966456132404D+00, &
                     0.999967360369669D+00, &
                     0.999968113566534D+00, &
                     0.999968707732561D+00, &
                     0.999969136618984D+00, &
                     0.999969395746500D+00, &
                     0.999969482421875D+00, &
                     0.999969482421875D+00, &
                     0.999969528246491D+00, &
                     0.999969665150031D+00, &
                     0.999969891429433D+00, &
                     0.999970204272348D+00, &
                     0.999970599795845D+00, &
                     0.999971073099757D+00, &
                     0.999971618333848D+00, &
                     0.999972228777745D+00, &
                     0.999972896932390D+00, &
                     0.999973614621581D+00, &
                     0.999974373101987D+00, &
                     0.999975163179889D+00, &
                     0.999975975332774D+00, &
                     0.999976799833804D+00, &
                     0.999977626877127D+00, &
                     0.999978446701969D+00, &
                     0.999979249713454D+00, &
                     0.999980026598159D+00, &
                     0.999980768432555D+00, &
                     0.999981466782641D+00, &
                     0.999982113793372D+00, &
                     0.999982702266828D+00, &
                     0.999983225728506D+00, &
                     0.999983678481607D+00, &
                     0.999984055649754D+00, &
                     0.999984353209074D+00, &
                     0.999984568011030D+00, &
                     0.999984697797664D+00, &
                     0.999984741210938D+00, &
                     0.999984741210938D+00, &
                     0.999984764093537D+00, &
                     0.999984832457465D+00, &
                     0.999984945454992D+00, &
                     0.999985101686110D+00, &
                     0.999985299217634D+00, &
                     0.999985535609544D+00, &
                     0.999985807948153D+00, &
                     0.999986112885601D+00, &
                     0.999986446685064D+00, &
                     0.999986805270981D+00, &
                     0.999987184283516D+00, &
                     0.999987579136394D+00, &
                     0.999987985077205D+00, &
                     0.999988397249208D+00, &
                     0.999988810753634D+00, &
                     0.999989220711480D+00, &
                     0.999989622323781D+00, &
                     0.999990010929386D+00, &
                     0.999990382059305D+00, &
                     0.999990731486800D+00, &
                     0.999991055272507D+00, &
                     0.999991349804054D+00, &
                     0.999991611829840D+00, &
                     0.999991838486901D+00, &
                     0.999992027323014D+00, &
                     0.999992176313496D+00, &
                     0.999992283873322D+00, &
                     0.999992348865351D+00, &
                     0.999992370605469D+00, &
                     0.999992370605469D+00, &
                     0.999992382033552D+00, &
                     0.999992416176436D+00, &
                     0.999992472611940D+00, &
                     0.999992550642790D+00, &
                     0.999992649306065D+00, &
                     0.999992767386215D+00, &
                     0.999992903431470D+00, &
                     0.999993055773372D+00, &
                     0.999993222549158D+00, &
                     0.999993401726630D+00, &
                     0.999993591131141D+00, &
                     0.999993788474282D+00, &
                     0.999993991383811D+00, &
                     0.999994197434362D+00, &
                     0.999994404178435D+00, &
                     0.999994609177182D+00, &
                     0.999994810030477D+00, &
                     0.999995004405795D+00, &
                     0.999995190065447D+00, &
                     0.999995364891736D+00, &
                     0.999995526909701D+00, &
                     0.999995674307166D+00, &
                     0.999995805451906D+00, &
                     0.999995918905908D+00, &
                     0.999996013436760D+00, &
                     0.999996088026390D+00, &
                     0.999996141877444D+00, &
                     0.999996174417676D+00, &
                     0.999996185302734D+00, &
                     0.999996185302734D+00, &
                     0.999996191010839D+00, &
                     0.999996208064724D+00, &
                     0.999996236254053D+00, &
                     0.999996275231412D+00, &
                     0.999996324516983D+00, &
                     0.999996383504990D+00, &
                     0.999996451471822D+00, &
                     0.999996527585705D+00, &
                     0.999996610917793D+00, &
                     0.999996700454495D+00, &
                     0.999996795110861D+00, &
                     0.999996893744817D+00, &
                     0.999996995172032D+00, &
                     0.999997098181186D+00, &
                     0.999997201549388D+00, &
                     0.999997304057523D+00, &
                     0.999997404505250D+00, &
                     0.999997501725444D+00, &
                     0.999997594597829D+00, &
                     0.999997682061617D+00, &
                     0.999997763126953D+00, &
                     0.999997836885045D+00, &
                     0.999997902516884D+00, &
                     0.999997959300507D+00, &
                     0.999998006616857D+00, &
                     0.999998043954301D+00, &
                     0.999998070911965D+00, &
                     0.999998087202060D+00, &
                     0.999998092651367D+00, &
                     0.999998092651367D+00, &
                     0.999998095502732D+00, &
                     0.999998104021728D+00, &
                     0.999998118103525D+00, &
                     0.999998137574969D+00, &
                     0.999998162196892D+00, &
                     0.999998191667310D+00, &
                     0.999998225625444D+00, &
                     0.999998263656517D+00, &
                     0.999998305297253D+00, &
                     0.999998350041992D+00, &
                     0.999998397349335D+00, &
                     0.999998446649213D+00, &
                     0.999998497350273D+00, &
                     0.999998548847476D+00, &
                     0.999998600529773D+00, &
                     0.999998651787750D+00, &
                     0.999998702021117D+00, &
                     0.999998750645923D+00, &
                     0.999998797101388D+00, &
                     0.999998840856239D+00, &
                     0.999998881414471D+00, &
                     0.999998918320464D+00, &
                     0.999998951163393D+00, &
                     0.999998979580934D+00, &
                     0.999999003262257D+00, &
                     0.999999021950363D+00, &
                     0.999999035443811D+00, &
                     0.999999043597947D+00, &
                     0.999999046325684D+00, &
                     0.999999046325684D+00, &
                     0.999999047750142D+00, &
                     0.999999052006021D+00, &
                     0.999999059041059D+00, &
                     0.999999068768930D+00, &
                     0.999999081070388D+00, &
                     0.999999095794849D+00, &
                     0.999999112762393D+00, &
                     0.999999131766134D+00, &
                     0.999999152574958D+00, &
                     0.999999174936552D+00, &
                     0.999999198580706D+00, &
                     0.999999223222828D+00, &
                     0.999999248567613D+00, &
                     0.999999274312827D+00, &
                     0.999999300153126D+00, &
                     0.999999325783873D+00, &
                     0.999999350904868D+00, &
                     0.999999375223963D+00, &
                     0.999999398460475D+00, &
                     0.999999420348368D+00, &
                     0.999999440639148D+00, &
                     0.999999459104446D+00, &
                     0.999999475538248D+00, &
                     0.999999489758774D+00, &
                     0.999999501610016D+00, &
                     0.999999510962931D+00, &
                     0.999999517716339D+00, &
                     0.999999521797563D+00, &
                     0.999999523162842D+00, &
                     0.999999523162842D+00, &
                     0.999999523874511D+00, &
                     0.999999526000794D+00, &
                     0.999999529515630D+00, &
                     0.999999534375972D+00, &
                     0.999999540522350D+00, &
                     0.999999547879659D+00, &
                     0.999999556358152D+00, &
                     0.999999565854619D+00, &
                     0.999999576253740D+00, &
                     0.999999587429596D+00, &
                     0.999999599247306D+00, &
                     0.999999611564777D+00, &
                     0.999999624234528D+00, &
                     0.999999637105573D+00, &
                     0.999999650025323D+00, &
                     0.999999662841491D+00, &
                     0.999999675403956D+00, &
                     0.999999687566564D+00, &
                     0.999999699188839D+00, &
                     0.999999710137581D+00, &
                     0.999999720288319D+00, &
                     0.999999729526610D+00, &
                     0.999999737749171D+00, &
                     0.999999744864831D+00, &
                     0.999999750795311D+00, &
                     0.999999755475838D+00, &
                     0.999999758855612D+00, &
                     0.999999760898134D+00, &
                     0.999999761581421D+00, &
                     0.999999761581421D+00, &
                     0.999999761936998D+00, &
                     0.999999762999378D+00, &
                     0.999999764755563D+00, &
                     0.999999767184082D+00, &
                     0.999999770255270D+00, &
                     0.999999773931662D+00, &
                     0.999999778168481D+00, &
                     0.999999782914228D+00, &
                     0.999999788111354D+00, &
                     0.999999793697007D+00, &
                     0.999999799603850D+00, &
                     0.999999805760931D+00, &
                     0.999999812094587D+00, &
                     0.999999818529386D+00, &
                     0.999999824989073D+00, &
                     0.999999831397519D+00, &
                     0.999999837679653D+00, &
                     0.999999843762362D+00, &
                     0.999999849575347D+00, &
                     0.999999855051924D+00, &
                     0.999999860129755D+00, &
                     0.999999864751499D+00, &
                     0.999999868865388D+00, &
                     0.999999872425704D+00, &
                     0.999999875393184D+00, &
                     0.999999877735324D+00, &
                     0.999999879426627D+00, &
                     0.999999880448768D+00, &
                     0.999999880790710D+00, &
                     0.999999880790710D+00, &
                     0.999999880968380D+00, &
                     0.999999881499219D+00, &
                     0.999999882376743D+00, &
                     0.999999883590240D+00, &
                     0.999999885124911D+00, &
                     0.999999886962062D+00, &
                     0.999999889079351D+00, &
                     0.999999891451077D+00, &
                     0.999999894048515D+00, &
                     0.999999896840290D+00, &
                     0.999999899792782D+00, &
                     0.999999902870557D+00, &
                     0.999999906036821D+00, &
                     0.999999909253884D+00, &
                     0.999999912483639D+00, &
                     0.999999915688027D+00, &
                     0.999999918829509D+00, &
                     0.999999921871511D+00, &
                     0.999999924778856D+00, &
                     0.999999927518164D+00, &
                     0.999999930058216D+00, &
                     0.999999932370289D+00, &
                     0.999999934428439D+00, &
                     0.999999936209747D+00, &
                     0.999999937694523D+00, &
                     0.999999938866461D+00, &
                     0.999999939712768D+00, &
                     0.999999940224246D+00, &
                     0.999999940395355D+00, &
                     0.999999940395355D+00, &
                     0.999999940484135D+00, &
                     0.999999940749392D+00, &
                     0.999999941187891D+00, &
                     0.999999941794286D+00, &
                     0.999999942561194D+00, &
                     0.999999943479286D+00, &
                     0.999999944537412D+00, &
                     0.999999945722743D+00, &
                     0.999999947020941D+00, &
                     0.999999948416342D+00, &
                     0.999999949892156D+00, &
                     0.999999951430689D+00, &
                     0.999999953013558D+00, &
                     0.999999954621934D+00, &
                     0.999999956236769D+00, &
                     0.999999957839039D+00, &
                     0.999999959409972D+00, &
                     0.999999960931272D+00, &
                     0.999999962385339D+00, &
                     0.999999963755465D+00, &
                     0.999999965026018D+00, &
                     0.999999966182611D+00, &
                     0.999999967212245D+00, &
                     0.999999968103433D+00, &
                     0.999999968846301D+00, &
                     0.999999969432673D+00, &
                     0.999999969856130D+00, &
                     0.999999970112059D+00, &
                     0.999999970197678D+00, &
                     0.999999970197678D+00, &
                     0.999999970242042D+00, &
                     0.999999970374595D+00, &
                     0.999999970593722D+00, &
                     0.999999970896756D+00, &
                     0.999999971280012D+00, &
                     0.999999971738833D+00, &
                     0.999999972267655D+00, &
                     0.999999972860073D+00, &
                     0.999999973508930D+00, &
                     0.999999974206404D+00, &
                     0.999999974944111D+00, &
                     0.999999975713212D+00, &
                     0.999999976504525D+00, &
                     0.999999977308639D+00, &
                     0.999999978116037D+00, &
                     0.999999978917207D+00, &
                     0.999999979702762D+00, &
                     0.999999980463551D+00, &
                     0.999999981190768D+00, &
                     0.999999981876050D+00, &
                     0.999999982511572D+00, &
                     0.999999983090127D+00, &
                     0.999999983605204D+00, &
                     0.999999984051046D+00, &
                     0.999999984422704D+00, &
                     0.999999984716077D+00, &
                     0.999999984927947D+00, &
                     0.999999985055999D+00, &
                     0.999999985098839D+00, &
                     0.999999985098839D+00, &
                     0.999999985121009D+00, &
                     0.999999985187250D+00, &
                     0.999999985296757D+00, &
                     0.999999985448198D+00, &
                     0.999999985639733D+00, &
                     0.999999985869039D+00, &
                     0.999999986133338D+00, &
                     0.999999986429432D+00, &
                     0.999999986753748D+00, &
                     0.999999987102379D+00, &
                     0.999999987471139D+00, &
                     0.999999987855613D+00, &
                     0.999999988251212D+00, &
                     0.999999988653235D+00, &
                     0.999999989056925D+00, &
                     0.999999989457526D+00, &
                     0.999999989850345D+00, &
                     0.999999990230804D+00, &
                     0.999999990594498D+00, &
                     0.999999990937241D+00, &
                     0.999999991255115D+00, &
                     0.999999991544514D+00, &
                     0.999999991802173D+00, &
                     0.999999992025210D+00, &
                     0.999999992211143D+00, &
                     0.999999992357918D+00, &
                     0.999999992463919D+00, &
                     0.999999992527986D+00, &
                     0.999999992549419D+00, &
                     0.999999992549419D+00, &
                     0.999999992560499D+00, &
                     0.999999992593603D+00, &
                     0.999999992648330D+00, &
                     0.999999992724015D+00, &
                     0.999999992819739D+00, &
                     0.999999992934344D+00, &
                     0.999999993066441D+00, &
                     0.999999993214434D+00, &
                     0.999999993376539D+00, &
                     0.999999993550805D+00, &
                     0.999999993735142D+00, &
                     0.999999993927342D+00, &
                     0.999999994125115D+00, &
                     0.999999994326111D+00, &
                     0.999999994527951D+00, &
                     0.999999994728259D+00, &
                     0.999999994924688D+00, &
                     0.999999995114948D+00, &
                     0.999999995296834D+00, &
                     0.999999995468254D+00, &
                     0.999999995627244D+00, &
                     0.999999995772000D+00, &
                     0.999999995900886D+00, &
                     0.999999996012459D+00, &
                     0.999999996105474D+00, &
                     0.999999996178902D+00, &
                     0.999999996231934D+00, &
                     0.999999996263986D+00, &
                     0.999999996274710D+00, &
                     0.999999996274710D+00, &
                     0.999999996280247D+00, &
                     0.999999996296791D+00, &
                     0.999999996324142D+00, &
                     0.999999996361968D+00, &
                     0.999999996409810D+00, &
                     0.999999996467089D+00, &
                     0.999999996533113D+00, &
                     0.999999996607085D+00, &
                     0.999999996688113D+00, &
                     0.999999996775223D+00, &
                     0.999999996867371D+00, &
                     0.999999996963454D+00, &
                     0.999999997062328D+00, &
                     0.999999997162818D+00, &
                     0.999999997263736D+00, &
                     0.999999997363894D+00, &
                     0.999999997462117D+00, &
                     0.999999997557261D+00, &
                     0.999999997648223D+00, &
                     0.999999997733955D+00, &
                     0.999999997813475D+00, &
                     0.999999997885880D+00, &
                     0.999999997950349D+00, &
                     0.999999998006161D+00, &
                     0.999999998052691D+00, &
                     0.999999998089425D+00, &
                     0.999999998115955D+00, &
                     0.999999998131990D+00, &
                     0.999999998137355D+00, &
                     0.999999998137355D+00, &
                     0.999999998140122D+00, &
                     0.999999998148391D+00, &
                     0.999999998162061D+00, &
                     0.999999998180966D+00, &
                     0.999999998204877D+00, &
                     0.999999998233506D+00, &
                     0.999999998266507D+00, &
                     0.999999998303481D+00, &
                     0.999999998343983D+00, &
                     0.999999998387527D+00, &
                     0.999999998433591D+00, &
                     0.999999998481625D+00, &
                     0.999999998531056D+00, &
                     0.999999998581298D+00, &
                     0.999999998631756D+00, &
                     0.999999998681836D+00, &
                     0.999999998730952D+00, &
                     0.999999998778531D+00, &
                     0.999999998824020D+00, &
                     0.999999998866897D+00, &
                     0.999999998906669D+00, &
                     0.999999998942883D+00, &
                     0.999999998975131D+00, &
                     0.999999999003048D+00, &
                     0.999999999026324D+00, &
                     0.999999999044700D+00, &
                     0.999999999057972D+00, &
                     0.999999999065994D+00, &
                     0.999999999068677D+00, &
                     0.999999999068677D+00, &
                     0.999999999070061D+00, &
                     0.999999999074193D+00, &
                     0.999999999081025D+00, &
                     0.999999999090474D+00, &
                     0.999999999102426D+00, &
                     0.999999999116735D+00, &
                     0.999999999133230D+00, &
                     0.999999999151711D+00, &
                     0.999999999171957D+00, &
                     0.999999999193724D+00, &
                     0.999999999216752D+00, &
                     0.999999999240765D+00, &
                     0.999999999265478D+00, &
                     0.999999999290597D+00, &
                     0.999999999315825D+00, &
                     0.999999999340866D+00, &
                     0.999999999365426D+00, &
                     0.999999999389218D+00, &
                     0.999999999411967D+00, &
                     0.999999999433411D+00, &
                     0.999999999453302D+00, &
                     0.999999999471415D+00, &
                     0.999999999487545D+00, &
                     0.999999999501509D+00, &
                     0.999999999513152D+00, &
                     0.999999999522344D+00, &
                     0.999999999528983D+00, &
                     0.999999999532996D+00, &
                     0.999999999534339D+00, &
                     0.999999999534339D+00, &
                     0.999999999535030D+00, &
                     0.999999999537096D+00, &
                     0.999999999540510D+00, &
                     0.999999999545233D+00, &
                     0.999999999551207D+00, &
                     0.999999999558359D+00, &
                     0.999999999566604D+00, &
                     0.999999999575842D+00, &
                     0.999999999585962D+00, &
                     0.999999999596843D+00, &
                     0.999999999608355D+00, &
                     0.999999999620360D+00, &
                     0.999999999632715D+00, &
                     0.999999999645274D+00, &
                     0.999999999657888D+00, &
                     0.999999999670409D+00, &
                     0.999999999682689D+00, &
                     0.999999999694587D+00, &
                     0.999999999705964D+00, &
                     0.999999999716687D+00, &
                     0.999999999726636D+00, &
                     0.999999999735695D+00, &
                     0.999999999743762D+00, &
                     0.999999999750747D+00, &
                     0.999999999756571D+00, &
                     0.999999999761169D+00, &
                     0.999999999764490D+00, &
                     0.999999999766498D+00, &
                     0.999999999767169D+00, &
                     0.999999999767169D+00, &
                     0.999999999767515D+00, &
                     0.999999999768547D+00, &
                     0.999999999770254D+00, &
                     0.999999999772615D+00, &
                     0.999999999775600D+00, &
                     0.999999999779175D+00, &
                     0.999999999783297D+00, &
                     0.999999999787915D+00, &
                     0.999999999792973D+00, &
                     0.999999999798413D+00, &
                     0.999999999804168D+00, &
                     0.999999999810169D+00, &
                     0.999999999816346D+00, &
                     0.999999999822625D+00, &
                     0.999999999828932D+00, &
                     0.999999999835193D+00, &
                     0.999999999841334D+00, &
                     0.999999999847283D+00, &
                     0.999999999852972D+00, &
                     0.999999999858335D+00, &
                     0.999999999863311D+00, &
                     0.999999999867842D+00, &
                     0.999999999871877D+00, &
                     0.999999999875370D+00, &
                     0.999999999878283D+00, &
                     0.999999999880583D+00, &
                     0.999999999882245D+00, &
                     0.999999999883249D+00, &
                     0.999999999883585D+00, &
                     0.999999999883585D+00, &
                     0.999999999883757D+00, &
                     0.999999999884273D+00, &
                     0.999999999885126D+00, &
                     0.999999999886306D+00, &
                     0.999999999887799D+00, &
                     0.999999999889586D+00, &
                     0.999999999891646D+00, &
                     0.999999999893954D+00, &
                     0.999999999896483D+00, &
                     0.999999999899202D+00, &
                     0.999999999902079D+00, &
                     0.999999999905080D+00, &
                     0.999999999908168D+00, &
                     0.999999999911307D+00, &
                     0.999999999914461D+00, &
                     0.999999999917591D+00, &
                     0.999999999920662D+00, &
                     0.999999999923637D+00, &
                     0.999999999926482D+00, &
                     0.999999999929164D+00, &
                     0.999999999931652D+00, &
                     0.999999999933918D+00, &
                     0.999999999935936D+00, &
                     0.999999999937684D+00, &
                     0.999999999939141D+00, &
                     0.999999999940291D+00, &
                     0.999999999941122D+00, &
                     0.999999999941624D+00, &
                     0.999999999941792D+00, &
                     0.999999999941792D+00, &
                     0.999999999941879D+00, &
                     0.999999999942137D+00, &
                     0.999999999942563D+00, &
                     0.999999999943153D+00, &
                     0.999999999943899D+00, &
                     0.999999999944792D+00, &
                     0.999999999945822D+00, &
                     0.999999999946976D+00, &
                     0.999999999948240D+00, &
                     0.999999999949599D+00, &
                     0.999999999951037D+00, &
                     0.999999999952537D+00, &
                     0.999999999954081D+00, &
                     0.999999999955651D+00, &
                     0.999999999957228D+00, &
                     0.999999999958793D+00, &
                     0.999999999960328D+00, &
                     0.999999999961816D+00, &
                     0.999999999963239D+00, &
                     0.999999999964580D+00, &
                     0.999999999965824D+00, &
                     0.999999999966958D+00, &
                     0.999999999967967D+00, &
                     0.999999999968841D+00, &
                     0.999999999969570D+00, &
                     0.999999999970145D+00, &
                     0.999999999970561D+00, &
                     0.999999999970812D+00, &
                     0.999999999970896D+00, &
                     0.999999999970896D+00, &
                     0.999999999970985D+00, &
                     0.999999999971249D+00, &
                     0.999999999971686D+00, &
                     0.999999999972291D+00, &
                     0.999999999973055D+00, &
                     0.999999999973970D+00, &
                     0.999999999975023D+00, &
                     0.999999999976203D+00, &
                     0.999999999977493D+00, &
                     0.999999999978879D+00, &
                     0.999999999980344D+00, &
                     0.999999999981869D+00, &
                     0.999999999983435D+00, &
                     0.999999999985024D+00, &
                     0.999999999986616D+00, &
                     0.999999999988192D+00, &
                     0.999999999989732D+00, &
                     0.999999999991219D+00, &
                     0.999999999992634D+00, &
                     0.999999999993962D+00, &
                     0.999999999995186D+00, &
                     0.999999999996292D+00, &
                     0.999999999997269D+00, &
                     0.999999999998106D+00, &
                     0.999999999998795D+00, &
                     0.999999999999330D+00, &
                     0.999999999999709D+00, &
                     0.999999999999930D+00, &
                     0.100000000000000D+01, &
                     0.100000000000000D+01, &
                     0.100000000000005D+01, &
                     0.100000000000022D+01, &
                     0.100000000000051D+01, &
                     0.100000000000092D+01, &
                     0.100000000000145D+01, &
                     0.100000000000209D+01, &
                     0.100000000000283D+01, &
                     0.100000000000368D+01, &
                     0.100000000000461D+01, &
                     0.100000000000562D+01, &
                     0.100000000000670D+01, &
                     0.100000000000784D+01, &
                     0.100000000000901D+01, &
                     0.100000000001021D+01, &
                     0.100000000001143D+01, &
                     0.100000000001264D+01, &
                     0.100000000001384D+01, &
                     0.100000000001500D+01, &
                     0.100000000001612D+01, &
                     0.100000000001717D+01, &
                     0.100000000001816D+01, &
                     0.100000000001906D+01, &
                     0.100000000001986D+01, &
                     0.100000000002056D+01, &
                     0.100000000002114D+01, &
                     0.100000000002160D+01, &
                     0.100000000002193D+01, &
                     0.100000000002214D+01, &
                     0.100000000002220D+01, &
                     0.100000000002220D+01, &
                     0.100000000002227D+01, &
                     0.100000000002246D+01, &
                     0.100000000002278D+01, &
                     0.100000000002322D+01, &
                     0.100000000002377D+01, &
                     0.100000000002444D+01, &
                     0.100000000002521D+01, &
                     0.100000000002607D+01, &
                     0.100000000002702D+01, &
                     0.100000000002805D+01, &
                     0.100000000002913D+01, &
                     0.100000000003027D+01, &
                     0.100000000003144D+01, &
                     0.100000000003263D+01, &
                     0.100000000003384D+01, &
                     0.100000000003503D+01, &
                     0.100000000003621D+01, &
                     0.100000000003736D+01, &
                     0.100000000003845D+01, &
                     0.100000000003949D+01, &
                     0.100000000004045D+01, &
                     0.100000000004133D+01, &
                     0.100000000004212D+01, &
                     0.100000000004280D+01, &
                     0.100000000004337D+01, &
                     0.100000000004382D+01, &
                     0.100000000004415D+01, &
                     0.100000000004434D+01, &
                     0.100000000004441D+01, &
                     0.100000000004441D+01, &
                     0.100000000004454D+01, &
                     0.100000000004492D+01, &
                     0.100000000004555D+01, &
                     0.100000000004643D+01, &
                     0.100000000004754D+01, &
                     0.100000000004888D+01, &
                     0.100000000005042D+01, &
                     0.100000000005215D+01, &
                     0.100000000005404D+01, &
                     0.100000000005609D+01, &
                     0.100000000005826D+01, &
                     0.100000000006053D+01, &
                     0.100000000006287D+01, &
                     0.100000000006526D+01, &
                     0.100000000006767D+01, &
                     0.100000000007006D+01, &
                     0.100000000007242D+01, &
                     0.100000000007471D+01, &
                     0.100000000007690D+01, &
                     0.100000000007898D+01, &
                     0.100000000008091D+01, &
                     0.100000000008267D+01, &
                     0.100000000008424D+01, &
                     0.100000000008560D+01, &
                     0.100000000008674D+01, &
                     0.100000000008764D+01, &
                     0.100000000008829D+01, &
                     0.100000000008869D+01, &
                     0.100000000008882D+01, &
                     0.100000000008882D+01, &
                     0.100000000008907D+01, &
                     0.100000000008984D+01, &
                     0.100000000009111D+01, &
                     0.100000000009286D+01, &
                     0.100000000009508D+01, &
                     0.100000000009775D+01, &
                     0.100000000010083D+01, &
                     0.100000000010429D+01, &
                     0.100000000010808D+01, &
                     0.100000000011217D+01, &
                     0.100000000011652D+01, &
                     0.100000000012105D+01, &
                     0.100000000012574D+01, &
                     0.100000000013052D+01, &
                     0.100000000013533D+01, &
                     0.100000000014012D+01, &
                     0.100000000014483D+01, &
                     0.100000000014941D+01, &
                     0.100000000015380D+01, &
                     0.100000000015795D+01, &
                     0.100000000016181D+01, &
                     0.100000000016533D+01, &
                     0.100000000016848D+01, &
                     0.100000000017120D+01, &
                     0.100000000017348D+01, &
                     0.100000000017528D+01, &
                     0.100000000017658D+01, &
                     0.100000000017737D+01, &
                     0.100000000017763D+01, &
                     0.100000000017763D+01, &
                     0.100000000017815D+01, &
                     0.100000000017968D+01, &
                     0.100000000018221D+01, &
                     0.100000000018572D+01, &
                     0.100000000019016D+01, &
                     0.100000000019549D+01, &
                     0.100000000020165D+01, &
                     0.100000000020856D+01, &
                     0.100000000021615D+01, &
                     0.100000000022434D+01, &
                     0.100000000023302D+01, &
                     0.100000000024209D+01, &
                     0.100000000025146D+01, &
                     0.100000000026102D+01, &
                     0.100000000027064D+01, &
                     0.100000000028022D+01, &
                     0.100000000028965D+01, &
                     0.100000000029880D+01, &
                     0.100000000030759D+01, &
                     0.100000000031589D+01, &
                     0.100000000032361D+01, &
                     0.100000000033065D+01, &
                     0.100000000033694D+01, &
                     0.100000000034240D+01, &
                     0.100000000034696D+01, &
                     0.100000000035056D+01, &
                     0.100000000035316D+01, &
                     0.100000000035474D+01, &
                     0.100000000035527D+01, &
                     0.100000000035527D+01, &
                     0.100000000035629D+01, &
                     0.100000000035935D+01, &
                     0.100000000036442D+01, &
                     0.100000000037143D+01, &
                     0.100000000038032D+01, &
                     0.100000000039097D+01, &
                     0.100000000040328D+01, &
                     0.100000000041710D+01, &
                     0.100000000043228D+01, &
                     0.100000000044864D+01, &
                     0.100000000046600D+01, &
                     0.100000000048415D+01, &
                     0.100000000050289D+01, &
                     0.100000000052199D+01, &
                     0.100000000054124D+01, &
                     0.100000000056041D+01, &
                     0.100000000057926D+01, &
                     0.100000000059758D+01, &
                     0.100000000061514D+01, &
                     0.100000000063175D+01, &
                     0.100000000064719D+01, &
                     0.100000000066129D+01, &
                     0.100000000067388D+01, &
                     0.100000000068479D+01, &
                     0.100000000069391D+01, &
                     0.100000000070112D+01, &
                     0.100000000070633D+01, &
                     0.100000000070948D+01, &
                     0.100000000071054D+01, &
                     0.100000000071054D+01, &
                     0.100000000071258D+01, &
                     0.100000000071870D+01, &
                     0.100000000072883D+01, &
                     0.100000000074285D+01, &
                     0.100000000076062D+01, &
                     0.100000000078192D+01, &
                     0.100000000080652D+01, &
                     0.100000000083416D+01, &
                     0.100000000086450D+01, &
                     0.100000000089722D+01, &
                     0.100000000093193D+01, &
                     0.100000000096823D+01, &
                     0.100000000100570D+01, &
                     0.100000000104391D+01, &
                     0.100000000108240D+01, &
                     0.100000000112073D+01, &
                     0.100000000115845D+01, &
                     0.100000000119509D+01, &
                     0.100000000123023D+01, &
                     0.100000000126345D+01, &
                     0.100000000129435D+01, &
                     0.100000000132255D+01, &
                     0.100000000134773D+01, &
                     0.100000000136956D+01, &
                     0.100000000138780D+01, &
                     0.100000000140223D+01, &
                     0.100000000141265D+01, &
                     0.100000000141896D+01, &
                     0.100000000142107D+01, &
                     0.100000000142107D+01, &
                     0.100000000142517D+01, &
                     0.100000000143740D+01, &
                     0.100000000145765D+01, &
                     0.100000000148568D+01, &
                     0.100000000152118D+01, &
                     0.100000000156377D+01, &
                     0.100000000161296D+01, &
                     0.100000000166821D+01, &
                     0.100000000172889D+01, &
                     0.100000000179430D+01, &
                     0.100000000186370D+01, &
                     0.100000000193630D+01, &
                     0.100000000201124D+01, &
                     0.100000000208765D+01, &
                     0.100000000216464D+01, &
                     0.100000000224131D+01, &
                     0.100000000231674D+01, &
                     0.100000000239004D+01, &
                     0.100000000246034D+01, &
                     0.100000000252679D+01, &
                     0.100000000258860D+01, &
                     0.100000000264503D+01, &
                     0.100000000269539D+01, &
                     0.100000000273909D+01, &
                     0.100000000277558D+01, &
                     0.100000000280444D+01, &
                     0.100000000282530D+01, &
                     0.100000000283792D+01, &
                     0.100000000284214D+01, &
                     0.100000000284214D+01, &
                     0.100000000285033D+01, &
                     0.100000000287479D+01, &
                     0.100000000291525D+01, &
                     0.100000000297129D+01, &
                     0.100000000304226D+01, &
                     0.100000000312740D+01, &
                     0.100000000322575D+01, &
                     0.100000000333620D+01, &
                     0.100000000345752D+01, &
                     0.100000000358831D+01, &
                     0.100000000372709D+01, &
                     0.100000000387225D+01, &
                     0.100000000402212D+01, &
                     0.100000000417495D+01, &
                     0.100000000432893D+01, &
                     0.100000000448228D+01, &
                     0.100000000463316D+01, &
                     0.100000000477978D+01, &
                     0.100000000492041D+01, &
                     0.100000000505334D+01, &
                     0.100000000517700D+01, &
                     0.100000000528990D+01, &
                     0.100000000539066D+01, &
                     0.100000000547809D+01, &
                     0.100000000555110D+01, &
                     0.100000000560884D+01, &
                     0.100000000565058D+01, &
                     0.100000000567583D+01, &
                     0.100000000568429D+01, &
                     0.100000000568429D+01, &
                     0.100000000570064D+01, &
                     0.100000000574953D+01, &
                     0.100000000583042D+01, &
                     0.100000000594243D+01, &
                     0.100000000608431D+01, &
                     0.100000000625450D+01, &
                     0.100000000645112D+01, &
                     0.100000000667194D+01, &
                     0.100000000691449D+01, &
                     0.100000000717601D+01, &
                     0.100000000745351D+01, &
                     0.100000000774379D+01, &
                     0.100000000804350D+01, &
                     0.100000000834913D+01, &
                     0.100000000865712D+01, &
                     0.100000000896382D+01, &
                     0.100000000926562D+01, &
                     0.100000000955893D+01, &
                     0.100000000984024D+01, &
                     0.100000001010619D+01, &
                     0.100000001035359D+01, &
                     0.100000001057946D+01, &
                     0.100000001078107D+01, &
                     0.100000001095598D+01, &
                     0.100000001110208D+01, &
                     0.100000001121760D+01, &
                     0.100000001130113D+01, &
                     0.100000001135166D+01, &
                     0.100000001136857D+01, &
                     0.100000001136857D+01, &
                     0.100000001140126D+01, &
                     0.100000001149898D+01, &
                     0.100000001166066D+01, &
                     0.100000001188454D+01, &
                     0.100000001216815D+01, &
                     0.100000001250837D+01, &
                     0.100000001290141D+01, &
                     0.100000001334289D+01, &
                     0.100000001382782D+01, &
                     0.100000001435070D+01, &
                     0.100000001490557D+01, &
                     0.100000001548604D+01, &
                     0.100000001608539D+01, &
                     0.100000001669664D+01, &
                     0.100000001731262D+01, &
                     0.100000001792609D+01, &
                     0.100000001852977D+01, &
                     0.100000001911649D+01, &
                     0.100000001967926D+01, &
                     0.100000002021132D+01, &
                     0.100000002070628D+01, &
                     0.100000002115819D+01, &
                     0.100000002156157D+01, &
                     0.100000002191156D+01, &
                     0.100000002220390D+01, &
                     0.100000002243505D+01, &
                     0.100000002260219D+01, &
                     0.100000002270330D+01, &
                     0.100000002273714D+01, &
                     0.100000002273714D+01, &
                     0.100000002280248D+01, &
                     0.100000002299779D+01, &
                     0.100000002332094D+01, &
                     0.100000002376842D+01, &
                     0.100000002433531D+01, &
                     0.100000002501538D+01, &
                     0.100000002580108D+01, &
                     0.100000002668364D+01, &
                     0.100000002765313D+01, &
                     0.100000002869858D+01, &
                     0.100000002980805D+01, &
                     0.100000003096878D+01, &
                     0.100000003216735D+01, &
                     0.100000003338979D+01, &
                     0.100000003462179D+01, &
                     0.100000003584883D+01, &
                     0.100000003705637D+01, &
                     0.100000003823007D+01, &
                     0.100000003935589D+01, &
                     0.100000004042035D+01, &
                     0.100000004141064D+01, &
                     0.100000004231481D+01, &
                     0.100000004312195D+01, &
                     0.100000004382225D+01, &
                     0.100000004440723D+01, &
                     0.100000004486977D+01, &
                     0.100000004520424D+01, &
                     0.100000004540656D+01, &
                     0.100000004547428D+01, &
                     0.100000004547428D+01, &
                     0.100000004560486D+01, &
                     0.100000004599520D+01, &
                     0.100000004664106D+01, &
                     0.100000004753542D+01, &
                     0.100000004866849D+01, &
                     0.100000005002784D+01, &
                     0.100000005159841D+01, &
                     0.100000005336270D+01, &
                     0.100000005530091D+01, &
                     0.100000005739110D+01, &
                     0.100000005960945D+01, &
                     0.100000006193047D+01, &
                     0.100000006432733D+01, &
                     0.100000006677210D+01, &
                     0.100000006923616D+01, &
                     0.100000007169046D+01, &
                     0.100000007410595D+01, &
                     0.100000007645386D+01, &
                     0.100000007870615D+01, &
                     0.100000008083579D+01, &
                     0.100000008281714D+01, &
                     0.100000008462628D+01, &
                     0.100000008624131D+01, &
                     0.100000008764263D+01, &
                     0.100000008881322D+01, &
                     0.100000008973882D+01, &
                     0.100000009040815D+01, &
                     0.100000009081304D+01, &
                     0.100000009094856D+01, &
                     0.100000009094856D+01, &
                     0.100000009120952D+01, &
                     0.100000009198960D+01, &
                     0.100000009328035D+01, &
                     0.100000009506777D+01, &
                     0.100000009733239D+01, &
                     0.100000010004937D+01, &
                     0.100000010318873D+01, &
                     0.100000010671554D+01, &
                     0.100000011059028D+01, &
                     0.100000011476915D+01, &
                     0.100000011920459D+01, &
                     0.100000012384568D+01, &
                     0.100000012863877D+01, &
                     0.100000013352808D+01, &
                     0.100000013845633D+01, &
                     0.100000014336544D+01, &
                     0.100000014819724D+01, &
                     0.100000015289420D+01, &
                     0.100000015740015D+01, &
                     0.100000016166098D+01, &
                     0.100000016562535D+01, &
                     0.100000016924532D+01, &
                     0.100000017247704D+01, &
                     0.100000017528124D+01, &
                     0.100000017762378D+01, &
                     0.100000017947610D+01, &
                     0.100000018081560D+01, &
                     0.100000018162591D+01, &
                     0.100000018189712D+01, &
                     0.100000018189712D+01, &
                     0.100000018241861D+01, &
                     0.100000018397745D+01, &
                     0.100000018655686D+01, &
                     0.100000019012894D+01, &
                     0.100000019465487D+01, &
                     0.100000020008516D+01, &
                     0.100000020636002D+01, &
                     0.100000021340980D+01, &
                     0.100000022115564D+01, &
                     0.100000022951014D+01, &
                     0.100000023837828D+01, &
                     0.100000024765839D+01, &
                     0.100000025724325D+01, &
                     0.100000026702134D+01, &
                     0.100000027687811D+01, &
                     0.100000028669740D+01, &
                     0.100000029636281D+01, &
                     0.100000030575918D+01, &
                     0.100000031477405D+01, &
                     0.100000032329906D+01, &
                     0.100000033123139D+01, &
                     0.100000033847501D+01, &
                     0.100000034494203D+01, &
                     0.100000035055376D+01, &
                     0.100000035524179D+01, &
                     0.100000035894889D+01, &
                     0.100000036162971D+01, &
                     0.100000036325145D+01, &
                     0.100000036379424D+01, &
                     0.100000036379424D+01, &
                     0.100000036483626D+01, &
                     0.100000036795114D+01, &
                     0.100000037310543D+01, &
                     0.100000038024358D+01, &
                     0.100000038928828D+01, &
                     0.100000040014092D+01, &
                     0.100000041268229D+01, &
                     0.100000042677355D+01, &
                     0.100000044225735D+01, &
                     0.100000045895930D+01, &
                     0.100000047668967D+01, &
                     0.100000049524538D+01, &
                     0.100000051441222D+01, &
                     0.100000053396723D+01, &
                     0.100000055368137D+01, &
                     0.100000057332225D+01, &
                     0.100000059265696D+01, &
                     0.100000061145500D+01, &
                     0.100000062949117D+01, &
                     0.100000064654848D+01, &
                     0.100000066242091D+01, &
                     0.100000067691612D+01, &
                     0.100000068985792D+01, &
                     0.100000070108862D+01, &
                     0.100000071047109D+01, &
                     0.100000071789056D+01, &
                     0.100000072325615D+01, &
                     0.100000072650207D+01, &
                     0.100000072758849D+01, &
                     0.100000072758849D+01, &
                     0.100000072967045D+01, &
                     0.100000073589410D+01, &
                     0.100000074619286D+01, &
                     0.100000076045616D+01, &
                     0.100000077853002D+01, &
                     0.100000080021803D+01, &
                     0.100000082528268D+01, &
                     0.100000085344716D+01, &
                     0.100000088439764D+01, &
                     0.100000091778619D+01, &
                     0.100000095323406D+01, &
                     0.100000099033569D+01, &
                     0.100000102866305D+01, &
                     0.100000106777052D+01, &
                     0.100000110720008D+01, &
                     0.100000114648683D+01, &
                     0.100000118516469D+01, &
                     0.100000122277224D+01, &
                     0.100000125885855D+01, &
                     0.100000129298897D+01, &
                     0.100000132475077D+01, &
                     0.100000135375847D+01, &
                     0.100000137965897D+01, &
                     0.100000140213613D+01, &
                     0.100000142091500D+01, &
                     0.100000143576545D+01, &
                     0.100000144650521D+01, &
                     0.100000145300234D+01, &
                     0.100000145517697D+01, &
                     0.100000145517697D+01, &
                     0.100000145933639D+01, &
                     0.100000147177041D+01, &
                     0.100000149234659D+01, &
                     0.100000152084487D+01, &
                     0.100000155695879D+01, &
                     0.100000160029725D+01, &
                     0.100000165038712D+01, &
                     0.100000170667675D+01, &
                     0.100000176854044D+01, &
                     0.100000183528406D+01, &
                     0.100000190615174D+01, &
                     0.100000198033359D+01, &
                     0.100000205697451D+01, &
                     0.100000213518381D+01, &
                     0.100000221404562D+01, &
                     0.100000229262993D+01, &
                     0.100000237000401D+01, &
                     0.100000244524409D+01, &
                     0.100000251744713D+01, &
                     0.100000258574242D+01, &
                     0.100000264930291D+01, &
                     0.100000270735602D+01, &
                     0.100000275919386D+01, &
                     0.100000280418255D+01, &
                     0.100000284177068D+01, &
                     0.100000287149667D+01, &
                     0.100000289299493D+01, &
                     0.100000290600077D+01, &
                     0.100000291035394D+01, &
                     0.100000291035394D+01, &
                     0.100000291866295D+01, &
                     0.100000294350204D+01, &
                     0.100000298460782D+01, &
                     0.100000304154262D+01, &
                     0.100000311369670D+01, &
                     0.100000320029163D+01, &
                     0.100000330038530D+01, &
                     0.100000341287868D+01, &
                     0.100000353652459D+01, &
                     0.100000366993867D+01, &
                     0.100000381161260D+01, &
                     0.100000395992942D+01, &
                     0.100000411318096D+01, &
                     0.100000426958709D+01, &
                     0.100000442731645D+01, &
                     0.100000458450857D+01, &
                     0.100000473929671D+01, &
                     0.100000488983136D+01, &
                     0.100000503430385D+01, &
                     0.100000517096969D+01, &
                     0.100000529817132D+01, &
                     0.100000541436000D+01, &
                     0.100000551811624D+01, &
                     0.100000560816878D+01, &
                     0.100000568341152D+01, &
                     0.100000574291844D+01, &
                     0.100000578595597D+01, &
                     0.100000581199299D+01, &
                     0.100000582070788D+01, &
                     0.100000582070788D+01, &
                     0.100000583730442D+01, &
                     0.100000588691932D+01, &
                     0.100000596902908D+01, &
                     0.100000608276366D+01, &
                     0.100000622691052D+01, &
                     0.100000639992104D+01, &
                     0.100000659992006D+01, &
                     0.100000682471885D+01, &
                     0.100000707183224D+01, &
                     0.100000733850008D+01, &
                     0.100000762171323D+01, &
                     0.100000791824394D+01, &
                     0.100000822468038D+01, &
                     0.100000853746500D+01, &
                     0.100000885293599D+01, &
                     0.100000916737137D+01, &
                     0.100000947703496D+01, &
                     0.100000977822344D+01, &
                     0.100001006731376D+01, &
                     0.100001034081018D+01, &
                     0.100001059539012D+01, &
                     0.100001082794812D+01, &
                     0.100001103563721D+01, &
                     0.100001121590706D+01, &
                     0.100001136653833D+01, &
                     0.100001148567260D+01, &
                     0.100001157183761D+01, &
                     0.100001162396719D+01, &
                     0.100001164141577D+01, &
                     0.100001164141577D+01, &
                     0.100001167456187D+01, &
                     0.100001177365317D+01, &
                     0.100001193764990D+01, &
                     0.100001216482351D+01, &
                     0.100001245276408D+01, &
                     0.100001279839240D+01, &
                     0.100001319797786D+01, &
                     0.100001364716349D+01, &
                     0.100001414099907D+01, &
                     0.100001467398305D+01, &
                     0.100001524011364D+01, &
                     0.100001583294884D+01, &
                     0.100001644567498D+01, &
                     0.100001707118295D+01, &
                     0.100001770215111D+01, &
                     0.100001833113340D+01, &
                     0.100001895065151D+01, &
                     0.100001955328935D+01, &
                     0.100002013178842D+01, &
                     0.100002067914244D+01, &
                     0.100002118868976D+01, &
                     0.100002165420204D+01, &
                     0.100002206996770D+01, &
                     0.100002243086905D+01, &
                     0.100002273245156D+01, &
                     0.100002297098455D+01, &
                     0.100002314351201D+01, &
                     0.100002324789313D+01, &
                     0.100002328283153D+01, &
                     0.100002328283153D+01, &
                     0.100002334902098D+01, &
                     0.100002354690067D+01, &
                     0.100002387440679D+01, &
                     0.100002432810739D+01, &
                     0.100002490321579D+01, &
                     0.100002559361283D+01, &
                     0.100002639188048D+01, &
                     0.100002728934954D+01, &
                     0.100002827616358D+01, &
                     0.100002934136059D+01, &
                     0.100003047297317D+01, &
                     0.100003165814696D+01, &
                     0.100003288327650D+01, &
                     0.100003413415694D+01, &
                     0.100003539614935D+01, &
                     0.100003665435715D+01, &
                     0.100003789381073D+01, &
                     0.100003909965731D+01, &
                     0.100004025735269D+01, &
                     0.100004135285193D+01, &
                     0.100004237279560D+01, &
                     0.100004330468876D+01, &
                     0.100004413706962D+01, &
                     0.100004485966530D+01, &
                     0.100004546353206D+01, &
                     0.100004594117799D+01, &
                     0.100004628666601D+01, &
                     0.100004649569578D+01, &
                     0.100004656566307D+01, &
                     0.100004656566307D+01, &
                     0.100004669781801D+01, &
                     0.100004709291716D+01, &
                     0.100004774686710D+01, &
                     0.100004865285865D+01, &
                     0.100004980139066D+01, &
                     0.100005118031024D+01, &
                     0.100005277487538D+01, &
                     0.100005456784524D+01, &
                     0.100005653960292D+01, &
                     0.100005866831405D+01, &
                     0.100006093012281D+01, &
                     0.100006329938527D+01, &
                     0.100006574893845D+01, &
                     0.100006825040188D+01, &
                     0.100007077450741D+01, &
                     0.100007329145218D+01, &
                     0.100007577126883D+01, &
                     0.100007818420686D+01, &
                     0.100008050111862D+01, &
                     0.100008269384348D+01, &
                     0.100008473558372D+01, &
                     0.100008660126597D+01, &
                     0.100008826788220D+01, &
                     0.100008971480478D+01, &
                     0.100009092407041D+01, &
                     0.100009188062857D+01, &
                     0.100009257255027D+01, &
                     0.100009299119400D+01, &
                     0.100009313132614D+01, &
                     0.100009313132614D+01, &
                     0.100009339515215D+01, &
                     0.100009418392394D+01, &
                     0.100009548952850D+01, &
                     0.100009729846576D+01, &
                     0.100009959188925D+01, &
                     0.100010234567788D+01, &
                     0.100010553055082D+01, &
                     0.100010911223719D+01, &
                     0.100011305171066D+01, &
                     0.100011730549625D+01, &
                     0.100012182605305D+01, &
                     0.100012656223333D+01, &
                     0.100013145981470D+01, &
                     0.100013646209936D+01, &
                     0.100014151057201D+01, &
                     0.100014654560588D+01, &
                     0.100015150720530D+01, &
                     0.100015633577211D+01, &
                     0.100016097288282D+01, &
                     0.100016536206328D+01, &
                     0.100016944954758D+01, &
                     0.100017318500856D+01, &
                     0.100017652224771D+01, &
                     0.100017941983296D+01, &
                     0.100018184167401D+01, &
                     0.100018375752579D+01, &
                     0.100018514341192D+01, &
                     0.100018598196126D+01, &
                     0.100018626265228D+01, &
                     0.100018626265228D+01, &
                     0.100018678927736D+01, &
                     0.100018836379350D+01, &
                     0.100019097013177D+01, &
                     0.100019458154364D+01, &
                     0.100019916066767D+01, &
                     0.100020465965428D+01, &
                     0.100021102037366D+01, &
                     0.100021817473187D+01, &
                     0.100022604511665D+01, &
                     0.100023454498889D+01, &
                     0.100024357962853D+01, &
                     0.100025304703635D+01, &
                     0.100026283898595D+01, &
                     0.100027284221425D+01, &
                     0.100028293973347D+01, &
                     0.100029301224391D+01, &
                     0.100030293962354D+01, &
                     0.100031260246900D+01, &
                     0.100032188366094D+01, &
                     0.100033066992675D+01, &
                     0.100033885337362D+01, &
                     0.100034633296573D+01, &
                     0.100035301592055D+01, &
                     0.100035881900080D+01, &
                     0.100036366968056D+01, &
                     0.100036750716623D+01, &
                     0.100037028325557D+01, &
                     0.100037196302080D+01, &
                     0.100037252530456D+01, &
                     0.100037252530456D+01, &
                     0.100037357644939D+01, &
                     0.100037671927613D+01, &
                     0.100038192197165D+01, &
                     0.100038913155859D+01, &
                     0.100039827399832D+01, &
                     0.100040925440104D+01, &
                     0.100042195739608D+01, &
                     0.100043624771554D+01, &
                     0.100045197103750D+01, &
                     0.100046895512333D+01, &
                     0.100048701126929D+01, &
                     0.100050593607681D+01, &
                     0.100052551353183D+01, &
                     0.100054551736988D+01, &
                     0.100056571369346D+01, &
                     0.100058586379968D+01, &
                     0.100060572716974D+01, &
                     0.100062506456813D+01, &
                     0.100064364119660D+01, &
                     0.100066122984736D+01, &
                     0.100067761400023D+01, &
                     0.100069259080983D+01, &
                     0.100070597393157D+01, &
                     0.100071759613817D+01, &
                     0.100072731168254D+01, &
                     0.100073499836755D+01, &
                     0.100074055928806D+01, &
                     0.100074392421649D+01, &
                     0.100074505060911D+01, &
                     0.100074505060911D+01, &
                     0.100074714886959D+01, &
                     0.100075342265153D+01, &
                     0.100076380896326D+01, &
                     0.100077820286428D+01, &
                     0.100079645760946D+01, &
                     0.100081838498991D+01, &
                     0.100084375598256D+01, &
                     0.100087230182097D+01, &
                     0.100090371558645D+01, &
                     0.100093765439429D+01, &
                     0.100097374221997D+01, &
                     0.100101157337804D+01, &
                     0.100105071663616D+01, &
                     0.100109071991891D+01, &
                     0.100113111553434D+01, &
                     0.100117142583757D+01, &
                     0.100121116923357D+01, &
                     0.100124986641183D+01, &
                     0.100128704670096D+01, &
                     0.100132225442918D+01, &
                     0.100135505517737D+01, &
                     0.100138504181458D+01, &
                     0.100141184021078D+01, &
                     0.100143511452836D+01, &
                     0.100145457200233D+01, &
                     0.100146996712809D+01, &
                     0.100148110518703D+01, &
                     0.100148784505081D+01, &
                     0.100149010121822D+01, &
                     0.100149010121822D+01, &
                     0.100149429113482D+01, &
                     0.100150681926102D+01, &
                     0.100152756071426D+01, &
                     0.100155630736539D+01, &
                     0.100159276801092D+01, &
                     0.100163656890664D+01, &
                     0.100168725489906D+01, &
                     0.100174429139323D+01, &
                     0.100180706736733D+01, &
                     0.100187489959484D+01, &
                     0.100194703817234D+01, &
                     0.100202267338407D+01, &
                     0.100210094387045D+01, &
                     0.100218094601043D+01, &
                     0.100226174438055D+01, &
                     0.100234238311630D+01, &
                     0.100242189797431D+01, &
                     0.100249932887576D+01, &
                     0.100257373270139D+01, &
                     0.100264419610451D+01, &
                     0.100270984811107D+01, &
                     0.100276987228216D+01, &
                     0.100282351822560D+01, &
                     0.100287011225745D+01, &
                     0.100290906703112D+01, &
                     0.100293988997180D+01, &
                     0.100296219037538D+01, &
                     0.100297568505410D+01, &
                     0.100298020243645D+01, &
                     0.100298020243645D+01, &
                     0.100298857590165D+01, &
                     0.100301361351724D+01, &
                     0.100305506687722D+01, &
                     0.100311252185226D+01, &
                     0.100318539874654D+01, &
                     0.100327295313264D+01, &
                     0.100337427786293D+01, &
                     0.100348830675945D+01, &
                     0.100361382042667D+01, &
                     0.100374945452612D+01, &
                     0.100389371072032D+01, &
                     0.100404497035260D+01, &
                     0.100420151079455D+01, &
                     0.100436152427346D+01, &
                     0.100452313889348D+01, &
                     0.100468444148787D+01, &
                     0.100484350188495D+01, &
                     0.100499839813490D+01, &
                     0.100514724222625D+01, &
                     0.100528820581640D+01, &
                     0.100541954550811D+01, &
                     0.100553962722075D+01, &
                     0.100564694923018D+01, &
                     0.100574016348216D+01, &
                     0.100581809482091D+01, &
                     0.100587975781456D+01, &
                     0.100592437090370D+01, &
                     0.100595136764524D+01, &
                     0.100596040487289D+01, &
                     0.100596040487289D+01, &
                     0.100597716627364D+01, &
                     0.100602728468665D+01, &
                     0.100611026256724D+01, &
                     0.100622527029612D+01, &
                     0.100637114631855D+01, &
                     0.100654639863266D+01, &
                     0.100674920867118D+01, &
                     0.100697743862518D+01, &
                     0.100722864313084D+01, &
                     0.100750008601400D+01, &
                     0.100778876250505D+01, &
                     0.100809142703738D+01, &
                     0.100840462645811D+01, &
                     0.100872473823122D+01, &
                     0.100904801301117D+01, &
                     0.100937062081313D+01, &
                     0.100968869990192D+01, &
                     0.100999840746024D+01, &
                     0.101029597107207D+01, &
                     0.101057774006124D+01, &
                     0.101084023575438D+01, &
                     0.101108019978336D+01, &
                     0.101129463960389D+01, &
                     0.101148087047733D+01, &
                     0.101163655324184D+01, &
                     0.101175972728230D+01, &
                     0.101184883819617D+01, &
                     0.101190275974165D+01, &
                     0.101192080974579D+01, &
                     0.101192080974579D+01, &
                     0.101195446057276D+01, &
                     0.101205507665813D+01, &
                     0.101222164853247D+01, &
                     0.101245249324978D+01, &
                     0.101274525512469D+01, &
                     0.101309690947116D+01, &
                     0.101350377150737D+01, &
                     0.101396151257463D+01, &
                     0.101446518551773D+01, &
                     0.101500926056794D+01, &
                     0.101558767245067D+01, &
                     0.101619387879617D+01, &
                     0.101682092933216D+01, &
                     0.101746154482776D+01, &
                     0.101810820436123D+01, &
                     0.101875323920339D+01, &
                     0.101938893143841D+01, &
                     0.102000761536699D+01, &
                     0.102060177973906D+01, &
                     0.102116416892463D+01, &
                     0.102168788123711D+01, &
                     0.102216646276009D+01, &
                     0.102259399518391D+01, &
                     0.102296517632439D+01, &
                     0.102327539216643D+01, &
                     0.102352077944395D+01, &
                     0.102369827793388D+01, &
                     0.102380567180175D+01, &
                     0.102384161949158D+01, &
                     0.102384161949158D+01, &
                     0.102390951714174D+01, &
                     0.102411251211690D+01, &
                     0.102444851070429D+01, &
                     0.102491402502249D+01, &
                     0.102550417855974D+01, &
                     0.102621271950402D+01, &
                     0.102703204623041D+01, &
                     0.102795324915137D+01, &
                     0.102896617236760D+01, &
                     0.103005949736791D+01, &
                     0.103122084963508D+01, &
                     0.103243692762291D+01, &
                     0.103369365233101D+01, &
                     0.103497633471390D+01, &
                     0.103626985745727D+01, &
                     0.103755886723504D+01, &
                     0.103882797339283D+01, &
                     0.104006194904210D+01, &
                     0.104124593074224D+01, &
                     0.104236561324820D+01, &
                     0.104340743616575D+01, &
                     0.104435875975053D+01, &
                     0.104520802748490D+01, &
                     0.104594491344819D+01, &
                     0.104656045284981D+01, &
                     0.104704715441344D+01, &
                     0.104739909358129D+01, &
                     0.104761198575277D+01, &
                     0.104768323898315D+01, &
                     0.104768323898315D+01, &
                     0.104782133463030D+01, &
                     0.104823412126563D+01, &
                     0.104891710251107D+01, &
                     0.104986280140041D+01, &
                     0.105106079199071D+01, &
                     0.105249775450095D+01, &
                     0.105415756207908D+01, &
                     0.105602140644556D+01, &
                     0.105806796753855D+01, &
                     0.106027362938381D+01, &
                     0.106261274127182D+01, &
                     0.106505792041863D+01, &
                     0.106758038994200D+01, &
                     0.107015034436819D+01, &
                     0.107273733402188D+01, &
                     0.107531065946586D+01, &
                     0.107783976751535D+01, &
                     0.108029464109827D+01, &
                     0.108264617621765D+01, &
                     0.108486654036642D+01, &
                     0.108692950784536D+01, &
                     0.108881076846885D+01, &
                     0.109048820706262D+01, &
                     0.109194215193640D+01, &
                     0.109315559114380D+01, &
                     0.109411435582429D+01, &
                     0.109480727026987D+01, &
                     0.109522626858838D+01, &
                     0.109536647796631D+01, &
                     0.109536647796631D+01, &
                     0.109565068055070D+01, &
                     0.109649989493231D+01, &
                     0.109790396904982D+01, &
                     0.109984607945475D+01, &
                     0.110230288633927D+01, &
                     0.110524476755142D+01, &
                     0.110863614229226D+01, &
                     0.111243589142621D+01, &
                     0.111659787533130D+01, &
                     0.112107154332923D+01, &
                     0.112580262230830D+01, &
                     0.113073386718843D+01, &
                     0.113580585290523D+01, &
                     0.114095778668146D+01, &
                     0.114612832023820D+01, &
                     0.115125634381824D+01, &
                     0.115628174693491D+01, &
                     0.116114613414939D+01, &
                     0.116579348753928D+01, &
                     0.117017077058586D+01, &
                     0.117422847081695D+01, &
                     0.117792108062658D+01, &
                     0.118120751724512D+01, &
                     0.118405148389083D+01, &
                     0.118642177475713D+01, &
                     0.118829252674974D+01, &
                     0.118964342085509D+01, &
                     0.119045983576227D+01, &
                     0.119073295593262D+01, &
                     0.119073295593262D+01, &
                     0.119132698543447D+01, &
                     0.119310089249382D+01, &
                     0.119603027382615D+01, &
                     0.120007492816088D+01, &
                     0.120517956316194D+01, &
                     0.121127478070291D+01, &
                     0.121827832553995D+01, &
                     0.122609656852371D+01, &
                     0.123462618103399D+01, &
                     0.124375594563308D+01, &
                     0.125336864143105D+01, &
                     0.126334294227058D+01, &
                     0.127355527107242D+01, &
                     0.128388156305859D+01, &
                     0.129419890222035D+01, &
                     0.130438700755100D+01, &
                     0.131432955684288D+01, &
                     0.132391534538505D+01, &
                     0.133303928429119D+01, &
                     0.134160324840566D+01, &
                     0.134951678698428D+01, &
                     0.135669771195690D+01, &
                     0.136307257892185D+01, &
                     0.136857707544672D+01, &
                     0.137315633005887D+01, &
                     0.137676515373815D+01, &
                     0.137936822395133D+01, &
                     0.138094021940997D+01, &
                     0.138146591186523D+01, &
                     0.138146591186523D+01, &
                     0.138272810558301D+01, &
                     0.138649367568476D+01, &
                     0.139270024291729D+01, &
                     0.140124598169592D+01, &
                     0.141199269281076D+01, &
                     0.142476984050373D+01, &
                     0.143937931323955D+01, &
                     0.145560063568618D+01, &
                     0.147319635466327D+01, &
                     0.149191734389366D+01, &
                     0.151150781595592D+01, &
                     0.153170988612211D+01, &
                     0.155226759223763D+01, &
                     0.157293032946533D+01, &
                     0.159345570340066D+01, &
                     0.161361183757945D+01, &
                     0.163317919195240D+01, &
                     0.165195195921019D+01, &
                     0.166973910830179D+01, &
                     0.168636514153925D+01, &
                     0.170167062545724D+01, &
                     0.171551254775478D+01, &
                     0.172776454437017D+01, &
                     0.173831703278350D+01, &
                     0.174707728043116D+01, &
                     0.175396943084028D+01, &
                     0.175893450477628D+01, &
                     0.176193038927058D+01, &
                     0.176293182373047D+01, &
                     0.176293182373047D+01, &
                     0.176564286694454D+01, &
                     0.177372005645614D+01, &
                     0.178699837216681D+01, &
                     0.180521178444865D+01, &
                     0.182800567156310D+01, &
                     0.185495225211076D+01, &
                     0.188556753788460D+01, &
                     0.191932839187586D+01, &
                     0.195568853835757D+01, &
                     0.199409272966140D+01, &
                     0.203398864215525D+01, &
                     0.207483638976021D+01, &
                     0.211611577541990D+01, &
                     0.215733154432817D+01, &
                     0.219801696989445D+01, &
                     0.223773611446484D+01, &
                     0.227608508224824D+01, &
                     0.231269253875014D+01, &
                     0.234721972151142D+01, &
                     0.237936011876969D+01, &
                     0.240883894992880D+01, &
                     0.243541254610567D+01, &
                     0.245886770067967D+01, &
                     0.247902103804551D+01, &
                     0.249571843266572D+01, &
                     0.250883449896018D+01, &
                     0.251827216455803D+01, &
                     0.252396233410609D+01, &
                     0.252586364746094D+01, &
                     0.252586364746094D+01, &
                     0.253167688860071D+01, &
                     0.254896975540049D+01, &
                     0.257731221456166D+01, &
                     0.261602085129962D+01, &
                     0.266420303074566D+01, &
                     0.272080885678178D+01, &
                     0.278468433772325D+01, &
                     0.285462045565742D+01, &
                     0.292939476825765D+01, &
                     0.300780410266841D+01, &
                     0.308868841556147D+01, &
                     0.317094684638664D+01, &
                     0.325354743206306D+01, &
                     0.333553202180912D+01, &
                     0.341601778229227D+01, &
                     0.349419643693284D+01, &
                     0.356933211849552D+01, &
                     0.364075847429041D+01, &
                     0.370787546624720D+01, &
                     0.377014615659739D+01, &
                     0.382709365915031D+01, &
                     0.387829835864702D+01, &
                     0.392339544868338D+01, &
                     0.396207280543774D+01, &
                     0.399406919449255D+01, &
                     0.401917279728607D+01, &
                     0.403722003917193D+01, &
                     0.404809470058214D+01, &
                     0.405172729492188D+01, &
                     0.405172729492188D+01, &
                     0.406402082977662D+01, &
                     0.410053778948581D+01, &
                     0.416022062358699D+01, &
                     0.424141073028447D+01, &
                     0.434198166664407D+01, &
                     0.445948854744914D+01, &
                     0.459131164063186D+01, &
                     0.473477902630729D+01, &
                     0.488726123050600D+01, &
                     0.504623735793827D+01, &
                     0.520933635967275D+01, &
                     0.537435883071529D+01, &
                     0.553928485229371D+01, &
                     0.570227263615322D+01, &
                     0.586165166093097D+01, &
                     0.601591294620142D+01, &
                     0.616369823788700D+01, &
                     0.630378921725333D+01, &
                     0.643509737681630D+01, &
                     0.655665489197294D+01, &
                     0.666760661690361D+01, &
                     0.676720321302336D+01, &
                     0.685479535163237D+01, &
                     0.692982890044344D+01, &
                     0.699184099328238D+01, &
                     0.704045688480234D+01, &
                     0.707538750197330D+01, &
                     0.709642761796854D+01, &
                     0.710345458984375D+01, &
                     0.710345458984375D+01, &
                     0.712898332413999D+01, &
                     0.720473123783497D+01, &
                     0.732827303408514D+01, &
                     0.749584253080731D+01, &
                     0.770267739888486D+01, &
                     0.794339241549997D+01, &
                     0.821232199207819D+01, &
                     0.850379621062915D+01, &
                     0.881233865152432D+01, &
                     0.913279153349242D+01, &
                     0.946038243491060D+01, &
                     0.979074896402138D+01, &
                     0.101199360132357D+02, &
                     0.104443769908930D+02, &
                     0.107608670822496D+02, &
                     0.110665337928737D+02, &
                     0.113588079404209D+02, &
                     0.116353968231664D+02, &
                     0.118942603647849D+02, &
                     0.121335904701419D+02, &
                     0.123517935070912D+02, &
                     0.125474756664057D+02, &
                     0.127194308851829D+02, &
                     0.128666310083941D+02, &
                     0.129882178834287D+02, &
                     0.130834971181631D+02, &
                     0.131519332753360D+02, &
                     0.131931463199519D+02, &
                     0.132069091796875D+02, &
                     0.132069091796875D+02, &
                     0.132591194635718D+02, &
                     0.134139272905853D+02, &
                     0.136660758305702D+02, &
                     0.140074531941431D+02, &
                     0.144278909137090D+02, &
                     0.149160094499923D+02, &
                     0.154599716304101D+02, &
                     0.160480674716776D+02, &
                     0.166691131732133D+02, &
                     0.173126859831722D+02, &
                     0.179692334974645D+02, &
                     0.186300972562267D+02, &
                     0.192874839533251D+02, &
                     0.199344087500137D+02, &
                     0.205646270599299D+02, &
                     0.211725648484586D+02, &
                     0.217532530464901D+02, &
                     0.223022687845992D+02, &
                     0.228156843898291D+02, &
                     0.232900240866142D+02, &
                     0.237222278292245D+02, &
                     0.241096214738796D+02, &
                     0.244498924491712D+02, &
                     0.247410701254447D+02, &
                     0.249815101704095D+02, &
                     0.251698822826361D+02, &
                     0.253051608022386D+02, &
                     0.253866178019575D+02, &
                     0.254138183593750D+02, &
                     0.254138183593750D+02, &
                     0.255195267992590D+02, &
                     0.258328352681913D+02, &
                     0.263427584579998D+02, &
                     0.270324057024198D+02, &
                     0.278807087211529D+02, &
                     0.288642280840487D+02, &
                     0.299587361586212D+02, &
                     0.311404194628242D+02, &
                     0.323866753768065D+02, &
                     0.336765603552101D+02, &
                     0.349909787594548D+02, &
                     0.363127000154435D+02, &
                     0.376262748757950D+02, &
                     0.389179012258676D+02, &
                     0.401752720930330D+02, &
                     0.413874251780706D+02, &
                     0.425446041482958D+02, &
                     0.436381361995771D+02, &
                     0.446603270190587D+02, &
                     0.456043724487356D+02, &
                     0.464642852734054D+02, &
                     0.472348352404610D+02, &
                     0.479115004130314D+02, &
                     0.484904281100403D+02, &
                     0.489684039083995D+02, &
                     0.493428274250085D+02, &
                     0.496116938345814D+02, &
                     0.497735803026837D+02, &
                     0.498276367187500D+02, &
                     0.498276367187500D+02, &
                     0.500404185653889D+02, &
                     0.506709451241949D+02, &
                     0.516967354516298D+02, &
                     0.530832902126245D+02, &
                     0.547876911080444D+02, &
                     0.567623388473297D+02, &
                     0.589581992199333D+02, &
                     0.613272400021151D+02, &
                     0.638240200551633D+02, &
                     0.664065605099864D+02, &
                     0.690366893958513D+02, &
                     0.716800435219156D+02, &
                     0.743058733129875D+02, &
                     0.768867527509485D+02, &
                     0.793982594510897D+02, &
                     0.818186625637876D+02, &
                     0.841286378910208D+02, &
                     0.863110182434260D+02, &
                     0.883505804990993D+02, &
                     0.902338673542942D+02, &
                     0.919490401686329D+02, &
                     0.934857588058918D+02, &
                     0.948350844569658D+02, &
                     0.959894018059147D+02, &
                     0.969423573925416D+02, &
                     0.976888115436852D+02, &
                     0.982248017449617D+02, &
                     0.985475157865710D+02, &
                     0.986552734375000D+02, &
                     0.986552734375000D+02, &
                     0.990822433881177D+02, &
                     0.100347322073693D+03, &
                     0.102405016171537D+03, &
                     0.105185581277241D+03, &
                     0.108602371958246D+03, &
                     0.112559447931337D+03, &
                     0.116958148487609D+03, &
                     0.121701998090751D+03, &
                     0.126699878495252D+03, &
                     0.131867743820740D+03, &
                     0.137129275010602D+03, &
                     0.142415849895469D+03, &
                     0.147666124386204D+03, &
                     0.152825430261725D+03, &
                     0.157845119200052D+03, &
                     0.162681927388201D+03, &
                     0.167297398331762D+03, &
                     0.171657378883509D+03, &
                     0.175731590585769D+03, &
                     0.179493271683294D+03, &
                     0.182918882161275D+03, &
                     0.185987863294641D+03, &
                     0.188682443465482D+03, &
                     0.190987482825228D+03, &
                     0.192890350412873D+03, &
                     0.194380828411863D+03, &
                     0.195451039250058D+03, &
                     0.196095392185899D+03, &
                     0.196310546875000D+03, &
                     0.196310546875000D+03, &
                     0.197165914414290D+03, &
                     0.199700157345842D+03, &
                     0.203821746551970D+03, &
                     0.209390433034766D+03, &
                     0.216232103038247D+03, &
                     0.224154123334171D+03, &
                     0.232958573417629D+03, &
                     0.242452088232637D+03, &
                     0.252452195390785D+03, &
                     0.262790716945047D+03, &
                     0.273315042593022D+03, &
                     0.283888035454548D+03, &
                     0.294387165578700D+03, &
                     0.304703283212886D+03, &
                     0.314739290656954D+03, &
                     0.324408860277427D+03, &
                     0.333635272730160D+03, &
                     0.342350404375929D+03, &
                     0.350493867238099D+03, &
                     0.358012291585652D+03, &
                     0.364858735397280D+03, &
                     0.370992203361916D+03, &
                     0.376377258713162D+03, &
                     0.380983712906265D+03, &
                     0.384786380265389D+03, &
                     0.387764886905412D+03, &
                     0.389903525298682D+03, &
                     0.391191147749383D+03, &
                     0.391621093750000D+03, &
                     0.391621093750000D+03, &
                     0.393333267347447D+03, &
                     0.398405869289398D+03, &
                     0.406655293223968D+03, &
                     0.417800273585348D+03, &
                     0.431491752810837D+03, &
                     0.447343706219330D+03, &
                     0.464959690282331D+03, &
                     0.483952559463633D+03, &
                     0.503957133152309D+03, &
                     0.524636970282076D+03, &
                     0.545686879557213D+03, &
                     0.566832696329901D+03, &
                     0.587829520532816D+03, &
                     0.608459240805830D+03, &
                     0.628527861954215D+03, &
                     0.647862929765958D+03, &
                     0.666311200075487D+03, &
                     0.683736608975723D+03, &
                     0.700018550031233D+03, &
                     0.715050438024817D+03, &
                     0.728738527294404D+03, &
                     0.741000949653365D+03, &
                     0.751766938273165D+03, &
                     0.760976207401328D+03, &
                     0.768578462075733D+03, &
                     0.774533016383188D+03, &
                     0.778808502964932D+03, &
                     0.781382660271209D+03, &
                     0.782242187500000D+03, &
                     0.782242187500000D+03, &
                     0.785667978702648D+03, &
                     0.795817314057572D+03, &
                     0.812322429889949D+03, &
                     0.834620023768783D+03, &
                     0.862011146904421D+03, &
                     0.893722988907681D+03, &
                     0.928962058479187D+03, &
                     0.966953648403494D+03, &
                     0.100696716166369D+04, &
                     0.104832963147092D+04, &
                     0.109043070530138D+04, &
                     0.113272216380572D+04, &
                     0.117471436749431D+04, &
                     0.121597128252476D+04, &
                     0.125610511934670D+04, &
                     0.129477117112498D+04, &
                     0.133166314449157D+04, &
                     0.136650909536301D+04, &
                     0.139906798067640D+04, &
                     0.142912678447536D+04, &
                     0.145649815400263D+04, &
                     0.148101847546883D+04, &
                     0.150254632203864D+04, &
                     0.152096121363644D+04, &
                     0.153616263679928D+04, &
                     0.154806928161230D+04, &
                     0.155661846109447D+04, &
                     0.156176568601542D+04, &
                     0.156348437500000D+04, &
                     0.156348437500000D+04, &
                     0.157033740416962D+04, &
                     0.159064021408069D+04, &
                     0.162365672497687D+04, &
                     0.166825955882024D+04, &
                     0.172304998255386D+04, &
                     0.178648161296558D+04, &
                     0.185696686235031D+04, &
                     0.193295589977556D+04, &
                     0.201298729543331D+04, &
                     0.209571503135048D+04, &
                     0.217991843292792D+04, &
                     0.226450117183321D+04, &
                     0.234848413013796D+04, &
                     0.243099542940256D+04, &
                     0.251125969168354D+04, &
                     0.258858770516681D+04, &
                     0.266236707830019D+04, &
                     0.273205410682725D+04, &
                     0.279716687457548D+04, &
                     0.285727950422681D+04, &
                     0.291201742892684D+04, &
                     0.296105354375488D+04, &
                     0.300410510192081D+04, &
                     0.304093123474894D+04, &
                     0.307133099181045D+04, &
                     0.309514181521440D+04, &
                     0.311223837875522D+04, &
                     0.312253173785492D+04, &
                     0.312596875000000D+04, &
                     0.312596875000000D+04, &
                     0.313967625648517D+04, &
                     0.318028601937929D+04, &
                     0.324632532604472D+04, &
                     0.333553864629481D+04, &
                     0.344512767762386D+04, &
                     0.357199889047001D+04, &
                     0.371297650388395D+04, &
                     0.386496043932085D+04, &
                     0.402502760139934D+04, &
                     0.419048586991262D+04, &
                     0.435889392629539D+04, &
                     0.452805922447528D+04, &
                     0.469602368983688D+04, &
                     0.486104375492158D+04, &
                     0.502156886517289D+04, &
                     0.517622079894589D+04, &
                     0.532377496843506D+04, &
                     0.546314414912483D+04, &
                     0.559336467869843D+04, &
                     0.571358495717123D+04, &
                     0.582305598954208D+04, &
                     0.592112368866449D+04, &
                     0.600722266786800D+04, &
                     0.608087128130021D+04, &
                     0.614166770461798D+04, &
                     0.618928688399240D+04, &
                     0.622347821477834D+04, &
                     0.624406384170965D+04, &
                     0.625093750000000D+04, &
                     0.625093750000000D+04, &
                     0.627835396180777D+04, &
                     0.635957763260931D+04, &
                     0.649166253364183D+04, &
                     0.667009682994918D+04, &
                     0.688928307967390D+04, &
                     0.714303346020215D+04, &
                     0.742499580388151D+04, &
                     0.772896953684780D+04, &
                     0.804910823258752D+04, &
                     0.838002756647795D+04, &
                     0.871684493213928D+04, &
                     0.905517534807460D+04, &
                     0.939110282643707D+04, &
                     0.972114042184661D+04, &
                     0.100421872265593D+05, &
                     0.103514869993551D+05, &
                     0.106465907599641D+05, &
                     0.109253242434066D+05, &
                     0.111857602951074D+05, &
                     0.114261958697820D+05, &
                     0.116451331161567D+05, &
                     0.118412639826529D+05, &
                     0.120134578028544D+05, &
                     0.121607513765659D+05, &
                     0.122823411316261D+05, &
                     0.123775770223353D+05, &
                     0.124459578871755D+05, &
                     0.124871280495070D+05, &
                     0.125008750000000D+05, &
                     0.125008750000000D+05, &
                     0.125557093727928D+05, &
                     0.127181608604211D+05, &
                     0.129823369516715D+05, &
                     0.133392132016916D+05, &
                     0.137775938898073D+05, &
                     0.142851026071009D+05, &
                     0.148490344123996D+05, &
                     0.154569877411764D+05, &
                     0.160972695046301D+05, &
                     0.167591109693681D+05, &
                     0.174327469533976D+05, &
                     0.181094076044732D+05, &
                     0.187812611082965D+05, &
                     0.194413337636850D+05, &
                     0.200834239565779D+05, &
                     0.207020194066333D+05, &
                     0.212922223486821D+05, &
                     0.218496844368383D+05, &
                     0.223705515320276D+05, &
                     0.228514176983811D+05, &
                     0.232892873720912D+05, &
                     0.236815445727244D+05, &
                     0.240259280743806D+05, &
                     0.243205115681841D+05, &
                     0.245636879863419D+05, &
                     0.247541572994163D+05, &
                     0.248909172321462D+05, &
                     0.249732564651458D+05, &
                     0.250007500000000D+05, &
                     0.250007500000000D+05, &
                     0.251104201949490D+05, &
                     0.254353273165597D+05, &
                     0.259636857887237D+05, &
                     0.266774459469898D+05, &
                     0.275542155126487D+05, &
                     0.285692409041413D+05, &
                     0.296971116331728D+05, &
                     0.309130241539382D+05, &
                     0.321935920529177D+05, &
                     0.335172777794490D+05, &
                     0.348645509999082D+05, &
                     0.362178721216079D+05, &
                     0.375615776765313D+05, &
                     0.388817204513680D+05, &
                     0.401658974203403D+05, &
                     0.414030842244438D+05, &
                     0.425834855290145D+05, &
                     0.436984048261591D+05, &
                     0.447401340079593D+05, &
                     0.457018613572892D+05, &
                     0.465775958853328D+05, &
                     0.473621057539328D+05, &
                     0.480508686182164D+05, &
                     0.486400319519764D+05, &
                     0.491263816961243D+05, &
                     0.495073178537817D+05, &
                     0.497808359221753D+05, &
                     0.499455132964462D+05, &
                     0.500005000000000D+05, &
                     0.500005000000000D+05, &
                     0.502198418393502D+05, &
                     0.508696602293366D+05, &
                     0.519263834638734D+05, &
                     0.533539114391611D+05, &
                     0.551074587604351D+05, &
                     0.571375175008067D+05, &
                     0.593932660778249D+05, &
                     0.618250969827884D+05, &
                     0.643862371533144D+05, &
                     0.670336114033129D+05, &
                     0.697281590973150D+05, &
                     0.724348011585752D+05, &
                     0.751222108142610D+05, &
                     0.777624938283821D+05, &
                     0.803308443491094D+05, &
                     0.828052138613731D+05, &
                     0.851660118907048D+05, &
                     0.873958456057811D+05, &
                     0.894792989605956D+05, &
                     0.914027486757764D+05, &
                     0.931542129123466D+05, &
                     0.947232281167544D+05, &
                     0.961007497062068D+05, &
                     0.972790727197631D+05, &
                     0.982517691158391D+05, &
                     0.990136389625832D+05, &
                     0.995606733022726D+05, &
                     0.998900269590546D+05, &
                     0.100000000000000D+06  /

do int=1,nints-1
if(x .lt. ab(2,int)) exit
end do
a = ab(1,int)
b = ab(2,int)
xx   = (2*x - (b+a) ) /(b-a)
sum1=0
sum2=0
dd1 = 1.0d0
do i=1,k
 dd=1.0d0
 if (i .eq. 1 .OR. i .eq. k) dd = 0.5d0
diff = xx-xscheb(i)
if(abs(diff) .le. eps0 ) then
val = vals(i,int)
return
 endif
dd   = (dd1*dd)/diff
dd1  = - dd1
sum1 = sum1+dd*vals(i,int)
sum2 = sum2+dd
dd   = - dd
end do
val = sum1/sum2
end subroutine



end module prolates
