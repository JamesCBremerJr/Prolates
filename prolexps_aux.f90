subroutine chebpw_eval2(nints,ab,k,coefs1,coefs2,x,val1,val2)
implicit double precision (a-h,o-z)

integer, intent(in)           :: nints,k
double precision, intent(in)  :: ab(2,nints),coefs1(k,nints),coefs2(k,nints)
double precision, intent(out) :: val1,val2

!
!  Evaluate two functions represented via piecewise chebyshev expansions on
!  the same collection of subintervals.
!
!  Input parameters:
!    (nints,ab) - the 
!    k - an integer specifying the order of the Chebyshev expansions; on
!    
!    coefs1 - a (k,nints) array whose jth column specified the coefficients
!     of the first function's Chebyshev expansion on the jth subinterval
!    coefs2 - a (k,nints) array whose jth column specified the coefficients
!     of the second function's Chebyshev expansion on the jth subinterval

!    x - the point at which to evaluate the functions
!
!
!  Output parameters:
!    val1 - the value of the first function at the point x
!    val2 - the value of the second function at the point x

! 


double precision :: pols(k)

intl   = 1
intr   = nints

do int = intl,intr-1
b = ab(2,int)
if (x .le. b) exit
end do

a = ab(1,int)
b = ab(2,int)



!
!  Evaluate the Chebyshev expansions
!


xx = (x - (b+a)/2.0d0) * 2.0d0/(b-a)
call chebs(xx,k-1,pols)

val1 = 0
val2 = 0
do i=1,k
val1 = val1 + coefs1(i,int)*pols(i)
val2 = val2 + coefs2(i,int)*pols(i)
end do


end subroutine



subroutine chebs(x,n,pols)
implicit double precision (a-h,o-z)

integer          :: n 
double precision :: pols(0:n),x

!
!  Evaluate the Chebyshev polynomials of degree 0 through n at a specified point
!  using the standard 3-term recurrence relation.
!
!  Input parameters:
!
!    x - point at which to evaluate the polynomials
!    n - the order of the polynomials to evaluate
!
!  Output parameters:
!
!    pols - this user-supplied and allocated array of length n+1 will
!      contain the values of the polynomials of order 0 through n upon
!      return
!

if (x .eq. 1.0d0) then
do i=0,n
pols(i) = 1.0d0
end do
return
endif

if (x .eq. -1.0d0) then
pols(0) = 1.0d0
do i=1,n
pols(i) = -pols(i-1)
end do
return
endif

pols(0) = 1.0d0
if (n .eq. 0) return

pols(1) = x
if (n .eq. 1) return

xx1 = 1.0d0
xx2 = x

do i=1,n-1
xx        = 2*x*xx2 - xx1
pols(i+1) = xx
xx1       = xx2
xx2       = xx
end do

end subroutine


subroutine tensorpw_eval4(n,nintsab,ab,nintscd,cd,coefs1,coefs2,coefs3,coefs4,x,y, &
  val1,val2,val3,val4)
implicit double precision (a-h,o-z)

integer, intent(in)           :: n, nintsab, nintscd
double precision, intent(in)  :: ab(2,nintsab), cd(2,nintscd), x, y
double precision, intent(in)  :: coefs1(n*(n+1)/2,nintsab,nintscd)
double precision, intent(in)  :: coefs2(n*(n+1)/2,nintsab,nintscd)
double precision, intent(in)  :: coefs3(n*(n+1)/2,nintsab,nintscd)
double precision, intent(in)  :: coefs4(n*(n+1)/2,nintsab,nintscd)
double precision, intent(out) :: val1, val2, val3, val4
!
!  Evaluate 4 piecewise bivariate Chebyshev expansions on a collection of rectangles 
!  of the form (2) at a specified point.
!
!  Input parameters:
!
!    n - the number of points in the Chebyshev grid, so that the order of the
!       expansion is n-1
!    ab - a (2,nintsab) array specifying the intervals (a_i,b_i) in (2)
!    cd - a (2,nintscd) array specifying the intervals (c_j,d_j) in (2)
!    coefs? - the (ncoefs,nintsab,nintscd) arrays which specify the coefficient
!      coefficients in each of the four expansions (See tensorpw_eval for details)
!    (x,y) - the point at which to evaluate the expansions
!
!  Output parameters:
!   val? - the value of the expansions at the specified point
!

double precision :: polsx(n),polsy(n)


do intab=1,nintsab-1
if (ab(1,intab+1) .ge. x) exit
end do

do intcd=1,nintscd-1
if (cd(1,intcd+1) .ge. y) exit
end do


a = ab(1,intab)
b = ab(2,intab)

c = cd(1,intcd)
d = cd(2,intcd)


xx = (x-(b+a)/2)*2/(b-a)
yy = (y-(d+c)/2)*2/(d-c)


call chebs(xx,n-1,polsx)
call chebs(yy,n-1,polsy)


val1  = 0
val2  = 0
val3  = 0
val4  = 0

idx  = 1
do j=1,n
do i=1,n-j+1

val1  = val1  + coefs1(idx,intab,intcd) * polsx(i)*polsy(j)
val2  = val2  + coefs2(idx,intab,intcd) * polsx(i)*polsy(j)
val3  = val3  + coefs3(idx,intab,intcd) * polsx(i)*polsy(j)
val4  = val4  + coefs4(idx,intab,intcd) * polsx(i)*polsy(j)

idx   = idx+1
end do
end do

end subroutine



