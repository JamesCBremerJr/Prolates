!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  
!  This module contains code for representing functions of two variables given on
!  a rectangle [a,b] x [c,d] as a bivariate Chebyshev expansions -- that is,
!  and an expansion of the form
!                    
!                n-1  n-1-j                  
!     f(x,y) =   sum   sum     a_ij T_{1,i}(x) T_{2,j}(y)                                     (1)
!                j=0   i=0
!
!
!  where T_{1,i} denotes the Chebyshev polynomial of degree i on the interval [a,b]
!  and T_{2,i} denote the Chebyshev polynomial of degree i on the interval [c,d].
!
!  It also contains code representing functions f(x,y) via piecwise bivariate
!  Chebyshev expansions given over collections of rectangles of the form
! 
!    [a_i, b_i] x [c_j, d_j]   i=1,...,m1  j=1,..,m2.                                        (2)
!
!  The following subroutines should be regarded as public:
!
!    tensor_umatrix - return the matrix which takes the values of a funtion f at
!      nodes of a tensor product of Clenshaw-Curtis quadrature formulas to the 
!      coefficients in the expansion (1)
!
!    tensor_coefs - a utility routine for applying the matrix returned by
!      tensor_umatrix to the vector of values; its raison d'etre is that
!      it provides a mechanism for reshaping the array of values implicitly
!
!    tensor_eval - calculate the value of an expansion of the form (1)  at 
!      a specified point
!
!    tensor_evalder - calculate the value of an expansion of the form (1) and
!      its derivative at a specified point
!
!    tensorpw_eval - evaluate a function represented as a piecewise bivariate
!      Chebyshev expansion over a collection of rectanges of the form (2) at
!      a specified point
!
!    tensorpw_evalder - evaluate a function represented as a piecewise bivariate
!      Chebyshev expansion and over a collection of rectanges of the form (2) and
!      its first derivatives at a specified point
!
!    tensorpw_eval4 - this is a rather specialized routine which evaluates four
!      different functions represented via piecewise bivariate Chebyshev expansions
!      given on the same collection of rectangles (2) at a specified point;
!      its raise d'etre is that it is faster than calling tensorpw_eval four times
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module tensor

use chebyshev

contains

subroutine tensor_umatrix(n,ncoefs,u)
implicit double precision (a-h,o-z)

integer, intent(in)                         :: n
integer, intent(out)                        :: ncoefs
double precision, allocatable, intent(out)  :: u(:,:)

!
!  Return the matrix which takes the vector of values
!
!      f( xs(1)    , xs(1) )
!      f( xs(2)    , xs(1) )
!               .
!               .
!      f( xs(n)    , xs(1) )
!      f( xs(1)    , xs(2) )
!      f( xs(2)    , xs(2) )                                                      (2)
!               .
!               .
!
!      f( xs(1)    , xs(n) )
!      f( xs(2)    , xs(n) )
!               .
!               .
!      f( xs(n )   , xs(n) ),
!
!  where xs(1), ..., xs(n) are the nodes of the n-point Clenshaw-Curtis quadrature
!  formula to the vector of coefficients
!
!      a_00
!      a_01
!      a_10
!      a_02                                                                       (3)

!      a_11
!      a_20
!       .
!       .
!       .
!    
!
!  in the expansion (1).  Note that there are (n^2 + n)/2 such coefficients, so that
!  the matrix u is dimensioned ( (n**2 + n)/2 , n*n ).  Also note tht the coefficients
!  in (3) are sorted by the degree of the corresponding term in the expansion (1).
!
!  Input parameters:
!    n - the number of points in the Chebyshev grid on [-1,1] 
!    
!  Output parameters:
!    u - the ((n+1)*n/2, n*n) matrix which takes the vector of values (2) 
!     to the vector of coefficients (3)
!
!

double precision, allocatable :: polsx(:),polsy(:),xs(:)

ncoefs = n*(n+1)/2
allocate(u(ncoefs,n*n))


!
!  Fetch the Chebyshev nodes and allocate memory for the values of the
!  Chebyshev polynomials
!
allocate(polsx(n+1), polsy(n+1),xs(n))
call chebpts(n,xs)


!
!  Use the identity
!
!
!              n  ''                        { \delta_ij * (n)   if i = 0 or i = n
!             \sum   T_i(x_k) T_j(x_k)  =   {
!              k=0                          { \delta_ij * (n/2) if 0 < i < n
!
!  to form the matrix u0 which takes the values of an expansion of the form
!
!               n-1  n-1-i
!      f(x,y) = sum  sum   a_{ij} T_i(x) T_j(x) 
!               i=0  j=0
!
!  to the coefficients in (1).
!


do i1=1,n
x    = xs(i1)
call chebs(x,n,polsx)
if (i1 .eq. 1 .OR. i1 .eq. n) polsx = polsx/2.0d0

do i2=1,n
y    = xs(i2)
call chebs(y,n,polsy)
if (i2 .eq. 1 .OR. i2 .eq. n) polsy = polsy/2.0d0

idx = 1

do j2=1,n
do j1=1,n-j2+1
val    = polsx(j1)*polsy(j2)

dscaley = 1.0d0/(n-1)
if (j2 .gt. 1 .AND. j2 .lt. n) dscaley=dscaley*2

dscalex = 1.0d0/(n-1)
if (j1 .gt. 1 .AND. j1 .lt. n) dscalex=dscalex*2

! u(j1 + nx*(j2-1),i1+(i2-1)*nx) = val * dscaley * dscalex
u(idx,i1+(i2-1)*n) = val * dscaley * dscalex
idx = idx + 1

end do
end do

end do
end do


end subroutine


subroutine tensor_coefs(n,u,vals,coefs)
implicit double precision (a-h,o-z)

integer          ::  n
double precision ::  vals(n*n),coefs(n*(n+1)/2),u(n*(n+1)/2,n*n)

!
!  This is a simple utility code which applies the matrix u returned by
!  tensor_umatrix to a vector.  Its raison d'etre is that it allows
!  the user to shape the input vector vals as an (n,n) matix.  That is,
!  this routine is used to reshape the input vector vals.
!
!  Input parameters:
!    n - the number of points in the Chebyshev grid
!    u - the ((n+1)*n/2,n*n) matrix returned by tensor_umatrix  
!    vals - a (n,n) matrix containing the values of the expansion
!
!  Output parameters:
!    coefs- the coefficients in the expansion
!

coefs = matmul(u,vals)

end subroutine




subroutine tensor_eval(n,a,b,c,d,coefs,x,y,val)
implicit double precision (a-h,o-z)

integer, intent(in)           :: n
double precision, intent(in)  :: a,b,c,d,x,y,coefs(n*(n+1)/2)
double precision, intent(out) :: val

!
!  Evaluate an expansion of the form (1) at a specified point.
!
!  Input parameters:
!
!    n - the number of points in the Chebyshev grid on [-1,1]
!    (a,b) - the interval over which x varies
!    (c,d) - the interval over which y varies
!    coefs - the array of coefficients in the expansion
!    (x,y) - the point at which to evaluate the expansion
!
!  Output parameters:
!
!   val - the value of the expansion at the desired point
!

double precision :: polsx(n),polsy(n)

xx = (x-(b+a)/2)*2/(b-a)
yy = (y-(d+c)/2)*2/(d-c)

call chebs(xx,n-1,polsx)
call chebs(yy,n-1,polsy)

val  = 0
idx  = 1
do j=1,n
do i=1,n-j+1
val  = val  + coefs(idx) * polsx(i)*polsy(j)
idx  = idx+1
end do
end do

end subroutine


subroutine tensor_evalder(n,a,b,c,d,coefs,x,y,val,derx,dery)
implicit double precision (a-h,o-z)
integer, intent(in)           :: n
double precision, intent(in)  :: a,b,c,d,x,y,coefs(n*(n+1)/2)
double precision, intent(out) :: val
!
!  Evaluate an expansion of the form (1) at a specified point as well
!  as its derivative.
!
!  Input parameters:
!
!    n - the number of points in the Chebyshev grid on [-1,1]
!    (a,b) - the interval over which x varies
!    (c,d) - the interval over which y varies
!    coefs - the array of coefficients in the expansion
!    (x,y) - the point at which to evaluate the expansion
!
!  Output parameters:
!
!   val - the value of the expansion at the desired point
!   derx - the value of the derivative of the expansion w.r.t. x
!   dery - the value of the derivative of the expansion w.r.t. y
!

double precision :: polsx(n+10),polsy(n+10),dersx(n+10),dersy(n+10)

xx = (x-(b+a)/2)*2/(b-a)
yy = (y-(d+c)/2)*2/(d-c)

call chebders(xx,n-1,polsx,dersx)
call chebders(yy,n-1,polsy,dersy)

val   = 0
derx  = 0 
dery  = 0
idx   = 1
do j=1,n
do i=1,n-j+1
dd    = coefs(idx)
idx   = idx+1
val   = val   + dd * polsx(i)*polsy(j)
derx  = derx  + dd * dersx(i)*polsy(j)
dery  = dery  + dd * polsx(i)*dersy(j)
end do
end do

derx = derx * 2/(b-a)
dery = dery * 2/(d-c)

end subroutine


subroutine tensorpw_eval(n,nintsab,ab,nintscd,cd,coefs,x,y,val)
implicit double precision (a-h,o-z)

integer, intent(in)           :: n, nintsab, nintscd
double precision, intent(in)  :: ab(2,nintsab), cd(2,nintscd)
double precision, intent(in)  :: x,y,coefs(n*(n+1)/2,nintsab,nintscd)
double precision, intent(out) :: val
!
!  Evaluate a  piecewise bivariate Chebyshev expansion on a collection of rectangles 
!  of the form (2) at a specified point.
!
!  Input parameters:
!
!    n - the number of points in the Chebyshev grid, so that the order of the
!       expansion is n-1
!    ab - a (2,nintsab) array specifying the intervals (a_i,b_i) in (2)
!    cd - a (2,nintscd) array specifying the intervals (c_j,d_j) in (2)
!    coefs - the (ncoefs,nintsab,nintscd) array whose entries (:,i,j)
!     specify the coefficients for the rectangle (a_i,b_i) x (c_j, d_j)
!    (x,y) - the point at which to evaluate the expansion
!
!  Output parameters:
!   val - the value of the expansion at the specified point
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

val  = 0
idx  = 1
do j=1,n
do i=1,n-j+1
val  = val  + coefs(idx,intab,intcd) * polsx(i)*polsy(j)
idx  = idx+1
end do
end do

end subroutine



subroutine tensorpw_evalder(n,nintsab,ab,nintscd,cd,coefs,x,y,val,derx,dery)
implicit double precision (a-h,o-z)

integer, intent(in)           :: n, nintsab, nintscd
double precision, intent(in)  :: ab(2,nintsab), cd(2,nintscd)
double precision, intent(in)  :: x,y,coefs(n*(n+1)/2,nintsab,nintscd)
double precision, intent(out) :: val
!
!  Evaluate a  piecewise bivariate Chebyshev expansion on a collection of rectangles 
!  of the form (2) and its first order derivatives at a specified point.
!
!  Input parameters:
!
!    n - the number of points in the Chebyshev grid, so that the order of the
!       expansion is n-1
!    ab - a (2,nintsab) array specifying the intervals (a_i,b_i) in (2)
!    cd - a (2,nintscd) array specifying the intervals (c_j,d_j) in (2)
!    coefs - the (ncoefs,nintsab,nintscd) array whose entries (:,i,j)
!     specify the coefficients for the rectangle (a_i,b_i) x (c_j, d_j)
!    (x,y) - the point at which to evaluate the expansion
!
!  Output parameters:
!   val - the value of the expansion at the specified point
!   derx - the value of the deriv. w.r.t x of the expansion at the specified point
!   dery - the value of the deriv. w.r.t y of the expansion at the specified point

!

double precision :: polsx(n+10),polsy(n+10),dersx(n+10),dersy(n+10)


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

call chebders(xx,n-1,polsx,dersx)
call chebders(yy,n-1,polsy,dersy)

val   = 0
derx  = 0
dery  = 0

idx  = 1
do j=1,n
do i=1,n-j+1
dd    = coefs(idx,intab,intcd)
val   = val   + dd * polsx(i)*polsy(j)
derx  = derx  + dd * dersx(i)*polsy(j)
dery  = dery  + dd * polsx(i)*dersy(j)
idx  = idx+1
end do
end do

derx = derx * 2/(b-a)
dery = dery * 2/(d-c)


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

val1  = val1  + coefs1(idx1,intab,intcd) * polsx(i)*polsy(j)
val2  = val2  + coefs2(idx2,intab,intcd) * polsx(i)*polsy(j)
val3  = val3  + coefs3(idx3,intab,intcd) * polsx(i)*polsy(j)
val4  = val4  + coefs4(idx4,intab,intcd) * polsx(i)*polsy(j)

idx   = idx+1
end do
end do

end subroutine





end module
