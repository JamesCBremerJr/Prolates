!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains a simple codes for adaptive integration.
!
!  The following subroutines are user-callable:
!
!    adapint - adaptively compute the integral of a user-specified univariate function
!      on an interval
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module adapquad

use utils
use iso_c_binding


implicit double precision (a-h,o-z)

interface     

subroutine adapfun1d(n,xs,vals,userptr)
import c_ptr
implicit double precision (a-h,o-z)
double precision :: xs(n)
double precision :: vals(n)
type(c_ptr)      :: userptr
end subroutine




end interface 

contains



subroutine adapint(ier,eps,x1_0,x2_0,fun,userptr,val,ntotal)
implicit double precision (a-h,o-z)
procedure(adapfun1d)      :: fun
type(c_ptr)               :: userptr
double precision          :: eps,x1,x2
double precision          :: val

!
!  Adaptively compute the integral of a user-specified function on
!  the interval (x1,x2).
!
!  Input parameters:
!    (x1,y2) - the interval on which to evaluate the function
!    fun - the user-specified external function, which must conform to
!      the adapfun1d interface
!    userptr - a user-supplied pointer which is passed to the 
!       external function
!
!  Output parameters:
!     val - the value of the integral
!     ntotal - the number of nodes in the adaptively determined quadrature
!       rule used to evaluate the integral
!

double precision, allocatable :: stack(:,:)
double precision, allocatable :: vals(:)
double precision              :: sum1,sum2,sum0,sum00


! this is a 30 point Gauss-Legendre quadrature rule
double precision :: xs(30),whts(30)
integer          :: nquad
data nquad / 30 /
data xs  /                                     &
  -0.996893484074649540271630050918695241D+00, &
  -0.983668123279747209970032581605662842D+00, &
  -0.960021864968307512216871025581797646D+00, &
  -0.926200047429274325879324277080474015D+00, &
  -0.882560535792052681543116462530225587D+00, &
  -0.829565762382768397442898119732501916D+00, &
  -0.767777432104826194917977340974503127D+00, &
  -0.697850494793315796932292388026640086D+00, &
  -0.620526182989242861140477556431189243D+00, &
  -0.536624148142019899264169793311072772D+00, &
  -0.447033769538089176780609900322854008D+00, &
  -0.352704725530878113471037207089373868D+00, &
  -0.254636926167889846439805129817805114D+00, &
  -0.153869913608583546963794672743255924D+00, &
  -0.514718425553176958330252131667225730D-01, &
   0.514718425553176958330252131667225730D-01, &
   0.153869913608583546963794672743255924D+00, &
   0.254636926167889846439805129817805114D+00, &
   0.352704725530878113471037207089373868D+00, &
   0.447033769538089176780609900322854008D+00, &
   0.536624148142019899264169793311072772D+00, &
   0.620526182989242861140477556431189243D+00, &
   0.697850494793315796932292388026640086D+00, &
   0.767777432104826194917977340974503127D+00, &
   0.829565762382768397442898119732501916D+00, &
   0.882560535792052681543116462530225587D+00, &
   0.926200047429274325879324277080474015D+00, &
   0.960021864968307512216871025581797646D+00, &
   0.983668123279747209970032581605662842D+00, &
   0.996893484074649540271630050918695241D+00  /
data whts /                                    &
   0.796819249616660561546588347467375358D-02, &
   0.184664683110909591423021319120472140D-01, &
   0.287847078833233693497191796112920951D-01, &
   0.387991925696270495968019364463476587D-01, &
   0.484026728305940529029381404228076120D-01, &
   0.574931562176190664817216894020561829D-01, &
   0.659742298821804951281285151159623159D-01, &
   0.737559747377052062682438500221908157D-01, &
   0.807558952294202153546949384605296923D-01, &
   0.868997872010829798023875307151257485D-01, &
   0.921225222377861287176327070876187878D-01, &
   0.963687371746442596394686263518098668D-01, &
   0.995934205867952670627802821035694149D-01, &
   0.101762389748405504596428952168553985D+00, &
   0.102852652893558840341285636705414965D+00, &
   0.102852652893558840341285636705414965D+00, &
   0.101762389748405504596428952168553985D+00, &
   0.995934205867952670627802821035694149D-01, &
   0.963687371746442596394686263518098668D-01, &
   0.921225222377861287176327070876187878D-01, &
   0.868997872010829798023875307151257485D-01, &
   0.807558952294202153546949384605296923D-01, &
   0.737559747377052062682438500221908157D-01, &
   0.659742298821804951281285151159623159D-01, &
   0.574931562176190664817216894020561829D-01, &
   0.484026728305940529029381404228076120D-01, &
   0.387991925696270495968019364463476587D-01, &
   0.287847078833233693497191796112920951D-01, &
   0.184664683110909591423021319120472140D-01, &
   0.796819249616660561546588347467375358D-02  /

ier = 0


!
!  Allocate and initialize the stack
!

maxstack = 10000
allocate(stack(4,maxstack),vals(nquad))

istack          = 1
stack(1,istack) = x1_0
stack(2,istack) = x2_0

call fun(nquad,xs*(x2_0-x1_0)/2+(x2_0+x1_0)/2,vals,userptr)
sum1 = sum(vals*whts*(x2_0-x1_0)/2)
stack(3,istack) = sum1


val    = 0 
ntotal = 0

do while (istack .gt. 0)

!
!  Pop an interval off the stack
!

x1      = stack(1,istack)
x2      = stack(2,istack)
sum0    = stack(3,istack)
istack  = istack -1



!
!  Check to see if subdividing the interval yields a different result
!
x3 = (x1+x2)/2

xx1 = x1
xx2 = x3
call fun(nquad,xs*(xx2-xx1)/2+(xx2+xx1)/2,vals,userptr)
sum1 = sum(vals*whts*(xx2-xx1)/2)

xx1 = x3
xx2 = x2
call fun(nquad,xs*(xx2-xx1)/2+(xx2+xx1)/2,vals,userptr)
sum2 = sum(vals*whts*(xx2-xx1)/2)


sum00 = sum1+sum2
diff  = abs(sum0-sum00)


!
!  If not, add the contribution of the interval to the total
!


if (diff .lt. eps) then
val    = val + sum0
ntotal = ntotal + nquad
else

!
!  Otherwise, push the child intervals onto the stack
!

if (istack+2 .gt. maxstack) then
ier = 4
return
endif

istack = istack+1
stack(1,istack) = x1
stack(2,istack) = x3
stack(3,istack) = sum1

istack = istack+1
stack(1,istack) = x3
stack(2,istack) = x2
stack(3,istack) = sum2


end if

end do

end subroutine




end module
