!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  These are utility routines.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module experiments2_functions

use utils
use prolates
use prolates_phase
use prolates_exps

implicit double precision (a-h,o-z)

integer, parameter                :: ntimes = 1               ! number of times to repeat timings
integer, parameter                :: ncount = 100             ! number of samples per sigma / gamma
integer, parameter                :: nx     = 100             ! number of samples in the argument

integer, parameter                :: ntimes2 = 10             ! number of times to repeat timings
integer, parameter                :: ncount2 = 20             ! number of samples per sigma / gamma
integer, parameter                :: nx2     = 100            ! number of samples in the argument

contains

subroutine build_phase(chebdata,exp3,phase,gamma,sigma)
implicit double precision (a-h,o-z)
type(chebexps_data)                         :: chebdata
type(exp3_data)                             :: exp3
type(prolates_phase_data), intent(out)      :: phase
!
!
!
data pi / 3.14159265358979323846264338327950288d0 / 

call prolates_exp3_eval(exp3,gamma,sigma,chi,apval,apppval)

aval          = -pi/2*(xi+1)
appval        = -apval

phase%gamma   = gamma
phase%chi     = chi
phase%aval    = -pi/2*(xi+1)
phase%apval   = apval
phase%appval  = appval
phase%apppval = apppval

! Build the phase function

call prolates_phase2(chebdata,phase)
end subroutine


subroutine build_phase2(chebdata,exp3,phase,gamma,n)
implicit double precision (a-h,o-z)
type(chebexps_data)                         :: chebdata
type(exp3_data)                             :: exp3
type(prolates_phase_data), intent(out)      :: phase
!
!
!
double precision, allocatable               :: coefsps(:)
data pi / 3.14159265358979323846264338327950288d0 / 


!xi    = n
!sigma = xi/gamma


! Build the Legendre expansion
inormal = 1
call prolates_integer(inormal,gamma,n,chi,n1,n2,coefsps)

! call prolates_exp3_eval(exp3,gamma,sigma,chi,apval,apppval)

! aval          = -pi/2*(xi+1)
! appval        = -apval

! phase%gamma   = gamma
! phase%chi     = chi
! phase%aval    = -pi/2*(xi+1)
! phase%apval   = apval
! phase%appval  = appval
! phase%apppval = apppval

! Build the phase function

call prolates_phase1(chebdata,gamma,chi,phase)
call prolates_phase2(chebdata,phase)

end subroutine

subroutine phase_cost(phase,ncost)
implicit double precision (a-h,o-z)
type(prolates_phase_data)      :: phase
!
!  Return the number of coefficients actually required to represent the 
!  phase function 
!
double precision, allocatable :: coefs(:)

! Figure out how many coefficients are actually required to represent it
k     = phase%k
nints = phase%nints
! ncost = k*nints
! return

nn    = 0
allocate(coefs(k))

do int=1,nints
coefs = phase%apcoefs(:,int)
coefs = abs(coefs)
coefs = coefs / coefs(1)
do i=1,k
if (coefs(i) .gt. 1.0d-14) nn = max(nn,i)
end do
end do

ncost          = nn*nints
end subroutine


subroutine ps_accuracy(chebdata,gamma,n,derr)
implicit double precision (a-h,o-z)
!
!  Check the accuracy of PS_n(z;\gamma) at 100 equispaced points
!  on the intervla (0,1)
!
type(chebexps_data)                         :: chebdata
type(prolates_phase_data)                   :: phase
double precision, allocatable               :: coefsps(:)
double precision, allocatable               :: derrs(:)


! Build the Legendre expansion
inormal = 1
call prolates_integer(inormal,gamma,n,chi0,n1,n2,coefsps)
call prolates_phase1(chebdata,gamma,chi0,phase)

! Compute the normalization constant
call prolates_phase_eval0(phase,valps,derps)
call prolates_legen0(n,val0,der0)
if (mod(n,2) .eq. 1) then
dnorm = der0/derps
else
dnorm = val0/valps
endif

x1 = 0.0d0
x2 = 1.0d0
allocate(derrs(nx))

do i=1,nx
dd = (i+0.0d0)/(nx+1.0d0)
x  = x1 + (x2-x1)*dd

call prolates_integer_eval(n,n1,n2,coefsps,x,valps0)
call prolates_phase_eval(phase,x,valps)
valps       = valps*dnorm
!derr        =  abs(valps-valps0)/(abs(valps0)+1)
derr        =  abs(valps-valps0)

derrs(i)    = derr
end do

derr = maxval(derrs)

end subroutine


subroutine ps_time1(chebdata,exp3,gamma,n,tphase1,tphase2,trokhlin)
implicit double precision (a-h,o-z)
!
!  Return the time required to construct the phase function with the
!  accelerated and unaccelerated algorithm and the time to construct
!  the Legendre expansions with Osipov-Xiao-Rokhlin.
!
type(chebexps_data)                         :: chebdata
type(exp3_data)                             :: exp3
type(prolates_phase_data)                   :: phase
double precision, allocatable               :: coefsps(:)
double precision, allocatable               :: derrs(:)


! Build the Legendre expansion
inormal = 1
call elapsed(t1)
do ii=1,ntimes
call prolates_integer(inormal,gamma,n,chi0,n1,n2,coefsps)
end do
call elapsed(t2)
trokhlin = (t2-t1) / (ntimes+0.0d0)


call elapsed(t1)
do ii=1,ntimes
call prolates_phase1(chebdata,gamma,chi0,phase)
end do
call elapsed(t2)
tphase1 = (t2-t1)/(ntimes+0.0d0)


call elapsed(t1)
do i=1,ntimes
xi = n
sigma = xi/gamma
call prolates_exp3_eval(exp3,gamma,sigma,chi,apval,apppval)
aval          = -pi/2*(xi+1)
appval        = -apval
phase%gamma   = gamma
phase%chi     = chi
phase%aval    = -pi/2*(xi+1)
phase%apval   = apval
phase%appval  = appval
phase%apppval = apppval
call prolates_phase2(chebdata,phase)
end do
call elapsed(t2)

tphase2=(t2-t1)/(ntimes+0.0d0)

end subroutine


subroutine ps_time2(chebdata,exp3,gamma,n,tphase1,tphase2)
implicit double precision (a-h,o-z)
!
!  Return the time required to construct the phase function using the
!  accelerated and unaccelerated algorithms.
!
type(chebexps_data)                         :: chebdata
type(exp3_data)                             :: exp3
type(prolates_phase_data)                   :: phase
double precision, allocatable               :: coefsps(:)
double precision, allocatable               :: derrs(:)


! ! Build the Legendre expansion
! inormal = 1
! call elapsed(t1)
! do ii=1,ntimes
! call prolates_integer(inormal,gamma,n,chi0,n1,n2,coefsps)
! end do
! call elapsed(t2)
! trokhlin = (t2-t1) / (ntimes+0.0d0)

! Phase function
call elapsed(t1)
do i=1,ntimes
xi = n
sigma = xi/gamma
call prolates_exp3_eval(exp3,gamma,sigma,chi,apval,apppval)
aval          = -pi/2*(xi+1)
appval        = -apval
phase%gamma   = gamma
phase%chi     = chi
phase%aval    = -pi/2*(xi+1)
phase%apval   = apval
phase%appval  = appval
phase%apppval = apppval
call prolates_phase2(chebdata,phase)
end do
call elapsed(t2)

tphase2=(t2-t1)/(ntimes+0.0d0)

call elapsed(t1)
do ii=1,ntimes
call prolates_phase1(chebdata,gamma,chi,phase)
end do
call elapsed(t2)
tphase1 = (t2-t1)/(ntimes+0.0d0)

end subroutine


subroutine ps_time3(chebdata,exp3,gamma,n,tphase)
implicit double precision (a-h,o-z)
!
!  Compare the time required to evaluate PS_n(z;\gamma) with the 
!  phase function
!
type(chebexps_data)                         :: chebdata
type(exp3_data)                             :: exp3
type(prolates_phase_data)                   :: phase


xi    = n
sigma = xi/gamma
x1    = 0.0d0
x2    = 1.0d0
call build_phase(chebdata,exp3,phase,gamma,sigma)

call elapsed(t1)
do ii=1,ntimes2
do i=1,nx2
dd = (i-1.0d0)/(nx2-0.0d0)
call random_number(dd)
x  = x1 + dd * (x2-x1)*dd
call prolates_phase_eval_ps(phase,x,valps)
end do
end do
call elapsed(t2)

call elapsed(t1)
do ii=1,ntimes2
do i=1,nx2
dd = (i-1.0d0)/(nx2-0.0d0)
x  = x1 + dd * (x2-x1)*dd
call prolates_phase_eval_ps(phase,x,valps)
end do
end do
call elapsed(t2)

tphase = (t2-t1) / (ntimes2+0.0d0) * 1.0d0/(nx2+0.0d0)

end subroutine


subroutine ps_time4(chebdata,exp3,gamma,n,tphase,trokhlin)
implicit double precision (a-h,o-z)
!
!  Return the time required to evaluate the phase
!

type(chebexps_data)                         :: chebdata
type(exp3_data)                             :: exp3
type(prolates_phase_data)                   :: phase
double precision, allocatable               :: coefsps(:)


! Build the Legendre expansion
inormal = 1
call prolates_integer(inormal,gamma,n,chi0,n1,n2,coefsps)


x1 = 0.0d0
x2 = 1.0d0

call elapsed(t1)
do ii=1,ntimes2
do i=1,nx2
dd = (i-1.0d0)/(nx2-0.0d0)
x  = x1 + dd * (x2-x1)*dd
call prolates_integer_eval(n,n1,n2,coefsps,x,valps)
end do
end do
call elapsed(t2)
trokhlin = (t2-t1) / (ntimes2+0.0d0) * 1.0d0/(nx2+0.0d0)



call prolates_phase1(chebdata,gamma,chi0,phase)
call elapsed(t1)
do ii=1,ntimes2
do i=1,nx2
dd = (i-1.0d0)/(nx2-0.0d0)
x  = x1 + dd * (x2-x1)*dd
call prolates_phase_eval_ps(phase,x,valps)
end do
end do
call elapsed(t2)
tphase = (t2-t1) / (ntimes2+0.0d0) * 1.0d0/(nx2+0.0d0)

end subroutine


end module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  The code for executing the experiments.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program experiments2

use utils
use prolates
use prolates_phase
use prolates_exps
use experiments2_functions

implicit double precision (a-h,o-z)

type(chebexps_data)            :: chebdata
type(prolates_phase_data)      :: phase
type(exp3_data)                :: exp3
double precision, allocatable  :: coefsps(:), coefs(:)
double precision, allocatable  :: xs(:), ys(:,:), gammas(:), sigmas(:)
double precision, allocatable  :: derrs1(:,:), derrs2(:,:)
double precision, allocatable  :: times1(:,:), times2(:,:), times3(:,:)
integer, allocatable           :: ns(:), ncosts1(:,:), ncosts2(:,:)

character(len=100)               :: label1, label2, label3, label4, label5
character(len=100)               :: label6, label7, label8, label9, label10
character(len=100), allocatable  :: range_labels(:)


pi     = acos(-1.0d0)
eps0   = epsilon(0.0d0)
k      = 30
call chebexps(k,chebdata)

iw = 20
open(iw,FILE="expansions.txt",STATUS='old')
call elapsed(t1)
call prolates_exp3_read(iw,exp3)
call elapsed(t2)
call prin2("expansions.txt read time = ",t2-t1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Calculate the cost of storing chi only
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nexps = exp3%nints
ncost = k + 2*nexps

do i=1,nexps
  do j=1,k
    nints = exp3%exps(i)%exps(j)%nints
    kk     = exp3%exps(i)%exps(j)%k
   allocate(coefs(kk))
    nn    = 0 
    do int=1,nints
      coefs = exp3%exps(i)%exps(j)%chi_coefs(:,int)
      coefs = abs(coefs)
      coefs = coefs / coefs(1)
      do l=1,kk
        if (coefs(l) .gt. 1.0d-12) nn = max(nn,l)      
      end do
    end do
    ncost = ncost + nn*nints
    deallocate(coefs)
  end do
end do
dd = ncost*8.0d0/(1024.0d0*1024.0d0)
call prin2("cost of chi expansion (in MB) = ",dd)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Cost to represent phase function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! print *,"-------- COST OF REPRESENTATIONS ------------ "

! gamma1 =     100.0d0
! gamma2 = 1000000.0d0
! iflogx = 1 
! nn     = 100

! ylimit1 = 0
! ylimit2 = 800

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the cost to represent the phase function as a function of sigma for several small 
! !  values of gamma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call elapsed(t1)
! nfuns = 4
! allocate(xs(nn), ys(nfuns,nn), gammas(nfuns) )

! gammas(1) = 100
! gammas(2) = 200
! gammas(3) = 500
! gammas(4) = 1000

! sigma1 = 0.0d0
! sigma2 = 1.0d0

! do i=1,nn

! dd    = (i-1.0d0)/(nn-1.0d0)
! sigma = sigma1 + (sigma2-sigma1)*dd

! !!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!

! do j=1,nfuns
! gamma = gammas(j)
! n     = gamma*sigma
! xi    = n
! sigma = xi/gamma
! call build_phase(chebdata,exp3,phase,gamma,sigma)
! call phase_cost(phase,ncost)
! ys(j,i) = ncost
! end do

! xs(i)   = sigma
! end do

! call plot_functions("plot101.py","plot101.pdf",nfuns,"$\sigma$","Number of Chebyshev Coefficients",0,0,   &
!    sigma1,sigma2,ylimit1,ylimit2,nn,xs,ys,"upper right",&
!    "$\gamma=100$*$\gamma=200$*$\gamma=500$*$\gamma=1000$*", &
!    "b-*r--*g-.*k:*mx*")


! deallocate(xs,ys,gammas)
! call elapsed(t2)
! print *,"plot 101: ",t2-t1


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the cost to represent the phase function on the interval as 
! !  a function of sigma for a range of values of gamma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call elapsed(t1)
! nfuns = 4
! allocate(xs(nn), ys(nfuns,nn), gammas(nfuns) )

! gammas(1) = 1000
! gammas(2) = 10000
! gammas(3) = 100000
! gammas(4) = 1000000

! sigma1 = 0.0d0
! sigma2 = 1.0d0

! do i=1,nn

! dd    = (i-1.0d0)/(nn-1.0d0)
! sigma = sigma1 + (sigma2-sigma1)*dd

! do j=1,nfuns
! gamma = gammas(j)
! n     = gamma*sigma
! xi    = n
! sigma = xi/gamma
! call build_phase(chebdata,exp3,phase,gamma,sigma)

! call phase_cost(phase,ncost)
! ys(j,i) = ncost


! end do

! xs(i)   = sigma

! end do

! call plot_functions("plot102.py","plot102.pdf",nfuns,"$\sigma$","Number of Chebyshev Coefficients",0,0,   &
!    sigma1,sigma2,ylimit1,ylimit2,nn,xs,ys,"upper right",&
!    "$\gamma=10^3$*$\gamma=10^4$*$\gamma=10^5$*$\gamma=10^6$*", &
!    "b-*r--*g-.*k:*mx*")


! deallocate(xs,ys,gammas)
! call elapsed(t2)
! print *,"plot 102: ",t2-t1


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the cost to represent the phase function on the interval as 
! !  a function of sigma for large gamma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call elapsed(t1)
! nfuns = 4
! allocate(xs(nn), ys(nfuns,nn), gammas(nfuns) )

! gammas(1) = 0.60d0*10**6
! gammas(2) = 0.70d0*10**6
! gammas(3) = 0.80d0*10**6
! gammas(4) = 1.00d0*10**6

! sigma1 = 0.0d0
! sigma2 = 1.0d0

! do i=1,nn

! dd    = (i-1.0d0)/(nn-1.0d0)
! sigma = sigma1 + (sigma2-sigma1)*dd

! do j=1,nfuns
! gamma = gammas(j)
! n     = gamma*sigma
! xi    = n
! sigma = xi/gamma
! call build_phase(chebdata,exp3,phase,gamma,sigma)
! call phase_cost(phase,ncost)
! ys(j,i) = ncost

! end do

! xs(i)   = sigma

! end do

! call plot_functions("plot103.py","plot103.pdf",nfuns,"$\sigma$","Number of Chebyshev Coefficients",0,0,   &
!    sigma1,sigma2,ylimit1,ylimit2,nn,xs,ys,"upper right",&
!    "$\gamma=0.7 \\times 10^6$*$\gamma=0.8 \\times 10^6$*$\gamma=0.9 \\times 10^6$*$\gamma = 1.0 \\times 10^6$*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,gammas)

! call elapsed(t2)
! print *,"plot 103: ",t2-t1


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the cost to represent the phase function as a function of gamma for several 
! !  values of sigma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call elapsed(t1)
! nfuns = 4
! allocate(xs(nn), ys(nfuns,nn), sigmas(nfuns) )

! gamma1   = 100.0d0
! gamma2   = 1000000.0d0

! sigmas(1) = 0.10d0
! sigmas(2) = 0.25d0
! sigmas(3) = 0.50d0
! sigmas(4) = 0.60d0

! do i=1,nn

! dd    = (i-1.0d0)/(nn-1.0d0)
! if (iflogx .eq. 0) then
! gamma = gamma1  + (gamma2-gamma1)*dd
! else
! dd2   = log(gamma1) + dd*(log(gamma2)-log(gamma1))
! gamma = exp(dd2)
! endif

! do j=1,nfuns
! sigma = sigmas(j)
! n     = gamma*sigma
! xi    = n
! sigma = xi/gamma
! call build_phase(chebdata,exp3,phase,gamma,sigma)
! call phase_cost(phase,ncost)
! ys(j,i) = ncost
! end do

! xs(i)   = gamma

! end do

! call plot_functions("plot104.py","plot104.pdf",nfuns,"$\gamma$","Number of Chebyshev Coefficients",iflogx,0,   &
!    gamma1,gamma2,ylimit1,ylimit2,nn,xs,ys,"lower right",&
!    "$\sigma=0.10$*$\sigma=0.25$*$\sigma=0.50$*$\sigma=0.60$*", &
!    "b-*r--*g-.*k:*mx*")


! deallocate(xs,ys,sigmas)
! call elapsed(t2)
! print *,"plot 104:", t2-t1


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the cost to represent the phase function as a function of gamma for several 
! !  values of sigma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! call elapsed(t1)
! nfuns = 4
! allocate(xs(nn), ys(nfuns,nn), sigmas(nfuns) )

! gamma1   = 100.0d0
! gamma2   = 1000000.0d0

! sigmas(1) = 0.62d0
! sigmas(2) = 0.63d0
! sigmas(3) = 2/pi
! sigmas(4) = 0.66d0

! do i=1,nn

! dd    = (i-1.0d0)/(nn-1.0d0)
! if (iflogx .eq. 0) then
! gamma = gamma1  + (gamma2-gamma1)*dd
! else
! dd2   = log(gamma1) + dd*(log(gamma2)-log(gamma1))
! gamma = exp(dd2)
! endif

! do j=1,nfuns
! sigma = sigmas(j)
! n     = gamma*sigma
! xi    = n
! sigma = xi/gamma
! call build_phase(chebdata,exp3,phase,gamma,sigma)
! call phase_cost(phase,ncost)
! ys(j,i) = ncost
! end do

! xs(i)   = gamma

! end do

! call plot_functions("plot105.py","plot105.pdf",nfuns,"$\gamma$","Number of Chebyshev Coefficients",iflogx,0,   &
!    gamma1,gamma2,ylimit1,ylimit2,nn,xs,ys,"upper left",&
!    "$\sigma=0.62$*$\sigma=0.63$*$\sigma=\\frac{2}{\pi}$*$\sigma=0.66$*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,sigmas)
! call elapsed(t2)
! print *,"plot 105: ",t2-t1


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the cost to represent the phase function as a function of gamma for several 
! !  values of sigma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! call elapsed(t1)
! nfuns = 4
! allocate(xs(nn), ys(nfuns,nn), sigmas(nfuns) )

! gamma1   = 100.0d0
! gamma2   = 1000000.0d0

! sigmas(1) = 0.70d0
! sigmas(2) = 0.80d0
! sigmas(3) = 0.90d0
! sigmas(4) = 1.00d0

! do i=1,nn

! dd    = (i-1.0d0)/(nn-1.0d0)
! if (iflogx .eq. 0) then
! gamma = gamma1  + (gamma2-gamma1)*dd
! else
! dd2   = log(gamma1) + dd*(log(gamma2)-log(gamma1))
! gamma = exp(dd2)
! endif

! ! dd    = (i-1.0d0)/(nn-1.0d0)
! ! gamma = gamma1  + (gamma2-gamma1)*dd

! do j=1,nfuns
! sigma = sigmas(j)
! n     = gamma*sigma
! xi    = n
! sigma = xi/gamma
! call build_phase(chebdata,exp3,phase,gamma,sigma)
! call phase_cost(phase,ncost)
! ys(j,i) = ncost
! end do

! xs(i)   = gamma


! end do

! call plot_functions("plot106.py","plot106.pdf",nfuns,"$\gamma$","Number of Chebyshev Coefficients",iflogx,0,   &
!    gamma1,gamma2,ylimit1,ylimit2,nn,xs,ys,"upper left",&
!    "$\sigma=0.70$*$\sigma=0.80$*$\sigma=0.9$*$\sigma=1.0$*", &
!    "b-*r--*g-.*k:*mx*")


! deallocate(xs,ys,sigmas)
! call elapsed(t2)
! print *,"plot 106: ",t2-t1

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Create a table listing the cost of representing PsiP
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call elapsed(t3)
! nsigma = 4
! ngamma = 8

! allocate(sigmas(nsigma+1),gammas(ngamma+1))
! allocate(range_labels(ngamma))

! sigmas(1)  = 0.00d0
! sigmas(2)  = 0.25d0
! sigmas(3)  = 0.50d0
! sigmas(4)  = 0.75d0
! sigmas(5)  = 1.00d0

! gammas(1)  = 100
! gammas(2)  = 500
! gammas(3)  = 1000
! gammas(4)  = 5000
! gammas(5)  = 10000
! gammas(6)  = 50000
! gammas(7)  = 100000
! gammas(8)  = 500000
! gammas(9)  = 1000000

! label1     = "100 to 500"
! label2     = "500 to 1,000"
! label3     = "1,000 to 5,000"
! label4     = "5,000 to 10,000"
! label5     = "10,000 to 50,000"
! label6     = "50,000 to 100,000"
! label7     = "100,000 to 500,000"
! label8     = "500,000 to 1,000,000"

! range_labels(1) = label1
! range_labels(2) = label2
! range_labels(3) = label3
! range_labels(4) = label4
! range_labels(5) = label5
! range_labels(6) = label6
! range_labels(7) = label7
! range_labels(8) = label8

! allocate(ncosts1(ncount,ncount))
! allocate(ncosts2(ncount,ncount))

! iw = 1001
! open(iw,FILE="table101.tex")
! write(iw,"(A)") "\begin{tabular}{ccr@{\hspace{2em}}ccr}"
! write(iw,"(A)") "\toprule"
! write(iw,"(A)",advance="no") "Range of $\gamma$ & Range of $\sigma$ & Max Coefs  &"
! write(iw,"(A)",advance="no") "Range of $\gamma$ & Range of $\sigma$ & Max Coefs \\"
! write(iw,"(A)") "\midrule"

! do i=1,ceiling(ngamma/2.0d0)
! do j=1,nsigma

! sigma1 = sigmas(j)
! sigma2 = sigmas(j+1)

! ! Find the first error
! gamma1 = gammas(i)
! gamma2 = gammas(i+1)


! !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,jj,dd1,dd2,sigma,gamma,phase,ncost)
! !$OMP DO
! do ii=1,ncount
! do jj=1,ncount

! dd1   = (ii-1.0d0)/(ncount-1.0d0)
! dd2   = (jj-1.0d0)/(ncount-1.0d0)
! sigma = sigma1 + (sigma2-sigma1)*dd1
! gamma = gamma1 + (gamma2-gamma1)*dd2

! n     = sigma*gamma
! xi    = n 
! sigma = xi / gamma
! call build_phase(chebdata,exp3,phase,gamma,sigma)
! call phase_cost(phase,ncost)

! ncosts1(ii,jj) = ncost

! end do
! end do
! !$OMP END DO
! !$OMP END PARALLEL

! ncost1 = maxval(ncosts1)

! if (j .eq. 1) then
! write(iw,"(A,A)",advance = 'no') range_labels(i)," &"
! else
! write(iw,"(A)") " & "
! endif


! write(iw,"(A,F4.2,A,F4.2,A)", advance='no')  "$",sigmas(j)," - ",sigmas(j+1), " $ & "
! write(iw,"(A,I5,A)",advance='no') "$",ncost1,"$ &"

! i2 = ceiling(ngamma/2.0d0)+i
! if (i2 .gt. ngamma) then
! ncost2 = -1
! else
! ! Find the second error
! gamma1 = gammas(i2)
! gamma2 = gammas(i2+1)


! !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,jj,dd1,dd2,sigma,gamma,phase,ncost)
! !$OMP DO
! do ii=1,ncount
! do jj=1,ncount

! dd1   = (ii-1.0d0)/(ncount-1.0d0)
! dd2   = (jj-1.0d0)/(ncount-1.0d0)
! sigma = sigma1 + (sigma2-sigma1)*dd1
! gamma = gamma1 + (gamma2-gamma1)*dd2

! call build_phase(chebdata,exp3,phase,gamma,sigma)
! call phase_cost(phase,ncost)

! ncosts2(ii,jj) = ncost


! end do
! end do
! !$OMP END DO
! !$OMP END PARALLEL

! ncost2 = maxval(ncosts2)
! endif



! if (ncost2 .eq. -1) then

! write(iw,"(A,I5.5,A)") "\\"

! else

! if (j .eq. 1) then
! nn1 = gamma1
! nn2 = gamma2
! write(iw,"(A,A)",advance = 'no') range_labels(i2)," &"

! else
! write(iw,"(A)") " & "
! endif

! write(iw,"(A,F4.2,A,F4.2,A)", advance='no')  "$",sigmas(j)," - ",sigmas(j+1), " $ & "
! write(iw,"(A,I5,A)") "$",ncost2,"$ \\"
! endif


! write(iw,"(A)")  "\addlinespace[.125em]"


! end do
! write(iw,"(A)")  "\addlinespace[.25em]"
! end do


! write(iw,"(A)") "\bottomrule"
! write(iw,"(A)") "\end{tabular}"

! close(iw)

! deallocate(gammas, sigmas, ncosts1, ncosts2, range_labels)

! call elapsed(t4)

! print *,"table101 : ",t4-t3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Accuracy of the phase function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! print *,"-------- ACCURACY OF REPRESENTATIONS ------------ "

! iflogx  = 1
! nn      = 100

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the accuracy of the evaluation of PS as a function of gamma for several sigma
! !  and smallish gamma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call elapsed(t1)
! nfuns   = 4
! allocate(xs(nn), ys(nfuns,nn), sigmas(nfuns) )

! gamma1    = 10.0d0**2
! gamma2    = 10.0d0**6

! sigmas(1) = 0.10d0
! sigmas(2) = 0.25d0
! sigmas(3) = 0.50d0
! sigmas(4) = 0.60d0

! do ifun=1,nfuns

! sigma = sigmas(ifun)

! !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,dd2,gamma,n,dmax)
! !$OMP DO 
! do i=1,nn

! dd    = (i-1.0d0)/(nn-1.0d0)
! if (iflogx .eq. 0) then
! gamma = gamma1  + (gamma2-gamma1)*dd
! else
! dd2   = log(gamma1) + dd*(log(gamma2)-log(gamma1))
! gamma = exp(dd2)
! endif
! n      = sigma*gamma
! call ps_accuracy(chebdata,gamma,n,dmax)

! xs(i)       = gamma
! ys(ifun,i)  = dmax


! end do
! !$OMP END DO
! !$OMP END PARALLEL

! end do

! call plot_functions("plot110.py","plot110.pdf",nfuns,"$\gamma$","Maximum Observed Error",iflogx,1,   &
!    gamma1,gamma2,1.0d-18,1.0d-10,nn,xs,ys,"lower right",&
!    "$\sigma=0.10$*$\sigma=0.25$*$\sigma=0.50$*$\sigma=0.60$*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,sigmas)
! call elapsed(t2)
! print *,"plot 110: ",t2-t1



! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the accuracy of the evaluation of PS as a function of gamma for several sigma
! !  and smallish gamma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call elapsed(t1)
! nfuns = 4
! allocate(xs(nn), ys(nfuns,nn), sigmas(nfuns) )

! gamma1 = 100.0d0
! gamma2 = 1000000.0d0

! sigmas(1) = 0.62d0
! sigmas(2) = 0.63d0
! sigmas(3) = 2/pi
! sigmas(4) = 0.66d0

! do ifun=1,nfuns

! sigma = sigmas(ifun)

! !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,dd2,gamma,n,dmax)
! !$OMP DO 
! do i=1,nn

! dd    = (i-1.0d0)/(nn-1.0d0)
! if (iflogx .eq. 0) then
! gamma = gamma1  + (gamma2-gamma1)*dd
! else
! dd2   = log(gamma1) + dd*(log(gamma2)-log(gamma1))
! gamma = exp(dd2)
! endif

! n      = sigma*gamma
! call ps_accuracy(chebdata,gamma,n,dmax)

! xs(i)       = gamma
! ys(ifun,i)  = dmax

! end do
! !$OMP END DO
! !$OMP END PARALLEL

! end do

! call plot_functions("plot111.py","plot111.pdf",nfuns,"$\gamma$","Maximum Observed Error",iflogx,1,   &
!    gamma1,gamma2,1.0d-18,1.0d-10,nn,xs,ys,"lower right",&
!    "$\sigma=0.62$*$\sigma=0.63$*$\sigma=\\frac{2}{\pi}$*$\sigma=0.66$*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,sigmas)
! call elapsed(t2)
! print *,"plot 111: ",t2-t1

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the accuracy of the evaluation of PS as a function of gamma for several sigma
! !  and smallish gamma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call elapsed(t1)
! nfuns = 4
! allocate(xs(nn), ys(nfuns,nn), sigmas(nfuns) )

! gamma1 = 100.0d0
! gamma2 = 1000000.0d0

! sigmas(1) = 0.70d0
! sigmas(2) = 0.80d0
! sigmas(3) = 0.90d0
! sigmas(4) = 1.00d0

! do ifun=1,nfuns

! sigma = sigmas(ifun)

! !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,dd2,gamma,n,dmax)
! !$OMP DO 
! do i=1,nn

! dd    = (i-1.0d0)/(nn-1.0d0)
! if (iflogx .eq. 0) then
! gamma = gamma1  + (gamma2-gamma1)*dd
! else
! dd2   = log(gamma1) + dd*(log(gamma2)-log(gamma1))
! gamma = exp(dd2)
! endif

! n      = sigma*gamma
! call ps_accuracy(chebdata,gamma,n,dmax)

! xs(i)       = gamma
! ys(ifun,i)  = dmax

! end do
! !$OMP END DO
! !$OMP END PARALLEL

! end do

! call plot_functions("plot112.py","plot112.pdf",nfuns,"$\gamma$","Maximum Observed Error",iflogx,1,   &
!    gamma1,gamma2,1.0d-18,1.0d-10,nn,xs,ys,"lower right",&
!    "$\sigma=0.70$*$\sigma=0.80$*$\sigma=0.90$*$\sigma=1.0$*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,sigmas)
! call elapsed(t2)
! print *,"plot 112: ",t2-t1


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the accuracy of the evaluation of PS as a function of sigma for several gamma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! call elapsed(t1)

! nfuns = 4
! allocate(xs(nn), ys(nfuns,nn), gammas(nfuns) )

! gammas(1) = 100
! gammas(2) = 200
! gammas(3) = 500
! gammas(4) = 1000

! sigma1 = 0.0d0
! sigma2 = 1.0d0

! do ifun=1,nfuns
! gamma = gammas(ifun)

! !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,sigma,xi,n,dmax)
! !$OMP DO 
! do i=1,nn

! dd     = (i-1.0d0)/(nn-1.0d0)
! sigma  = sigma1 + (sigma2-sigma1)*dd
! xi     = sigma*gamma
! n      = xi
! call ps_accuracy(chebdata,gamma,n,dmax)
! xs(i)      = sigma
! ys(ifun,i) = dmax

! end do
! !$OMP END DO
! !$OMP END PARALLEL
! end do

! call plot_functions("plot113.py","plot113.pdf",nfuns,"$\sigma$","Maximum Observed Error",0,1,   &
!    0.0d0,1.0d0,1.0d-18,1.0d-10,nn,xs,ys,"lower right",&
!    "$\gamma=100$*$\gamma=200$*$\gamma=500$*$\gamma=1000$*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,gammas)

! call elapsed(t2)
! print *,"plot 113: ",t2-t1


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the accuracy of the evaluation of PS as a function of sigma for several gamma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call elapsed(t1)

! nfuns = 4
! allocate(xs(nn), ys(nfuns,nn), gammas(nfuns) )

! gammas(1) = 10**3
! gammas(2) = 10**4
! gammas(3) = 10**5
! gammas(4) = 10**6

! sigma1 = 0.0d0
! sigma2 = 1.0d0

! do ifun=1,nfuns
! gamma = gammas(ifun)

! !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,sigma,xi,n,dmax)
! !$OMP DO 
! do i=1,nn

! dd     = (i-1.0d0)/(nn-1.0d0)
! sigma  = sigma1 + (sigma2-sigma1)*dd
! xi     = sigma*gamma
! n      = xi
! call ps_accuracy(chebdata,gamma,n,dmax)
! xs(i)      = sigma
! ys(ifun,i) = dmax

! end do
! !$OMP END DO
! !$OMP END PARALLEL
! end do

! call plot_functions("plot114.py","plot114.pdf",nfuns,"$\sigma$","Maximum Observed Error",0,1,   &
!    0.0d0,1.0d0,1.0d-18,1.0d-10,nn,xs,ys,"lower right",&
!    "$\gamma=10^3$*$\gamma=10^4$*$\gamma=10^5$*$\gamma=10^6$*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,gammas)

! call elapsed(t2)
! print *,"plot 114: ",t2-t1

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the accuracy of the evaluation of PS as a function of sigma for several gamma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call elapsed(t3)

! nfuns = 4
! allocate(xs(nn), ys(nfuns,nn), gammas(nfuns) )

! gammas(1) = 0.7d0 * 10**6
! gammas(2) = 0.8d0 * 10**6
! gammas(3) = 0.9d0 * 10**6
! gammas(4) = 1.0d0 * 10**6

! sigma1 = 0.0d0
! sigma2 = 1.0d0

! do ifun=1,nfuns
! gamma = gammas(ifun)

! !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,sigma,xi,n,dmax)
! !$OMP DO 
! do i=1,nn

! dd     = (i-1.0d0)/(nn-1.0d0)
! sigma  = sigma1 + (sigma2-sigma1)*dd
! xi     = sigma*gamma
! n      = xi
! call ps_accuracy(chebdata,gamma,n,dmax)
! xs(i)      = sigma
! ys(ifun,i) = dmax

! end do
! !$OMP END DO
! !$OMP END PARALLEL
! end do

! call plot_functions("plot115.py","plot115.pdf",nfuns,"$\sigma$","Maximum Observed Error",0,1,   &
!    0.0d0,1.0d0,1.0d-18,1.0d-10,nn,xs,ys,"lower right",&
!    "$\gamma=0.7 \\times 10^6$*$\gamma=0.8 \\times 10^6$*$\gamma=0.9 \\times 10^6$*$\gamma=10^6$*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,gammas)

! call elapsed(t2)
! print *,"plot 115: ",t2-t1




!  Create a table listing the accuracy with which PS is evaluated



! call elapsed(t3)
! nsigma = 4
! ngamma = 8

! allocate(sigmas(nsigma+1),gammas(ngamma+1),range_labels(8))

! sigmas(1)  = 0.00d0
! sigmas(2)  = 0.25d0
! sigmas(3)  = 0.50d0
! sigmas(4)  = 0.75d0
! sigmas(5)  = 1.00d0

! gammas(1)  = 100
! gammas(2)  = 500
! gammas(3)  = 1000
! gammas(4)  = 5000
! gammas(5)  = 10000
! gammas(6)  = 50000
! gammas(7)  = 100000
! gammas(8)  = 500000
! gammas(9)  = 1000000

! label1     = "100 to 500"
! label2     = "500 to 1,000"
! label3     = "1,000 to 5,000"
! label4     = "5,000 to 10,000"
! label5     = "10,000 to 50,000"
! label6     = "50,000 to 100,000"
! label7     = "100,000 to 500,000"
! label8     = "500,000 to 1,000,000"

! range_labels(1) = label1
! range_labels(2) = label2
! range_labels(3) = label3
! range_labels(4) = label4
! range_labels(5) = label5
! range_labels(6) = label6
! range_labels(7) = label7
! range_labels(8) = label8


! allocate(derrs1(ncount,ncount))
! allocate(derrs2(ncount,ncount))

! iw = 1001
! open(iw,FILE="table102.tex")
! write(iw,"(A)") "\begin{tabular}{ccc@{\hspace{1em}}ccc}"
! write(iw,"(A)") "\toprule"
! write(iw,"(A)",advance="no") "Range of $\gamma$ & Range of $\sigma$ & Max Error  &"
! write(iw,"(A)",advance="no") "Range of $\gamma$ & Range of $\sigma$ & Max Error \\"
! write(iw,"(A)") "\midrule"

! do i=1,ceiling(ngamma/2.0d0)
! do j=1,nsigma

! sigma1 = sigmas(j)
! sigma2 = sigmas(j+1)

! ! Find the first error
! gamma1 = gammas(i)
! gamma2 = gammas(i+1)


! derr1  = 0
! derr2  = 0

! derrs1 = 0
! derrs2 = 0

! !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,jj,dd1,dd2,sigma,gamma,n,dmax)
! !$OMP DO COLLAPSE(2)
! do ii=1,ncount
! do jj=1,ncount

! dmax  = 0
! dd1   = (ii-1.0d0)/(ncount-1.0d0)
! dd2   = (jj-1.0d0)/(ncount-1.0d0)

! sigma = sigma1 + (sigma2-sigma1)*dd1
! gamma = gamma1 + (gamma2-gamma1)*dd2
! n     = sigma*gamma

! call ps_accuracy(chebdata,gamma,n,dmax)
! derrs1(ii,jj) = dmax

! end do
! end do
! !$OMP END DO
! !$OMP END PARALLEL

! derr1 = maxval(derrs1)
! print *,"   ",gamma1,gamma2,sigma1,sigma2,derr1


! ! Find the second error
! i2 = ceiling(ngamma/2.0d0)+i

! if (i2 .gt. ngamma) then
! derr2 = -1
! else
! gamma1 = gammas(i2)
! gamma2 = gammas(i2+1)


! !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,jj,dd1,dd2,sigma,gamma,n,dmax)
! !$OMP DO COLLAPSE(2)
! do ii=1,ncount
! do jj=1,ncount

! dd1   = (ii-1.0d0)/(ncount-1.0d0)
! dd2   = (jj-1.0d0)/(ncount-1.0d0)
! sigma = sigma1 + (sigma2-sigma1)*dd1
! gamma = gamma1 + (gamma2-gamma1)*dd2
! n     = sigma*gamma
! call ps_accuracy(chebdata,gamma,n,dmax)

! derrs2(ii,jj) = dmax
! end do
! end do
! !$OMP END DO
! !$OMP END PARALLEL


! derr2 = maxval(derrs2)


! endif

! print *,"   ",gamma1,gamma2,sigma1,sigma2,derr2

! ! Write the first column to the table
! if (j .eq. 1) then
! write(iw,"(A,A)",advance = 'no') range_labels(i)," &"
! else
! write(iw,"(A)") " & "
! endif

! write(iw,"(A,F4.2,A,F4.2,A)", advance='no')  "$",sigmas(j)," - ",sigmas(j+1), " $ & "
! write(iw,"(A)",advance='no') "$"
! call write_table_double(iw,derr1)
! write(iw,"(A)",advance='no') "$ &"


! if (derr2 .eq. -1) then
! write(iw,"(A)") "\\"
! else
! if (j .eq. 1) then
! write(iw,"(A,A)",advance = 'no') range_labels(i2)," &"
! else
! write(iw,"(A)") " & "
! endif

! write(iw,"(A,F4.2,A,F4.2,A)", advance='no')  "$",sigmas(j)," - ",sigmas(j+1), " $ & "
! write(iw,"(A)",advance='no') "$"
! call write_table_double(iw,derr2)
! !write(iw,"(A,I5.5,A)",advance='no') "$",,"$ &"
! write(iw,"(A)",advance='no') "$ \\"

! endif
! write(iw,"(A)")  "\addlinespace[.125em]"


! end do
! write(iw,"(A)")  "\addlinespace[.25em]"
! end do

! write(iw,"(A)") "\bottomrule"
! write(iw,"(A)") "\end{tabular}"
! close(iw)
! deallocate(gammas, sigmas, derrs1, derrs2, range_labels)
! call elapsed(t4)

! print *,"table102 : = ",t4-t3



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Time required to construct PS_n
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! print *,"-------- TIME TO CONSTRUCT REPRESENTATIONS ------------ "

! ylimit1 = 0.0d0
! ylimit2 = 600.0d0

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the time to construct PS as a function of gamma for multiple sigma ...
! !  without reference to Rokhlin
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! call elapsed(t1)
! nfuns   = 4
! nn      = 100
! iflogx  = 1
! allocate(xs(nn), ys(nfuns,nn), sigmas(nn) )

! gamma1    = 10.0d0**2
! gamma2    = 10.0d0**6

! sigmas(1) = 0.10d0
! sigmas(2) = 0.25d0
! sigmas(3) = 0.50d0
! sigmas(4) = 0.60d0

! do ifun=1,nfuns

! sigma = sigmas(ifun)

! !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,dd2,gamma,n,tphase1,tphase2)
! !!$OMP DO 
! do i=1,nn

! dd    = (i-1.0d0)/(nn-1.0d0)
! if (iflogx .eq. 0) then
! gamma = gamma1  + (gamma2-gamma1)*dd
! else
! dd2   = log(gamma1) + dd*(log(gamma2)-log(gamma1))
! gamma = exp(dd2)
! endif

! n      = sigma*gamma
! call ps_time2(chebdata,exp3,gamma,n,tphase1,tphase2)

! tphase2 = tphase2*1000000
! xs(i)       = gamma
! ys(ifun,i)  = tphase2

! end do
! !!$OMP END DO
! !!$OMP END PARALLEL

! end do

! call plot_functions("plot120.py","plot120.pdf",nfuns,"$\gamma$","Time (microseconds)",iflogx,0,   &
!    gamma1,gamma2,ylimit1,ylimit2,nn,xs,ys,"lower right",&
!    "$\sigma=0.10$*$\sigma=0.25$*$\sigma=0.50$*$\sigma=0.60$*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,sigmas)
! call elapsed(t2)
! print *,"plot 120:",t2-t1



! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the time to construct  PS as a function of gamma for several sigma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call elapsed(t1)
! nfuns   = 4
! nn      = 100
! iflogx  = 1
! allocate(xs(nn), ys(nfuns,nn), sigmas(nn) )

! gamma1 = 10.0d0**2
! gamma2 = 10.0d0**6

! sigmas(1) = 0.62d0
! sigmas(2) = 0.63d0
! sigmas(3) = 2/pi
! sigmas(4) = 0.66d0

! do ifun=1,nfuns

! sigma = sigmas(ifun)

! !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,dd2,gamma,n,tphase1,tphase2)
! !!$OMP DO 
! do i=1,nn

! dd    = (i-1.0d0)/(nn-1.0d0)
! if (iflogx .eq. 0) then
! gamma = gamma1  + (gamma2-gamma1)*dd
! else
! dd2   = log(gamma1) + dd*(log(gamma2)-log(gamma1))
! gamma = exp(dd2)
! endif

! n      = sigma*gamma
! call ps_time2(chebdata,exp3,gamma,n,tphase1,tphase2)

! tphase2 = tphase2*1000000


! xs(i)       = gamma
! ys(ifun,i)  = tphase2

! end do
! !!$OMP END DO
! !!$OMP END PARALLEL

! end do

! call plot_functions("plot121.py","plot121.pdf",nfuns,"$\gamma$","Time (microseconds)",iflogx,0,   &
!    gamma1,gamma2,ylimit1,ylimit2,nn,xs,ys,"upper left",&
!    "$\sigma=0.62$*$\sigma=0.63$*$\sigma=\\frac{2}{\pi}$*$\sigma=0.66$*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,sigmas)
! call elapsed(t2)
! print *,"plot 121:",t2-t1

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the time to construct PS as a function of gamma for several sigma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call elapsed(t1)
! nfuns   = 4
! nn      = 100
! iflogx  = 1
! allocate(xs(nn), ys(nfuns,nn), sigmas(nn) )

! gamma1 = 10.0d0**2
! gamma2 = 10.0d0**6

! sigmas(1) = 0.7d0
! sigmas(2) = 0.8d0
! sigmas(3) = 0.9d0
! sigmas(4) = 1.0d0

! do ifun=1,nfuns

! sigma = sigmas(ifun)

! !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,dd2,gamma,n,tphase1,tphase2)
! !!$OMP DO 
! do i=1,nn

! dd    = (i-1.0d0)/(nn-1.0d0)
! if (iflogx .eq. 0) then
! gamma = gamma1  + (gamma2-gamma1)*dd
! else
! dd2   = log(gamma1) + dd*(log(gamma2)-log(gamma1))
! gamma = exp(dd2)
! endif

! n      = sigma*gamma
! call ps_time2(chebdata,exp3,gamma,n,tphase1,tphase2)

! tphase2 = tphase2*1000000

! xs(i)       = gamma
! ys(ifun,i)  = tphase2

! end do
! !!$OMP END DO
! !!$OMP END PARALLEL

! end do

! call plot_functions("plot122.py","plot122.pdf",nfuns,"$\gamma$","Time (microseconds)",iflogx,0,   &
!    gamma1,gamma2,ylimit1,ylimit2,nn,xs,ys,"upper left",&
!    "$\sigma=0.70$*$\sigma=0.80$*$\sigma=0.90$*$\sigma=1.00$*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,sigmas)
! call elapsed(t2)
! print *,"plot 122: ",t2-t1



! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the time to construct  PS as a function of sigma several gamma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call elapsed(t1)

! nfuns = 4
! allocate(xs(nn), ys(nfuns,nn), gammas(nfuns) )

! gammas(1) = 100
! gammas(2) = 200
! gammas(3) = 500
! gammas(4) = 1000

! sigma1 = 0.0d0
! sigma2 = 1.0d0

! do ifun=1,nfuns
! gamma = gammas(ifun)

! !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,sigma,xi,n,tphase1,tphase2)
! !!$OMP DO 
! do i=1,nn

! dd     = (i-1.0d0)/(nn-1.0d0)
! sigma  = sigma1 + (sigma2-sigma1)*dd

! xi     = sigma*gamma
! n      = xi
! n      = sigma*gamma
! call ps_time2(chebdata,exp3,gamma,n,tphase1,tphase2)
! tphase2 = tphase2*1000000

! xs(i)      = sigma
! ys(ifun,i) = tphase2


! end do
! !!$OMP END DO
! !!$OMP END PARALLEL
! end do

! call plot_functions("plot123.py","plot123.pdf",nfuns,"$\sigma$","Time (microseconds)",0,0,   &
!    0.0d0,1.0d0,ylimit1,ylimit2,nn,xs,ys,"upper right",&
!    "$\gamma=100$*$\gamma=200$*$\gamma=500$*$\gamma=1000$*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,gammas)

! call elapsed(t2)
! print *,"plot 123: ",t2-t1


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the time to construct  PS as a function of sigma several gamma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call elapsed(t1)

! nfuns = 4
! allocate(xs(nn), ys(nfuns,nn), gammas(nfuns) )

! gammas(1) = 10**3
! gammas(2) = 10**4
! gammas(3) = 10**5
! gammas(4) = 10**6

! sigma1 = 0.0d0
! sigma2 = 1.0d0

! do ifun=1,nfuns
! gamma = gammas(ifun)

! !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,sigma,xi,n,tphase1,tphase2)
! !!$OMP DO 
! do i=1,nn

! dd     = (i-1.0d0)/(nn-1.0d0)
! sigma  = sigma1 + (sigma2-sigma1)*dd

! xi     = sigma*gamma
! n      = xi
! n      = sigma*gamma
! call ps_time2(chebdata,exp3,gamma,n,tphase1,tphase2)
! tphase2 = tphase2*1000000

! xs(i)      = sigma
! ys(ifun,i) = tphase2


! end do
! !!$OMP END DO
! !!$OMP END PARALLEL
! end do

! call plot_functions("plot124.py","plot124.pdf",nfuns,"$\sigma$","Time (microseconds)",0,0,   &
!    0.0d0,1.0d0,ylimit1,ylimit2,nn,xs,ys,"upper right",&
!    "$\gamma=10^3$*$\gamma=10^4$*$\gamma=10^5$*$\gamma=10^6$*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,gammas)

! call elapsed(t2)
! print *,"plot 124: ",t2-t1

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the time to construct  PS as a function of sigma several gamma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call elapsed(t1)

! nfuns = 4
! allocate(xs(nn), ys(nfuns,nn), gammas(nfuns) )

! gammas(1) = 0.7d0*10**6
! gammas(2) = 0.8d0*10**6
! gammas(3) = 0.9d0*10**6
! gammas(4) = 1.0d0*10**6

! sigma1 = 0.0d0
! sigma2 = 1.0d0

! do ifun=1,nfuns
! gamma = gammas(ifun)

! !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,sigma,xi,n,tphase1,tphase2)
! !!$OMP DO 
! do i=1,nn

! dd     = (i-1.0d0)/(nn-1.0d0)
! sigma  = sigma1 + (sigma2-sigma1)*dd

! xi     = sigma*gamma
! n      = xi
! n      = sigma*gamma
! call ps_time2(chebdata,exp3,gamma,n,tphase1,tphase2)
! tphase2 = tphase2*1000000

! xs(i)      = sigma
! ys(ifun,i) = tphase2


! end do
! !!$OMP END DO
! !!$OMP END PARALLEL
! end do

! call plot_functions("plot125.py","plot125.pdf",nfuns,"$\sigma$","Time (microseconds)",0,0,   &
!    0.0d0,1.0d0,ylimit1,ylimit2,nn,xs,ys,"upper right",&
!    "$\gamma=0.7 \\times 10^6$*$\gamma=0.8 \\times 10^6$*$\gamma=0.09 \\times 10^6$*$\gamma=1.0 \\times 10^6$*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,gammas)

! call elapsed(t2)
! print *,"plot 125: ",t2-t1


! ylimit1 = 0.0d0
! ylimit2 = 2000.0d0

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Comparison with Ospiov-Rokhlin function of gamma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! nfuns = 2


! call elapsed(t1)
! nfuns   = 2
! nn      = 100
! iflogx  = 0
! allocate(xs(nn), ys(nfuns,nn), sigmas(nn) )

! gamma1 = 100.0d0
! gamma2 = 5000.0d0
! sigma  = 0.10d0

! !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,dd2,gamma,n,tphase1,tphase2,trokhlin)
! !!$OMP DO 
! do i=1,nn

! dd    = (i-1.0d0)/(nn-1.0d0)
! if (iflogx .eq. 0) then
! gamma = gamma1  + (gamma2-gamma1)*dd
! else
! dd2   = log(gamma1) + dd*(log(gamma2)-log(gamma1))
! gamma = exp(dd2)
! endif

! n      = sigma*gamma
! call ps_time1(chebdata,exp3,gamma,n,tphase1,tphase2,trokhlin)

! tphase1  =  tphase1*1000000
! tphase2  =  tphase2*1000000
! trokhlin = trokhlin*1000000

! xs(i)       = gamma
! ys(1,i)     = tphase2
! ys(2,i)     = trokhlin

! end do
! !!$OMP END DO
! !!$OMP END PARALLEL


! call plot_functions("plot130.py","plot130.pdf",nfuns,"$\gamma$","Time (microseconds)",iflogx,0,   &
!    gamma1,gamma2,ylimit1,ylimit2,nn,xs,ys,"upper left",&
!    "Phase function*Osipov-Xiao-Rokhlin*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,sigmas)
! call elapsed(t2)
! print *,"plot 130:",t2-t1


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Comparison with Ospiov-Rokhlin
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! nfuns = 2


! call elapsed(t1)
! nfuns   = 2
! nn      = 100
! iflogx  = 0
! allocate(xs(nn), ys(nfuns,nn), sigmas(nn) )

! gamma1 = 100.0d0
! gamma2 = 5000.0d0
! sigma  = 0.25d0

! !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,dd2,gamma,n,tphase1,tphase2,trokhlin)
! !!$OMP DO 
! do i=1,nn

! dd    = (i-1.0d0)/(nn-1.0d0)
! if (iflogx .eq. 0) then
! gamma = gamma1  + (gamma2-gamma1)*dd
! else
! dd2   = log(gamma1) + dd*(log(gamma2)-log(gamma1))
! gamma = exp(dd2)
! endif

! n      = sigma*gamma
! call ps_time1(chebdata,exp3,gamma,n,tphase1,tphase2,trokhlin)

! tphase1  = tphase1*1000000
! tphase2  = tphase2*1000000
! trokhlin = trokhlin*1000000

! xs(i)       = gamma
! ys(1,i)     = tphase2
! ys(2,i)     = trokhlin

! end do
! !!$OMP END DO
! !!$OMP END PARALLEL

! call plot_functions("plot131.py","plot131.pdf",nfuns,"$\gamma$","Time (mircoseconds)",iflogx,0,   &
!    gamma1,gamma2,ylimit1,ylimit2,nn,xs,ys,"upper left",&
!    "Phase function*Osipov-Xiao-Rokhlin*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,sigmas)
! call elapsed(t2)
! print *,"plot 131:",t2-t1


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Comparison with Ospiov-Rokhlin
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! nfuns = 2


! call elapsed(t1)
! nfuns   = 2
! nn      = 100
! iflogx  = 0
! allocate(xs(nn), ys(nfuns,nn), sigmas(nn) )

! gamma1 = 100.0d0
! gamma2 = 5000.0d0
! sigma  = 0.50d0

! !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,dd2,gamma,n,tphase1,tphase2,trokhlin)
! !!$OMP DO 
! do i=1,nn

! dd    = (i-1.0d0)/(nn-1.0d0)
! if (iflogx .eq. 0) then
! gamma = gamma1  + (gamma2-gamma1)*dd
! else
! dd2   = log(gamma1) + dd*(log(gamma2)-log(gamma1))
! gamma = exp(dd2)
! endif

! n      = sigma*gamma
! call ps_time1(chebdata,exp3,gamma,n,tphase1,tphase2,trokhlin)

! tphase1  = tphase1*1000000
! tphase2  = tphase2*1000000
! trokhlin = trokhlin*1000000

! xs(i)       = gamma
! ys(1,i)     = tphase2
! ys(2,i)     = trokhlin


! end do
! !!$OMP END DO
! !!$OMP END PARALLEL


! call plot_functions("plot132.py","plot132.pdf",nfuns,"$\gamma$","Time (microseconds)",iflogx,0,   &
!    gamma1,gamma2,ylimit1,ylimit2,nn,xs,ys,"upper left",&
!    "Phase function*Osipov-Xiao-Rokhlin*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,sigmas)
! call elapsed(t2)
! print *,"plot 132:",t2-t1




! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Comparison with Ospiov-Rokhlin function of sigma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ylimit1 = 0.0d0
! ylimit2 = 600.0d0


! call elapsed(t1)
! nfuns   = 2
! nn      = 100
! iflogx  = 0
! allocate(xs(nn), ys(nfuns,nn) )

! gamma  = 100
! sigma1 = 0
! sigma2 = 1.0d0

! !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,dd2,sigma,n,tphase1,tphase2,trokhlin)
! !!$OMP DO 
! do i=1,nn

! dd    = (i-1.0d0)/(nn-1.0d0)
! sigma = sigma1  + (sigma2-sigma1)*dd

! n      = sigma*gamma
! call ps_time1(chebdata,exp3,gamma,n,tphase1,tphase2,trokhlin)

! tphase1  =  tphase1*1000000
! tphase2  =  tphase2*1000000
! trokhlin = trokhlin*1000000

! xs(i)       = sigma
! ys(1,i)     = tphase2
! ys(2,i)     = trokhlin

! end do
! !!$OMP END DO
! !!$OMP END PARALLEL

! call plot_functions("plot133.py","plot133.pdf",nfuns,"$\sigma$","Time (microseconds)",iflogx,0,   &
!    sigma1,sigma2,ylimit1,ylimit2,nn,xs,ys,"upper right",&
!    "Phase function*Osipov-Xiao-Rokhlin*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys)
! call elapsed(t2)
! print *,"plot 133:",t2-t1



! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Comparison with Ospiov-Rokhlin function of sigma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! nfuns = 2


! call elapsed(t1)
! nfuns   = 2
! nn      = 100
! iflogx  = 0
! allocate(xs(nn), ys(nfuns,nn) )

! gamma  = 500
! sigma1 = 0
! sigma2 = 1.0d0

! !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,dd2,sigma,n,tphase1,tphase2,trokhlin)
! !!$OMP DO 
! do i=1,nn

! dd    = (i-1.0d0)/(nn-1.0d0)
! sigma = sigma1  + (sigma2-sigma1)*dd

! n      = sigma*gamma
! call ps_time1(chebdata,exp3,gamma,n,tphase1,tphase2,trokhlin)

! tphase1  = tphase1*1000000
! tphase2  = tphase2*1000000
! trokhlin = trokhlin*1000000

! xs(i)       = sigma
! ys(1,i)     = tphase2
! ys(2,i)     = trokhlin

! end do
! !!$OMP END DO
! !!$OMP END PARALLEL


! call plot_functions("plot134.py","plot134.pdf",nfuns,"$\sigma$","Time (microseconds)",iflogx,0,   &
!    sigma1,sigma2,ylimit1,ylimit2,nn,xs,ys,"upper right",&
!    "Phase function*Osipov-Xiao-Rokhlin*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys)
! call elapsed(t2)
! print *,"plot 134:",t2-t1


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Comparison with Ospiov-Rokhlin function of sigma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! nfuns = 2


! call elapsed(t1)
! nfuns   = 2
! nn      = 100
! iflogx  = 0
! allocate(xs(nn), ys(nfuns,nn) )

! gamma  = 1000
! sigma1 = 0
! sigma2 = 1.0d0

! !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,dd2,sigma,n,tphase1,tphase2,trokhlin)
! !!$OMP DO 
! do i=1,nn

! dd    = (i-1.0d0)/(nn-1.0d0)
! sigma = sigma1  + (sigma2-sigma1)*dd

! n      = sigma*gamma
! call ps_time1(chebdata,exp3,gamma,n,tphase1,tphase2,trokhlin)

! tphase1  = tphase1*1000000
! tphase2  = tphase2*1000000
! trokhlin = trokhlin*1000000

! xs(i)       = sigma
! ys(1,i)     = tphase2
! ys(2,i)     = trokhlin


! end do
! !!$OMP END DO
! !!$OMP END PARALLEL

! call plot_functions("plot135.py","plot135.pdf",nfuns,"$\sigma$","Time (microseconds)",iflogx,0,   &
!    sigma1,sigma2,ylimit1,ylimit2,nn,xs,ys,"upper left",&
!    "Phase function*Osipov-Xiao-Rokhlin*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys)
! call elapsed(t2)
! print *,"plot 135:",t2-t1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Table giving the time required to construct PS as a function of gamma and compared with
!  Osipov-Xiao-Rokhlin
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call elapsed(t3)
! nsigma = 4
! ngamma = 8

! allocate(sigmas(nsigma+1),gammas(ngamma+1),range_labels(8))

! sigmas(1)  = 0.00d0
! sigmas(2)  = 0.25d0
! sigmas(3)  = 0.50d0
! sigmas(4)  = 0.75d0
! sigmas(5)  = 1.00d0

! gammas(1)  = 100
! gammas(2)  = 500
! gammas(3)  = 1000
! gammas(4)  = 5000
! gammas(5)  = 10000
! gammas(6)  = 50000
! gammas(7)  = 100000
! gammas(8)  = 500000
! gammas(9)  = 1000000

! label1     = "100 to 500"
! label2     = "500 to 1,000"
! label3     = "1,000 to 5,000"
! label4     = "5,000 to 10,000"
! label5     = "10,000 to 50,000"
! label6     = "50,000 to 100,000"
! label7     = "100,000 to 500,000"
! label8     = "500,000 to 1,000,000"

! range_labels(1) = label1
! range_labels(2) = label2
! range_labels(3) = label3
! range_labels(4) = label4
! range_labels(5) = label5
! range_labels(6) = label6
! range_labels(7) = label7
! range_labels(8) = label8

! allocate(times1(ncount,ncount))
! allocate(times2(ncount,ncount))
! allocate(times3(ncount,ncount))

! iw = 1001
! open(iw,FILE="table103.tex")
! write(iw,"(A)") "\begin{tabular}{ccccccr}"
! write(iw,"(A)") "\toprule"
! write(iw,"(A)",advance="no") "Range of $\gamma$ & Range of $\sigma$ & Average Time   & Average time"
! write(iw,"(A)")              " & Average time        & Ratio\\"
! write(iw,"(A)") "              &                   & Unaccelerated  & Accelerated    & Rokhlin, et. al.\\"
! write(iw,"(A)") "\midrule"

! do i=1,ngamma
! do j=1,nsigma

! sigma1 = sigmas(j)
! sigma2 = sigmas(j+1)

! gamma1 = gammas(i)
! gamma2 = gammas(i+1)

! call elapsed(t1)
! !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,jj,dd1,dd2,sigma,gamma,n,tphase1,tphase2,tphase3,trokhlin) 
! !$OMP DO COLLAPSE(2)
! do ii=1,ncount
! do jj=1,ncount

! dd1   = (ii-1.0d0)/(ncount-1.0d0)
! dd2   = (jj-1.0d0)/(ncount-1.0d0)
! sigma = sigma1 + (sigma2-sigma1)*dd1
! gamma = gamma1 + (gamma2-gamma1)*dd2
! n     = sigma*gamma

! call ps_time1(chebdata,exp3,gamma,n,tphase1,tphase2,trokhlin)

! times1(ii,jj) = tphase1
! times2(ii,jj) = tphase2
! times3(ii,jj) = trokhlin

! end do
! end do
! !$OMP END DO
! !$OMP END PARALLEL
! call elapsed(t2)

! print *,"   ",gamma1,gamma2,sigma1,sigma2,t2-t1

! time1 = sum(times1) / (ncount*ncount)
! time2 = sum(times2) / (ncount*ncount)
! time3 = sum(times3) / (ncount*ncount)

! ! Write to the table
! if (j .eq. 1) then
! write(iw,"(A,A)",advance = 'no') range_labels(i)," &"
! else
! write(iw,"(A)") " & "
! endif

! write(iw,"(A,F4.2,A,F4.2,A)", advance='no')  "$",sigmas(j)," - ",sigmas(j+1), " $ & "

! write(iw,"(A)",advance='no') "$"
! call write_table_double(iw,time1)
! write(iw,"(A)",advance='no') "$ &"

! write(iw,"(A)",advance='no') "$"
! call write_table_double(iw,time2)
! write(iw,"(A)",advance='no') "$ &"

! write(iw,"(A)",advance='no') "$"
! call write_table_double(iw,time3)
! write(iw,"(A)",advance='no') "$ &"

! ratio = time3 / time2
! write(iw,"(A)",advance='no') "$"
! write(iw,"(F8.2)") ratio
! !call write_table_double(iw,ratio)
! write(iw,"(A)",advance='no') "$ \\"

! write(iw,"(A)")  "\addlinespace[.125em]"


! end do
! write(iw,"(A)")  "\addlinespace[.25em]"
! end do


! write(iw,"(A)") "\bottomrule"
! write(iw,"(A)") "\end{tabular}"

! close(iw)

! deallocate(gammas, sigmas, times1, times2, times3, range_labels)

! call elapsed(t4)

! print *,"table103:  = ",t4-t3




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Time required to evaluate PS_n
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



print *,"-------- TIME TO EVALUATE REPRESENTATIONS ------------ "

! ylimit1 =   0.0d0
! ylimit2 = 160.0d0

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the time to construct PS as a function of gamma for multiple sigma ...
! !  without reference to Rokhlin
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! call elapsed(t1)
! nfuns   = 4
! nn      = 100
! iflogx  = 1
! allocate(xs(nn), ys(nfuns,nn), sigmas(nn) )

! gamma1    = 10.0d0**2
! gamma2    = 10.0d0**6

! sigmas(1) = 0.10d0
! sigmas(2) = 0.25d0
! sigmas(3) = 0.50d0
! sigmas(4) = 0.60d0

! do ifun=1,nfuns

! sigma = sigmas(ifun)

! !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,dd2,gamma,n,tphase)
! !!$OMP DO 
! do i=1,nn

! dd    = (i-1.0d0)/(nn-1.0d0)
! if (iflogx .eq. 0) then
! gamma = gamma1  + (gamma2-gamma1)*dd
! else
! dd2   = log(gamma1) + dd*(log(gamma2)-log(gamma1))
! gamma = exp(dd2)
! endif

! n      = sigma*gamma
! call ps_time3(chebdata,exp3,gamma,n,tphase)

! tphase = tphase*1000000000
! xs(i)       = gamma
! ys(ifun,i)  = tphase

! end do
! !!$OMP END DO
! !!$OMP END PARALLEL

! end do

! call plot_functions("plot140.py","plot140.pdf",nfuns,"$\gamma$","Time (nanoseconds)",iflogx,0,   &
!    gamma1,gamma2,ylimit1,ylimit2,nn,xs,ys,"lower right",&
!    "$\sigma=0.10$*$\sigma=0.25$*$\sigma=0.50$*$\sigma=0.60$*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,sigmas)
! call elapsed(t2)
! print *,"plot 140:",t2-t1



! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the time to construct  PS as a function of gamma for several sigma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call elapsed(t1)
! nfuns   = 4
! nn      = 100
! iflogx  = 1
! allocate(xs(nn), ys(nfuns,nn), sigmas(nn) )

! gamma1 = 10.0d0**2
! gamma2 = 10.0d0**6

! sigmas(1) = 0.62d0
! sigmas(2) = 0.63d0
! sigmas(3) = 2/pi
! sigmas(4) = 0.66d0

! do ifun=1,nfuns

! sigma = sigmas(ifun)

! !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,dd2,gamma,n,tphase)
! !!$OMP DO 
! do i=1,nn

! dd    = (i-1.0d0)/(nn-1.0d0)
! if (iflogx .eq. 0) then
! gamma = gamma1  + (gamma2-gamma1)*dd
! else
! dd2   = log(gamma1) + dd*(log(gamma2)-log(gamma1))
! gamma = exp(dd2)
! endif

! n      = sigma*gamma
! call ps_time3(chebdata,exp3,gamma,n,tphase)

! tphase = tphase*1000000000


! xs(i)       = gamma
! ys(ifun,i)  = tphase

! end do
! !!$OMP END DO
! !!$OMP END PARALLEL

! end do

! call plot_functions("plot141.py","plot141.pdf",nfuns,"$\gamma$","Time (nanoseconds)",iflogx,0,   &
!    gamma1,gamma2,ylimit1,ylimit2,nn,xs,ys,"upper left",&
!    "$\sigma=0.62$*$\sigma=0.63$*$\sigma=\\frac{2}{\pi}$*$\sigma=0.66$*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,sigmas)
! call elapsed(t2)
! print *,"plot 141:",t2-t1


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the time to construct PS as a function of gamma for several sigma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call elapsed(t1)
! nfuns   = 4
! nn      = 100
! iflogx  = 1
! allocate(xs(nn), ys(nfuns,nn), sigmas(nn) )

! gamma1 = 10.0d0**2
! gamma2 = 10.0d0**6

! sigmas(1) = 0.7d0
! sigmas(2) = 0.8d0
! sigmas(3) = 0.9d0
! sigmas(4) = 1.0d0

! do ifun=1,nfuns

! sigma = sigmas(ifun)

! !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,dd2,gamma,n,tphase)
! !!$OMP DO 
! do i=1,nn

! dd    = (i-1.0d0)/(nn-1.0d0)
! if (iflogx .eq. 0) then
! gamma = gamma1  + (gamma2-gamma1)*dd
! else
! dd2   = log(gamma1) + dd*(log(gamma2)-log(gamma1))
! gamma = exp(dd2)
! endif

! n      = sigma*gamma
! call ps_time3(chebdata,exp3,gamma,n,tphase)

! tphase = tphase*1000000000

! xs(i)       = gamma
! ys(ifun,i)  = tphase

! end do
! !!$OMP END DO
! !!$OMP END PARALLEL

! end do

! call plot_functions("plot142.py","plot142.pdf",nfuns,"$\gamma$","Time (nanoseconds)",iflogx,0,   &
!    gamma1,gamma2,ylimit1,ylimit2,nn,xs,ys,"upper left",&
!    "$\sigma=0.70$*$\sigma=0.80$*$\sigma=0.90$*$\sigma=1.00$*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,sigmas)
! call elapsed(t2)
! print *,"plot 142: ",t2-t1



! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the time to construct  PS as a function of sigma several gamma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! ylimit1 =   0.0d0
! ! ylimit2 =  150.0d0

! call elapsed(t1)

! nfuns = 4
! allocate(xs(nn), ys(nfuns,nn), gammas(nfuns) )

! gammas(1) = 100
! gammas(2) = 200
! gammas(3) = 500
! gammas(4) = 1000

! sigma1 = 0.0d0
! sigma2 = 1.0d0

! do ifun=1,nfuns
! gamma = gammas(ifun)

! !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,sigma,xi,n,tphase)
! !!$OMP DO 
! do i=1,nn

! dd     = (i-1.0d0)/(nn-1.0d0)
! sigma  = sigma1 + (sigma2-sigma1)*dd

! xi     = sigma*gamma
! n      = xi
! n      = sigma*gamma
! call ps_time3(chebdata,exp3,gamma,n,tphase)
! tphase = tphase*1000000000

! xs(i)      = sigma
! ys(ifun,i) = tphase


! end do
! !!$OMP END DO
! !!$OMP END PARALLEL
! end do

! call plot_functions("plot143.py","plot143.pdf",nfuns,"$\sigma$","Time (nanoseconds)",0,0,   &
!    0.0d0,1.0d0,ylimit1,ylimit2,nn,xs,ys,"upper right",&
!    "$\gamma=100$*$\gamma=200$*$\gamma=500$*$\gamma=1000$*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,gammas)

! call elapsed(t2)
! print *,"plot 143: ",t2-t1


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the time to construct  PS as a function of sigma several gamma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call elapsed(t1)

! nfuns = 4
! allocate(xs(nn), ys(nfuns,nn), gammas(nfuns) )

! gammas(1) = 10**3
! gammas(2) = 10**4
! gammas(3) = 10**5
! gammas(4) = 10**6

! sigma1 = 0.0d0
! sigma2 = 1.0d0

! do ifun=1,nfuns
! gamma = gammas(ifun)

! !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,sigma,xi,n,tphase)
! !!$OMP DO 
! do i=1,nn

! dd     = (i-1.0d0)/(nn-1.0d0)
! sigma  = sigma1 + (sigma2-sigma1)*dd

! xi     = sigma*gamma
! n      = xi
! n      = sigma*gamma
! call ps_time3(chebdata,exp3,gamma,n,tphase)
! tphase = tphase*1000000000

! xs(i)      = sigma
! ys(ifun,i) = tphase


! end do
! !!$OMP END DO
! !!$OMP END PARALLEL
! end do

! call plot_functions("plot144.py","plot144.pdf",nfuns,"$\sigma$","Time (nanoseconds)",0,0,   &
!    0.0d0,1.0d0,ylimit1,ylimit2,nn,xs,ys,"upper right",&
!    "$\gamma=10^3$*$\gamma=10^4$*$\gamma=10^5$*$\gamma=10^6$*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,gammas)

! call elapsed(t2)
! print *,"plot 144: ",t2-t1


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !  Plot the time to construct  PS as a function of sigma several gamma
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call elapsed(t1)

! nfuns = 4
! allocate(xs(nn), ys(nfuns,nn), gammas(nfuns) )

! gammas(1) = 0.7d0*10**6
! gammas(2) = 0.8d0*10**6
! gammas(3) = 0.9d0*10**6
! gammas(4) = 1.0d0*10**6

! sigma1 = 0.0d0
! sigma2 = 1.0d0

! do ifun=1,nfuns
! gamma = gammas(ifun)

! !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,sigma,xi,n,tphase)
! !!$OMP DO 
! do i=1,nn

! dd     = (i-1.0d0)/(nn-1.0d0)
! sigma  = sigma1 + (sigma2-sigma1)*dd

! xi     = sigma*gamma
! n      = xi
! n      = sigma*gamma
! call ps_time3(chebdata,exp3,gamma,n,tphase)
! tphase = tphase*1000000000

! xs(i)      = sigma
! ys(ifun,i) = tphase


! end do
! !!$OMP END DO
! !!$OMP END PARALLEL
! end do

! call plot_functions("plot145.py","plot145.pdf",nfuns,"$\sigma$","Time (nanoseconds)",0,0,   &
!    0.0d0,1.0d0,ylimit1,ylimit2,nn,xs,ys,"upper right",&
!    "$\gamma=0.7 \\times 10^6$*$\gamma=0.8 \\times 10^6$*$\gamma=0.09 \\times 10^6$*$\gamma=1.0 \\times 10^6$*", &
!    "b-*r--*g-.*k:*mx*")

! deallocate(xs,ys,gammas)

! call elapsed(t2)
! print *,"plot 145: ",t2-t1





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Table giving the time required to evaluate PS as a function of gamma and compared with
!  Osipov-Xiao-Rokhlin
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call elapsed(t3)
nsigma = 4
ngamma = 8

allocate(sigmas(nsigma+1),gammas(ngamma+1),range_labels(8))

sigmas(1)  = 0.00d0
sigmas(2)  = 0.25d0
sigmas(3)  = 0.50d0
sigmas(4)  = 0.75d0
sigmas(5)  = 1.00d0

gammas(1)  = 100
gammas(2)  = 500
gammas(3)  = 1000
gammas(4)  = 5000
gammas(5)  = 10000
gammas(6)  = 50000
gammas(7)  = 100000
gammas(8)  = 500000
gammas(9)  = 1000000

label1     = "100 to 500"
label2     = "500 to 1,000"
label3     = "1,000 to 5,000"
label4     = "5,000 to 10,000"
label5     = "10,000 to 50,000"
label6     = "50,000 to 100,000"
label7     = "100,000 to 500,000"
label8     = "500,000 to 1,000,000"

range_labels(1) = label1
range_labels(2) = label2
range_labels(3) = label3
range_labels(4) = label4
range_labels(5) = label5
range_labels(6) = label6
range_labels(7) = label7
range_labels(8) = label8

allocate(times1(ncount2,ncount2))
allocate(times2(ncount2,ncount2))
allocate(times3(ncount2,ncount2))

iw = 1001
open(iw,FILE="table104.tex")
write(iw,"(A)") "\begin{tabular}{cccccr}"
write(iw,"(A)") "\toprule"
write(iw,"(A)",advance="no") "Range of $\gamma$ & Range of $\sigma$ &   Average time"
write(iw,"(A)")              " & Average time        & Ratio\\"
write(iw,"(A)") "              &                   & Phase  & Rokhlin, et. al.\\"
write(iw,"(A)") "\midrule"

do i=1,ngamma
do j=1,nsigma

sigma1 = sigmas(j)
sigma2 = sigmas(j+1)

gamma1 = gammas(i)
gamma2 = gammas(i+1)

call elapsed(t1)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,jj,dd1,dd2,sigma,gamma,n,tphase,trokhlin) 
!$OMP DO COLLAPSE(2)
do ii=1,ncount2
do jj=1,ncount2

dd1   = (ii-1.0d0)/(ncount2-1.0d0)
dd2   = (jj-1.0d0)/(ncount2-1.0d0)
sigma = sigma1 + (sigma2-sigma1)*dd1
gamma = gamma1 + (gamma2-gamma1)*dd2
n     = sigma*gamma

call ps_time4(chebdata,exp3,gamma,n,tphase,trokhlin)

times1(ii,jj) = tphase
times2(ii,jj) = trokhlin

end do
end do
!$OMP END DO
!$OMP END PARALLEL
call elapsed(t2)

print *,"   ",gamma1,gamma2,sigma1,sigma2,t2-t1

time1 = sum(times1) / (ncount2*ncount2)
time2 = sum(times2) / (ncount2*ncount2)

! Write to the table
if (j .eq. 1) then
write(iw,"(A,A)",advance = 'no') range_labels(i)," &"
else
write(iw,"(A)") " & "
endif

write(iw,"(A,F4.2,A,F4.2,A)", advance='no')  "$",sigmas(j)," - ",sigmas(j+1), " $ & "

write(iw,"(A)",advance='no') "$"
call write_table_double(iw,time1)
write(iw,"(A)",advance='no') "$ &"

write(iw,"(A)",advance='no') "$"
call write_table_double(iw,time2)
write(iw,"(A)",advance='no') "$ &"


ratio = time2 / time1
write(iw,"(A)",advance='no') "$"
write(iw,"(F10.2)") ratio
!call write_table_double(iw,ratio)
write(iw,"(A)",advance='no') "$ \\"

write(iw,"(A)")  "\addlinespace[.125em]"


end do
write(iw,"(A)")  "\addlinespace[.25em]"
end do


write(iw,"(A)") "\bottomrule"
write(iw,"(A)") "\end{tabular}"

close(iw)

deallocate(gammas, sigmas, times1, times2, times3, range_labels)

call elapsed(t4)

print *,"table104:  = ",t4-t3


end program

