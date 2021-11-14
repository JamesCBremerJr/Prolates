!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for constructing expansions which allow for the
!  calculation of the Sturm-Liouville eigenvalues for the reduced spheroidal wave
!  equation and certain related quantities.  Since the eigenvalue
!       
!    chi   (gamma).
!       nu
!
!  is discontinuous as a function of the characteristic exponent nu, we use
!  an alternate parameter sigma instead.  It is defined via the formula
!
!
!    sigma  (gamma^2) = xi  (gamma^2) / gamma,                                            (1)
!         nu              nu
!
!  where 
!                               
!    xi   (gamma) = -2/pi  Psi   (0, \gamma)  -1.
!        nu                    nu
!
!  When dnu is an integer, xi is equal to dnu.  Moreover, it is entire as a 
!  function of dnu and gamma.  The scaling in (1) is present because we
!  construct expansions as a function of gamma as well as sigma.
!
!  The "related quantities" are the values of the first and third derivatives
!  of Psi_nu(x;\gamma) at 0.  The second derivative is always equal to the
!  additive inverse of the first derivative.
!
!  The following routines should be regarded as publicly callable:
!
!    prolates_exp1 - adaptively construct piecewise polynomial expansions of 
!      the Sturm-Liouville eigenvalue and related quantities as a function of 
!      sigma with gamma fixed.
!
!    prolates_exp1_eval - evaluate the expansions constructed by prolates_exp1
!
!    prolates_exp1_check - check the accuracy of a set of  expansions produced
!       by prolates_exp1
!
!    write_prolate_exp1 - write the expansion data to a textfile on the disk
!
!    read_prolate_exp1 - read the expansion data from a textfile on the disk
!
!    prolates_exp2 - construct an expansions of the Sturm-Liouville eigenvalue
!      and related quantities as functions of gamma and sigma.  The expansions
!      are adaptive in sigma but not gamma.
!
!    prolates_exp2_eval - evaluate the expansions constructed by prolates_exp2
!
!    prolates_exp2_check - check the accuracy of a set of  expansions produced
!       by prolates_exp2
!
!    prolates_exp3 - adaptively construct piecewise polynomials expansions of
!      the eigenvalue and related quantities as functions of gamma and sigma.
!
!    prolates_exp3_eval - evaluate the expansions constructed by prolates_exp3
!
!    prolates_exp3_check - check the accuracy of the expansions constructed
!      by prolates_exp3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module prolates_exps

use utils
use chebyshev
use prolates
use prolates_phase

! The exp1_data structure stores expansions of the relevant quantities as functions
! of sigma for a fixed gamma.
type     exp1_data
double precision              :: gamma
double precision              :: sigma1
double precision              :: sigma2
integer                       :: ncost 

integer                       :: k
integer                       :: nints
double precision, allocatable :: ab(:,:)


double precision, allocatable :: chi_coefs(:,:)
double precision, allocatable :: ap_coefs(:,:)
double precision, allocatable :: appp_coefs(:,:)
end type exp1_data


!  The structure exp2_data stores an expansion of chi as a function of sigma and gamma
!  for a range of values of gamma and sigma.  The expansion is adaptive in sigma
!  BUT NOT gamma!  
type     exp2_data
double precision              :: gamma1
double precision              :: gamma2
double precision              :: sigma1
double precision              :: sigma2
integer                       :: ncost

integer                       :: k
double precision, allocatable :: xscheb(:)
type(exp1_data), allocatable  :: exps(:)
end type exp2_data


!  The structure exp2_data stores an expansion of chi as a function of sigma and gamma
!  for a range of values of gamma and sigma.  The expansion is adaptive in both
!  variables
type     exp3_data
double precision              :: gamma1
double precision              :: gamma2
double precision              :: sigma1
double precision              :: sigma2
integer                       :: ncost

integer                       :: nints
double precision, allocatable :: ab(:,:)
type(exp2_data), allocatable  :: exps(:)
end type exp3_data


contains

subroutine prolates_exp1_write(iw,exp1)
implicit double precision (a-h,o-z)
type(exp1_data)               :: exp1
!
!
!

write (iw,"(D24.16)")           exp1%gamma
write (iw,"(D24.16)")           exp1%sigma1
write (iw,"(D24.16)")           exp1%sigma2
write (iw,"(I8)")               exp1%k
write (iw,"(I8)")               exp1%nints
write (iw,"(I8)")               exp1%ncost
write (iw,"(D24.16)")           exp1%ab
write (iw,"(D24.16)")           exp1%chi_coefs
write (iw,"(D24.16)")           exp1%ap_coefs
write (iw,"(D24.16)")           exp1%appp_coefs

end subroutine

subroutine prolates_exp1_read(iw,exp1)
implicit double precision (a-h,o-z)
type(exp1_data), intent(out)               :: exp1
!
!
!

read (iw,"(D24.16)")           exp1%gamma
read (iw,"(D24.16)")           exp1%sigma1
read (iw,"(D24.16)")           exp1%sigma2
read (iw,"(I8)")               exp1%k
read (iw,"(I8)")               exp1%nints
read (iw,"(I8)")               exp1%ncost

k     = exp1%k
nints = exp1%nints

allocate(exp1%ab(2,nints))
allocate(exp1%chi_coefs(k,nints))
allocate(exp1%ap_coefs(k,nints))
allocate(exp1%appp_coefs(k,nints))

read (iw,"(D24.16)")           exp1%ab
read (iw,"(D24.16)")           exp1%chi_coefs
read (iw,"(D24.16)")           exp1%ap_coefs
read (iw,"(D24.16)")           exp1%appp_coefs

end subroutine


subroutine prolates_exp2_write(iw,exp2)
implicit double precision (a-h,o-z)
type(exp2_data)               :: exp2
!
!
!

write (iw,"(D24.16)")           exp2%gamma1
write (iw,"(D24.16)")           exp2%gamma2
write (iw,"(D24.16)")           exp2%sigma1
write (iw,"(D24.16)")           exp2%sigma2
write (iw,"(I8)")               exp2%ncost
write (iw,"(I8)")               exp2%k
write (iw,"(D24.16)")           exp2%xscheb

k = exp2%k
do i=1,k
call prolates_exp1_write(iw,exp2%exps(i))
end do

end subroutine

subroutine prolates_exp2_read(iw,exp2)
implicit double precision (a-h,o-z)
type(exp2_data), intent(out)               :: exp2
!
!
!

read (iw,"(D24.16)")           exp2%gamma1
read (iw,"(D24.16)")           exp2%gamma2
read (iw,"(D24.16)")           exp2%sigma1
read (iw,"(D24.16)")           exp2%sigma2
read (iw,"(I8)")               exp2%ncost
read (iw,"(I8)")               exp2%k

k = exp2%k
allocate(exp2%xscheb(k))
allocate(exp2%exps(k))

read (iw,"(D24.16)")           exp2%xscheb

do i=1,k
call prolates_exp1_read(iw,exp2%exps(i))
end do

end subroutine


subroutine prolates_exp3_write(iw,exp3)
implicit double precision (a-h,o-z)
type(exp3_data)               :: exp3
!
!
!

write (iw,"(D24.16)")           exp3%gamma1
write (iw,"(D24.16)")           exp3%gamma2
write (iw,"(D24.16)")           exp3%sigma1
write (iw,"(D24.16)")           exp3%sigma2
write (iw,"(I8)")               exp3%ncost
write (iw,"(I8)")               exp3%nints
write (iw,"(D24.16)")           exp3%ab

nints = exp3%nints
do i=1,nints
call prolates_exp2_write(iw,exp3%exps(i))
end do

end subroutine


subroutine prolates_exp3_read(iw,exp3)
implicit double precision (a-h,o-z)
type(exp3_data), intent(out)               :: exp3
!
!
!

read (iw,"(D24.16)")           exp3%gamma1
read (iw,"(D24.16)")           exp3%gamma2
read (iw,"(D24.16)")           exp3%sigma1
read (iw,"(D24.16)")           exp3%sigma2
read (iw,"(I8)")               exp3%ncost
read (iw,"(I8)")               exp3%nints


nints = exp3%nints
allocate(exp3%ab(2,nints))
allocate(exp3%exps(nints))
read (iw,"(D24.16)")           exp3%ab

do i=1,nints
call prolates_exp2_read(iw,exp3%exps(i))
end do

end subroutine


subroutine prolates_exp1(chebdata,sigma1,sigma2,gamma,expdata)
implicit double precision (a-h,o-z)
type(chebexps_data)                        :: chebdata
type(exp1_data), intent(out)               :: expdata
!
!  Construct piecewise Chebyshev expansions of the Sturm-Liouville eigenvalue chi 
!  and the related quantiites as a function of the parameter sigma.  The parameter 
!  sigma is allowed to vary over a user-specified range.
!
!  Input parameters:
!     chebdata - the structure returned by chebexps which, among other things,
!       specifies the order of the piecewise Chebyshev expansions to use
!    gamma - the (fixed) bandlimit
!    (sigma1,sigma2) - the range of values of sigma for the expansion
!
!  Output parameters:
!    expdata - a data structure describing the expansions

type(prolates_phase_data)     :: phase
double precision, allocatable :: xi_ab(:,:), xi_vals(:,:)
double precision, allocatable :: chi_ab(:,:), chi_vals(:,:), ap_vals(:,:),appp_vals(:,:)
double precision, allocatable :: ap_vals0(:,:),appp_vals0(:,:)
double precision, allocatable :: chi_coefs(:,:), ap_coefs(:,:), appp_coefs(:,:)
double precision, allocatable :: stack(:,:), vals0(:), coefs0(:), coefsps(:)
double precision, allocatable :: ap0(:),appp0(:)

pi       = acos(-1.0d0)
maxstack = 10000
maxints  = 10000
maxiters = 200           ! max # of iterations for bisection
k        = chebdata%k
inormal  = 1

allocate(stack(k,maxints),coefs0(k),vals0(k))
allocate(xi_ab(2,maxints),xi_vals(k,maxints))
allocate(chi_ab(2,maxints),chi_vals(k,maxints))
allocate(ap_vals0(k,maxints),appp_vals0(k,maxints))
allocate(ap_vals(k,maxints),appp_vals(k,maxints))

allocate(ap0(k),appp0(k))

nn    = k*0.25d0

eps0 = epsilon(0.0d0)
eps      = 1.0d-12
epscoefs = 1.0d-12

if (eps0 .lt. 1.0d-16) then
eps      = 1.0d-15
epscoefs = 1.0d-15
endif

if (eps0 .lt. 1.0d-30) then
eps      = 1.0d-20
epscoefs = 1.0d-20
endif

if (gamma .ge. 10**5) eps=eps*10

!
! Find lower and upper bounds for the range of chi corresponding to the range (0,gamma) of xi
!

inormal  = 1
idx1     = floor(sigma1*gamma)
idx2     = ceiling(sigma2*gamma)

call prolates_integer(inormal,gamma,idx1,chi1,n1,n2,coefsps)
call prolates_integer(inormal,gamma,idx2,chi2,n1,n2,coefsps)

chi1 = chi1*0.95d0
chi2 = chi2*1.05d0

!
! Construct an adaptive discretization of xi as a function of chi over that range
!
nints_xi   = 0
nstack     = 1
stack(1,1) = chi1
stack(2,1) = chi2

do while (nstack > 0)
a      = stack(1,nstack)
b      = stack(2,nstack)
nstack = nstack-1

do i=1,k
chi = chebdata%xs(i)*(b-a)/2 + (b+a)/2
call elapsed(t1)
call prolates_phase3(chebdata,gamma,chi,aval,apval,appval,apppval)
call elapsed(t2)
vals0(i) = -2/pi*aval-1
ap0(i)   = apval
appp0(i) = apppval
end do


ifsplit = 0

coefs0 = matmul(chebdata%u,vals0)
coefs0 = abs(coefs0) / abs(coefs0(1))
dd1    = maxval(coefs0(k-nn+1:k))

coefs0 = matmul(chebdata%u,ap0)
coefs0 = abs(coefs0) / abs(coefs0(1))
dd2    = maxval(coefs0(k-nn+1:k))

if (dd1 .gt. eps) ifsplit = 1
if (dd2 .gt. eps) ifsplit = 1

!print *,"            ",a,b,dd1,dd2,ifsplit

if (ifsplit .eq. 1) then

if (nstack+2 .gt. maxstack) then
print *,"prolates_exp1: stack overflow (xi exp)"
print *,"gamma  = ",gamma
print *,"sigma1 = ",sigma1
print *,"sigma2 = ",sigma2
stop
endif

nstack          = nstack+1
stack(1,nstack) = (a+b)/2
stack(2,nstack) = b

nstack          = nstack+1
stack(1,nstack) = a
stack(2,nstack) = (a+b)/2
else

if (nints_xi .eq. maxints) then
print *,"prolates_exp1: too many intervals (xi exp)"
print *,"gamma  = ",gamma
print *,"sigma1 = ",sigma1
print *,"sigma2 = ",sigma2
stop
endif

nints_xi               = nints_xi+1
xi_ab(1,nints_xi)      = a
xi_ab(2,nints_xi)      = b
xi_vals(:,nints_xi)    = vals0
ap_vals0(:,nints_xi)   = ap0
appp_vals0(:,nints_xi) = appp0
endif
end do

!call prini("nints_xi = ",nints_xi)
!call prin2("xi_ab   = ",xi_ab(:,1:nints_xi))
!call prin2("xi_vals = ",xi_vals(:,1:nints_xi))

!
! Construct an expansion of the inverse function via bisection
!
nints_chi  = 0 
nstack     = 1
stack(1,1) = sigma1
stack(2,1) = sigma2

if (sigma1 .lt. 2/pi .AND. sigma2 .gt. 2/pi) then
nstack     = 2
stack(1,1) = 2/pi
stack(2,1) = sigma2
stack(1,2) = sigma1
stack(2,2) = 2/pi
endif

do while (nstack > 0) 
a      = stack(1,nstack)
b      = stack(2,nstack)
nstack = nstack-1

do i=1,k
xi = (b-a)/2 * chebdata%xs(i) + (b+a)/2
xi = xi*gamma


  ! first find the interval containing the value of chi corresponding to 
  ! xi
  do int0=1,nints_xi-1
  if (xi .le. xi_vals(k,int0)) exit
  end do

  a0 = xi_ab(1,int0)
  b0 = xi_ab(2,int0)

  xi1  = xi_vals(1,int0)
  xi2  = xi_vals(k,int0)

  chi1 = xi_ab(1,int0)
  chi2 = xi_ab(2,int0)
  
  ! use the discretization to narrow the range before bisection
  do ii=1,k-1
  if (xi_vals(ii+1,int0) .ge. xi) then
  xi1  = xi_vals(ii,int0)
  xi2  = xi_vals(ii+1,int0)
  chi1 = (b0-a0)/2*chebdata%xs(ii)   + (b0+a0)/2
  chi2 = (b0-a0)/2*chebdata%xs(ii+1) + (b0+a0)/2
  exit
  endif
  end do

  ! use bisection to compute chi as a function of xi
  do iter=1,maxiters
  chi0 = (chi1+chi2)/2
  call chebeval(a0,b0,k,chebdata%xs,xi_vals(:,int0),chi0,xi0)
  if (xi0 .ge. xi) then
  chi2 = chi0
  xi2  = xi0
  else
  chi1 = chi0
  xi1  = xi0
  endif
  if (abs(chi2-chi1)/abs(chi2) .le. eps0*5) exit
  end do

  ! compute the corresponding values of ap and appp
  call chebpw_eval(nints_xi,xi_ab,k,chebdata%xs,ap_vals0,chi1,apval)
  call chebpw_eval(nints_xi,xi_ab,k,chebdata%xs,appp_vals0,chi1,apppval)

vals0(i) = chi1
ap0(i)   = apval
appp0(i) = apppval
end do

ifsplit = 0

coefs0 = matmul(chebdata%u,vals0)
coefs0 = abs(coefs0) / abs(coefs0(1))
dd1    = maxval(coefs0(k-nn+1:k))

coefs0 = matmul(chebdata%u,ap0)
coefs0 = abs(coefs0) / abs(coefs0(1))
dd2    = maxval(coefs0(k-nn+1:k))

if (dd1 .gt. eps) ifsplit = 1
if (dd2 .gt. eps) ifsplit = 1

if (ifsplit .eq. 1) then

if (nstack+2 .gt. maxstack) then
print *,"prolates_exp1: stack overflow while construct expansion of chi"
print *,"gamma  = ",gamma
print *,"sigma1 = ",sigma1
print *,"sigma2 = ",sigma2
stop
endif

nstack          = nstack+1
stack(1,nstack) = (a+b)/2
stack(2,nstack) = b

nstack          = nstack+1
stack(1,nstack) = a
stack(2,nstack) = (a+b)/2
else

if (nints_chi .eq. maxints) then
print *,"prolates_exp1: too many intervals while constructing expansion of chi"
print *,"gamma  = ",gamma
print *,"sigma1 = ",sigma1
print *,"sigma2 = ",sigma2
stop
endif

nints_chi               = nints_chi+1
chi_ab(1,nints_chi)     = a
chi_ab(2,nints_chi)     = b
chi_vals(:,nints_chi)   = vals0
ap_vals(:,nints_chi)    = ap0
appp_vals(:,nints_chi)  = appp0

endif

end do

!call prini("nints_chi = ",nints_chi)
!call prin2("chi_ab   = ",chi_ab(:,1:nints_chi))
!xcall prin2("chi_vals = ",chi_vals(:,1:nints_chi))

!
!  Create coefficient expansions and determine how many coefficients we really
!  need to store
!

allocate(chi_coefs(k,nints_chi))
allocate(ap_coefs(k,nints_chi))
allocate(appp_coefs(k,nints_chi))

nn  = 0

do int=1,nints_chi
coefs0            = matmul(chebdata%u,chi_vals(:,int))
chi_coefs(:,int)  = coefs0
coefs0            = abs(coefs0) / abs(coefs0(1))

do i=1,k
if (coefs0(i) .gt. epscoefs) nn = max(nn,i)
end do

coefs0            = matmul(chebdata%u,ap_vals(:,int))
ap_coefs(:,int)   = coefs0
coefs0            = abs(coefs0) / abs(coefs0(1))
do i=1,k
if (coefs0(i) .gt. epscoefs) nn = max(nn,i)
end do

coefs0            = matmul(chebdata%u,appp_vals(:,int))
appp_coefs(:,int) = coefs0

end do

!nn = k

!
!  Copy out the data
!

expdata%gamma  = gamma
expdata%sigma1 = sigma1
expdata%sigma2 = sigma2
expdata%nints  = nints_chi
expdata%k      = nn
expdata%ncost  = 2*nints_chi + 3*nn*nints_chi

allocate(expdata%ab(2,nints_chi))
expdata%ab     = chi_ab(:,1:nints_chi)

allocate(expdata%chi_coefs(nn,nints_chi))
allocate(expdata%ap_coefs(nn,nints_chi))
allocate(expdata%appp_coefs(nn,nints_chi))

do int=1,nints_chi
expdata%chi_coefs(1:nn,int)  = chi_coefs(1:nn,int)
expdata%ap_coefs(1:nn,int)   = ap_coefs(1:nn,int)
expdata%appp_coefs(1:nn,int) = appp_coefs(1:nn,int)
end do

end subroutine


subroutine prolates_exp1_eval(expdata,sigma,chi,apval,apppval)
implicit double precision (a-h,o-z)
type(exp1_data)               :: expdata
!
!  Evaluate one of the expansions produced by prolates_exp1.
!
!  Input parameters:
!    expdata - the data describing the expansions
!    sigma - the value of the parameter at which to evaluate the expansions
!
!  Output parameters:
!    chi - the Sturm-Louville eigenvalue 
!    apval - the value of the derivative of the phase function at 0
!    apppval - the value of the third derivative of the phase function at 0
!

call chebpw_eval2(expdata%nints,expdata%ab,expdata%k,expdata%chi_coefs,sigma,chi)
call chebpw_eval2(expdata%nints,expdata%ab,expdata%k,expdata%ap_coefs,sigma,apval)
call chebpw_eval2(expdata%nints,expdata%ab,expdata%k,expdata%appp_coefs,sigma,apppval)

end subroutine


subroutine prolates_exp1_check(chebdata,exp1,ncheck,dacc)
implicit double precision (a-h,o-z)
type(chebexps_data)                        :: chebdata
type(exp1_data)                            :: exp1
!
!  Check the accuracy of the expansions produced by 
! 
!  Input parameters:
!     chebdata - the structure returned by chebexps which, among other things,
!       specifies the order of the piecewise Chebyshev expansions to use
!    exp1 - the expansion generated by prolates_exp1
!    ncheck - the number of random points at which to test the expansions
!
!  Output parameters:
!    dacc - the maximum relative error which was observed while test the
!      expansions
!
data pi / 3.14159265358979323846264338327950288d0 /

gamma  = exp1%gamma
sigma1 = exp1%sigma1
sigma2 = exp1%sigma2

dmax = 0
do i=1,ncheck
!call random_number(dd)
dd = (i-1.0d0)/(ncheck-1.0d0)
sigma = sigma1 + (sigma2-sigma1)*dd
xi    = sigma*gamma

call prolates_exp1_eval(exp1,sigma,chi,apval,apppval)
call prolates_phase3(chebdata,gamma,chi,aval0,apval0,appval0,apppval0)

xi0     = -2/pi*aval0
xi      = xi+1

errrel1 = abs((xi-xi0)/xi0)
errrel2 = abs((apval-apval0)/apval0)

dmax   = max(errrel1,dmax)
dmax   = max(errrel2,dmax)

end do

dacc = dmax

end subroutine


subroutine prolates_exp2(chebdata,sigma1,sigma2,gamma1,gamma2,exp2)
implicit double precision (a-h,o-z)
type(chebexps_data)                        :: chebdata
type(exp2_data), intent(out)               :: exp2
!
!  Construct an expansion of the Sturm-Louville eigenvalue and the related quantities
!  for a specified range of values of gamma and sigma.  The expansion is adaptive in 
!  sigma but not gamma!
!
!  Input parameters:
!    (sigma1,sigma2) - range of values of the parameter sigma for which the expansion
!       will hold, 
!
!       OR, IF SIGMA2 < 0, THEN THIS ROUTINE SETS SIGMA1 = 0 AND LETS SIGMA2 BE THE
!       LEAST VALUE SUCH THAT 
!
!          \lambda_m (gamma1) < epsilon0
!
!       WHERE:
!
!          - \lambda_1(gamma) > \lambda_2(gamma)  > ... are the eigenvalues of the
!          restricted Fourier operator
!
!          - m   = floor(sigma2*gamma1)
!
!          - epsilon0 is 2^(-52) = IEEE double machine zero
!
!    (gamma1,gamma2) - the range of values of the parameter gamma for which the
!      expansion will hold
!
!  Output parameters:
!    exp2 - a data structure describing the expansion
!
type(exp1_data)               :: exp1
double precision, allocatable :: coefs(:)

data pi / 3.14159265358979323846264338327950288d0 /

eps0  = epsilon(0.0d0)

k = chebdata%k
allocate(exp2%exps(k))
allocate(exp2%xscheb(k))
exp2%xscheb = chebdata%xs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! if (sigma2 .lt. 0 ) then
! sigma1 = 0
! eps    = 2.0d0**(-52)
! call prolates_dimension(gamma1,eps,ndim)
! ndim       = ndim+1.0d0
! sigma2     = ndim/gamma1
! endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!  Construct an expansion for each value of gamma
!

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,gamma)
!$OMP DO SCHEDULE(DYNAMIC)
do i=1,k
gamma   = (gamma2-gamma1)/2 * chebdata%xs(i) + (gamma2+gamma1)/2
call prolates_exp1(chebdata,sigma1,sigma2,gamma,exp2%exps(i))
end do
!$OMP END DO
!$OMP END PARALLEL


!
!  Compte the cost of the expansion 
!

ncost = k
do i=1,k
ncost = ncost + exp2%exps(i)%ncost
end do

!
!  Copy out the data in the structure
! 

exp2%sigma1 = sigma1
exp2%sigma2 = sigma2
exp2%gamma1 = gamma1
exp2%gamma2 = gamma2
exp2%k      = k
exp2%ncost  = ncost

end subroutine


subroutine prolates_exp2_eval(exp2,gamma,sigma,chi,apval,apppval)
implicit double precision (a-h,o-z)
type(exp2_data)               :: exp2
!
!  Evaluate one of the expansions produced by prolates_exp2.
!
!  Input parameters:
!    exp2 - a data structure returned by prolates_exp2 which describes the expansionm
!      of chi as a function of gamma and sigma
!    (gamma,xi) - the values of the parameters at which to evaluate chi
!
!  Output parameters:
!    chi - the Sturm-Louville eigenvalue corresponding to xi
!    apval - the value of the derivative of the phase function at 0
!    apppval - the value of the third derivative of the phase function at 0
!

double precision, allocatable :: chis(:),apvals(:),apppvals(:)

!sigma = xi/gamma
k      = exp2%k
gamma1 = exp2%gamma1
gamma2 = exp2%gamma2

allocate(chis(k),apvals(k),apppvals(k))

do i=1,k
!call chebpw_eval(expdata%nints,expdata%ab,expdata%k,expdata%xscheb,expdata%vals,sigma,chi)
call prolates_exp1_eval(exp2%exps(i),sigma,chis(i),apvals(i),apppvals(i))
end do

! perform the barycentric interpolation
call chebeval(gamma1,gamma2,k,exp2%xscheb,chis,gamma,chi)
call chebeval(gamma1,gamma2,k,exp2%xscheb,apvals,gamma,apval)
call chebeval(gamma1,gamma2,k,exp2%xscheb,apppvals,gamma,apppval)

end subroutine


subroutine prolates_exp2_eval0(exp2,gamma,sigma,chi)
implicit double precision (a-h,o-z)
type(exp2_data)               :: exp2
!
!  Evaluate chi only.
!
!  Input parameters:
!    exp2 - a data structure returned by prolates_exp2 which describes the expansionm
!      of chi as a function of gamma and sigma
!    (gamma,xi) - the values of the parameters at which to evaluate chi
!
!  Output parameters:
!    chi - the Sturm-Louville eigenvalue corresponding to xi
!

double precision :: chis(100)

k      = exp2%k
gamma1 = exp2%gamma1
gamma2 = exp2%gamma2

do i=1,k
call chebpw_eval2(exp2%exps(i)%nints,exp2%exps(i)%ab,exp2%exps(i)%k,exp2%exps(i)%chi_coefs,sigma,chis(i))
end do

! perform the barycentric interpolation
call chebeval(gamma1,gamma2,k,exp2%xscheb,chis,gamma,chi)

end subroutine


subroutine prolates_exp2_check(chebdata,exp2,ncheck,dacc)
implicit double precision (a-h,o-z)
type(chebexps_data)                        :: chebdata
type(exp2_data)                            :: exp2
!
!  Check the accuracy of the expansions constructed by prolates_exp2.
!
!  Input parameters:
!     chebdata - the structure returned by chebexps which, among other things,
!       specifies the order of the piecewise Chebyshev expansions to use
!    exp2 - the expansion generated by prolates_exp2
!    ncheck - the number of points at which to check the expansionx
!
!  Output parameters:
!    dacc - the maximum relative error observed 


!double precision, allocatable              :: errs(:,:)
double precision, allocatable              :: errs(:,:)
data pi / 3.14159265358979323846264338327950288d0 /

gamma1 = exp2%gamma1
gamma2 = exp2%gamma2
sigma1 = exp2%sigma1
sigma2 = exp2%sigma2

!
!  Check the accuracy of the expansion
!
dmax   = 0.0d0

allocate(errs(ncheck,ncheck))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,dd1,dd2,sigma,gamma,xi,chi,apval,apppval, &
!$OMP    aval0,apval0,appval0,apppval0,xi0,errrel1,errrel2,errrel3)
!$OMP DO SCHEDULE(dynamic)
do i=1,ncheck
do j=1,ncheck

!call random_number(dd1)
!call random_number(dd2)
dd1 = (i-1.0d0)/(ncheck-1.0d0)
dd2 = (j-1.0d0)/(ncheck-1.0d0)

sigma = sigma1 + (sigma2-sigma1)*dd1
gamma = gamma1 + (gamma2-gamma1)*dd2
xi    = sigma*gamma

call prolates_exp2_eval(exp2,gamma,sigma,chi,apval,apppval)
call prolates_phase3(chebdata,gamma,chi,aval0,apval0,appval0,apppval0)

xi0       = -2/pi*aval0
xi        = xi+1

errrel1   = abs((xi-xi0)/xi0)
errrel2   = abs((apval-apval0)/apval0)
errs(i,j) = max(errrel1,errrel2)

! if (errs(i,j) .gt. 1.0d-14) then
! print *,gamma,sigma,errrel1,errrel2
! endif

end do
end do

!$OMP END DO
!$OMP END PARALLEL

dacc = maxval(errs)

end subroutine


subroutine prolates_exp3(chebdata,nints,ab,sigma1,sigma2,exp3)
implicit double precision (a-h,o-z)
type(chebexps_data)                            :: chebdata
type(exp3_data), intent(out)                   :: exp3
double precision, allocatable, intent(inout)   :: ab(:,:)
!
!  Construct expansions of the Sturm-Liouville eigenvalue and the related
!  parameters as functions of gamma and sigma.    The expansions are
!  adaptive in both variables.
!
!  An initial partition partition for the gamma variable is specified 
!  and it is adaptively refined until the resulting expansions are 
!  sufficiently accurate.
!
!  Input parameters:
!    chebdata - the structure returned by chebexps which, among other things,
!      specifies the order of the piecewise Chebyshev expansions to use
!    (nints,ab) - the initial list of intervals for the gamma variable
!    (sigma1,sigma2) - the range of values of sigma for which the
!      expansion will hold
!    
!
!  Output parameters:
!    exp3 - a data structure describing the 
!


double precision, allocatable :: ab0(:,:)
double precision, allocatable :: stack(:,:)
type(exp2_data), allocatable  :: exps(:)
type(exp2_data)               :: exp2

eps0     = epsilon(0.0d0)
eps      = 1.0d-10
if (eps0 .lt. 1.0d-16) eps = 1.0d-14
if (eps0 .lt. 1.0d-30) eps = 1.0d-20

k        = chebdata%k
maxstack = 1000
maxints  = 1000
ncheck   = 100
ncost    = 0

allocate(stack(2,maxstack), exps(maxints), ab0(2,maxints))

nints0     = 0 
nstack     = 0
do i=nints,1,-1
nstack          = nstack+1
a               = ab(1,i)
b               = ab(2,i)
stack(1,nstack) = ab(1,i)
stack(2,nstack) = ab(2,i)
end do

do while (nstack > 0 )
a      = stack(1,nstack)
b      = stack(2,nstack)

!print *,"exp2: ",a,b
call elapsed(t1)
call prolates_exp2(chebdata,sigma1,sigma2,a,b,exp2)
call prolates_exp2_check(chebdata,exp2,ncheck,dacc)
call elapsed(t2)

ifsplit = 0
if (dacc .gt. eps) ifsplit = 1
print *,a,b,dacc,ifsplit,t2-t1,exp2%ncost

if (ifsplit .eq. 1) then
if (nstack+2 .gt. maxstack) then
print *,"prolates_exp3: stack overflow"
stop
endif

c               = (a+b)/2
stack(1,nstack) = c
stack(2,nstack) = b

nstack = nstack+1
stack(1,nstack) = a
stack(2,nstack) = c

else

if (nints0 .ge. maxints) then
print *,"prolates_exp3: too many intervals"
stop
endif

nints0        = nints0+1
ab0(1,nints0) = a
ab0(2,nints0) = b
exps(nints0)  = exp2
nstack        = nstack-1

endif

end do

!
!  Calculate the cost of the expansions
!
ncost = 2*nints
do int=1,nints0
ncost = ncost + exps(int)%ncost
end do

call prini("nints0 = ",nints0)
call prin2("ab0 = ",ab0(:,1:nints0))
call prini("ncost = ",ncost)
dd = ncost*8.0d0/(1024.0d0*1024.0d0)
call prin2("size (in mb) = ",dd)


!
!  Copy the data out 
!

exp3%gamma1 = ab(1,1)
exp3%gamma2 = ab(2,nints)
exp3%sigma1 = sigma1
exp3%sigma2 = sigma2
exp3%nints  = nints0
exp3%ncost  = ncost

allocate(exp3%exps(nints0),exp3%ab(2,nints0))

exp3%ncost  = ncost

do i=1,nints0
exp3%exps(i) = exps(i)
exp3%ab(1,i) = ab0(1,i)
exp3%ab(2,i) = ab0(2,i)
end do

end subroutine


subroutine prolates_exp3_eval(exp3,gamma,sigma,chi,apval,apppval)
implicit double precision (a-h,o-z)
type(exp3_data)                   :: exp3
!
!  Evaluate the expansions constructed by prolates_exp3.
!
!  Input parameters:
!    (gamma,sigma) - the values of the input parameters
!
!  Output parameters:
!    chi - the Sturm-Liouville eigenvalue
!    apval - the value of the derivative of the phase function at 0
!    apppval - the value of the third derivative of the phase function at 0
!
!



nints = exp3%nints

do int=1,nints-1
b = exp3%ab(2,int)
if (int .eq. nints .OR. gamma .le. b) exit
end do

a = exp3%ab(1,int)
b = exp3%ab(2,int)
call prolates_exp2_eval(exp3%exps(int),gamma,sigma,chi,apval,apppval)

end subroutine


subroutine prolates_exp3_eval0(exp3,gamma,sigma,chi)
implicit double precision (a-h,o-z)
type(exp3_data)                   :: exp3
!
!  Evaluate chi only.
!
!  Input parameters:
!    (gamma,sigma) - the values of the input parameters
!
!  Output parameters:
!    chi - the Sturm-Liouville eigenvalue
!
!



nints = exp3%nints

do int=1,nints-1
b = exp3%ab(2,int)
if (int .eq. nints .OR. gamma .le. b) exit
end do

a = exp3%ab(1,int)
b = exp3%ab(2,int)
call prolates_exp2_eval0(exp3%exps(int),gamma,sigma,chi)

end subroutine


subroutine prolates_exp3_check(chebdata,exp3,ncheck,dacc)
implicit double precision (a-h,o-z)
type(chebexps_data)                            :: chebdata
type(exp3_data)                                :: exp3
!
!  Check the accuracy of the expansions constructed by prolates_exp3.
!
!  Input parameters:
!     chebdata - the structure returned by chebexps which, among other things,
!       specifies the order of the piecewise Chebyshev expansions to use
!    exp3 - the expansion generated by prolates_exp3
!    ncheck - the number of points at which to check the expansionx
!
!  Output parameters:
!    dacc - the maximum relative error observed 


double precision, allocatable              :: errs(:)
data pi / 3.14159265358979323846264338327950288d0 /

gamma1 = exp3%gamma1
gamma2 = exp3%gamma2
sigma1 = exp3%sigma1
sigma2 = exp3%sigma2

!
!  Check the accuracy of the expansion
!
dmax   = 0.0d0
allocate(errs(ncheck))

!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd1,dd2,sigma,gamma,xi,chi,apval,apppval, &
!!$OMP    aval0,apval0,appval0,apppval0,xi0,errrel1,errrel2)
!!$OMP DO SCHEDULE(dynamic)
do i=1,ncheck
call random_number(dd1)
call random_number(dd2)
sigma = sigma1 + (sigma2-sigma1)*dd1
gamma = gamma1 + (gamma2-gamma1)*dd2
xi    = sigma*gamma

call prolates_exp3_eval(exp3,gamma,sigma,chi,apval,apppval)
call prolates_phase3(chebdata,gamma,chi,aval0,apval0,appval0,apppval0)

xi0       = -2/pi*aval0-1
xi        = xi+1
xi0       = xi0+1

errrel1   = abs((xi-xi0)/xi0)
errrel2   = abs((apval-apval0)/apval0)
errs(i)   = max(errrel1,errrel2)

end do
!!$OMP END DO
!!$OMP END PARALLEL

dacc = maxval(errs)

end subroutine


subroutine prolates_exps_merge(nints,ab,nints2,ab2)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(inout)       :: ab(:,:)
double precision                                   :: ab2(:,:)
!
!  Merge two collections of discretization intervals.
!
!  Input parameters:
!    (nints,ab) - the first collection of discretization intervals
!    (nint2,ab2) - a second collection of discretization intervals
!
!  Output parameters:
!    (nints,ab) - the merged collection of intervals
!  
!
double precision, allocatable :: xs(:),xs2(:)

maxints = 10000
allocate(xs(maxints),xs2(maxints))

idx = 0
do int=1,nints
idx     = idx+1
xs(idx) = ab(1,int)

idx     = idx+1
xs(idx) = ab(2,int)
end do


do int=1,nints2
idx     = idx+1
xs(idx) = ab2(1,int)

idx     = idx+1
xs(idx) = ab2(2,int)
end do


call quicksort(idx,xs)

nn     = 1
xs2(1) = xs(1)

do i=2,idx
if (xs(i) .gt. xs2(nn)) then
nn      = nn+1
xs2(nn) = xs(i)
endif
end do

deallocate(ab)
nints = nn-1
allocate(ab(2,nints))
do i=1,nints
ab (1,i) = xs2(i)
ab (2,i) = xs2(i+1)
end do


end subroutine



end module
