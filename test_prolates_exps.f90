program test_prolates_exps

use utils
use prolates
use prolates_phase
use prolates_exps

implicit double precision (a-h,o-z)

type(chebexps_data)            :: chebdata
type(exp1_data)                :: exp1, exp11
type(exp2_data)                :: exp2, exp22
type(exp3_data)                :: exp3, exp33
double precision, allocatable  :: ab(:,:), coefs(:)

k      = 30
call chebexps(k,chebdata)
pi     = acos(-1.0d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Check an expansion written to the disk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

iw     = 20
ncheck = 10000
open(iw,FILE="expansions.txt",STATUS='old')
call elapsed(t1)
call prolates_exp3_read(iw,exp33)
call elapsed(t2)
close(iw)
call prin2("read time = ",t2-t1)


call prolates_exp3_check(chebdata,exp33,ncheck,dacc)
call prin2("max relative error  = ",dacc)
gamma = 100
sigma = 0.5d0
call prolates_exp3_eval0(exp33,gamma,sigma,chi)


! do i=100,200
! gamma = i
! sigma = 0.5d0
! call elapsed(t1)
! call prolates_exp3_eval0(exp33,gamma,sigma,chi)
! call elapsed(t2)
! end do

stop



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Build an expansion in one variable and check the results
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! gamma   = 1000000.0d0
! ncheck  = 1000
! sigma1  = 0.0d0
! sigma2  = 1.0d0

! call prin2("gamma  = ",gamma)
! call prin2("sigma1 = ",sigma1)
! call prin2("sigma2 = ",sigma2)

! call elapsed(t1)
! call prolates_exp1(chebdata,sigma1,sigma2,gamma,exp1)
! call elapsed(t2)

! call prin2("exp1 time  = ",t2-t1)
! call prini("nints      = ",exp1%nints)
! call prini("k          = ",exp1%k)
! call prini("ncost      = ",exp1%ncost)

! call elapsed(t1)
! call prolates_exp1_check(chebdata,exp1,ncheck,dacc)
! call elapsed(t2)

! call prin2("check time =",t2-t1)
! call prin2("max relative error  = ",dacc)

! iw = 20
! open(iw,FILE="expansions.txt")
! call prolates_exp1_write(iw,exp1)
! close(iw)
! open(iw,FILE="expansions.txt",STATUS='old')

! call prolates_exp1_read(iw,exp11)
! call prolates_exp1_check(chebdata,exp11,ncheck,dacc)
! call prin2("max relative error  = ",dacc)
! close(iw)

! call prina("")
! call prina("")
! stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Build an expansion for a range of gammas and sigmas which is adaptive only in
!  sigma
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! gamma1  =     16
! gamma2  =     64

! sigma1  =     0.0d0
! sigma2  =     2.1d0
! ncheck  =     20

! call prin2("gamma1 = ",gamma1)
! call prin2("gamma2 = ",gamma2)

! call elapsed(t1)
! call prolates_exp2(chebdata,sigma1,sigma2,gamma1,gamma2,exp2)
! call elapsed(t2)

! call prin2("sigma1 = ",sigma1)
! call prin2("sigma2 = ",sigma2)

! call prolates_exp2_check(chebdata,exp2,ncheck,dacc)

! ncost = exp2%ncost
! dd    = ncost*8.0d0/(1024.0d0*1024.0d0)

! call prin2("exp2 time = ",t2-t1)
! call prini("expansion cost = ",ncost)
! call prin2("size of expansion (in MB) = ",dd)
! call prin2("max relative error  = ",dacc)

! iw = 20
! open(iw,FILE="expansions.txt")
! call prolates_exp2_write(iw,exp2)
! close(iw)
! open(iw,FILE="expansions.txt",STATUS='old')

! call prolates_exp2_read(iw,exp22)
! call prolates_exp2_check(chebdata,exp22,ncheck,dacc)
! call prin2("max relative error  = ",dacc)
! close(iw)


! call prina("")
! call prina("")

! stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Build an expansion for a range of gammas and sigmas which is adaptive in gamma and
!  sigma
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! sigma1 = 0.00d0
! sigma2 = 1.10d0
! ncheck = 1000

! ! Powers of 4 starting at 16
! nints  = 7
! allocate(ab(2,nints))
! do i=1,nints
! ab(1,i) = 4.0d0**(i+2)
! ab(2,i) = 4.0d0**(i+3)
! end do

! ! nints = 1
! ! ab(1,1) = 4**6
! ! ab(2,1) = 4**7

! call prini("before exp3, nints = ",nints)
! call prin2("before exp3, ab = ",ab(:,1:nints))

! call elapsed(t1)
! call prolates_exp3(chebdata,nints,ab,sigma1,sigma2,exp3)
! call elapsed(t2)

! ncost = exp3%ncost
! dd    = ncost*8/(1024.0d0*1024.0d0)

! call prin2("exp3 time = ",t2-t1)
! call prini("ncost     = ",ncost)
! call prin2("size of expansions (in MB) = ",dd)

! iw = 20
! open(iw,FILE="expansions.txt")
! call prolates_exp3_write(iw,exp3)
! close(iw)

! call elapsed(t1)
! call prolates_exp3_check(chebdata,exp3,ncheck,dacc)
! call elapsed(t2)

! call prin2("check time = ",t2-t1)
! call prin2("max relative error = ",dacc)

! open(iw,FILE="expansions.txt",STATUS='old')

! call prolates_exp3_read(iw,exp33)
! call prolates_exp3_check(chebdata,exp33,ncheck,dacc)
! call prin2("max relative error  = ",dacc)
! close(iw)

! call prina("")
! call prina("")



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Calculate the cost of storing chi only
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! nexps = exp3%nints
! ncost = k + 2*nexps

! do i=1,nexps
!   do j=1,k
!     nints = exp3%exps(i)%exps(j)%nints
!     kk     = exp3%exps(i)%exps(j)%k
!     allocate(coefs(kk))
!     nn    = 0 
!     do int=1,nints
!       coefs = exp3%exps(i)%exps(j)%chi_coefs(:,int)
!       coefs = abs(coefs)
!       coefs = coefs / coefs(1)
!       do l=1,kk
!         if (coefs(l) .gt. 1.0d-14) nn = max(nn,l) 
!       end do
!     end do
! !    print *,nn,kk
!     ncost = ncost + nn*nints
!     deallocate(coefs)
!   end do
! end do
! dd = ncost*8.0d0/(1024.0d0*1024.0d0)


! call prin2("cost of chi expansion (in MB) = ",dd)


end program

