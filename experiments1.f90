!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Some utility routines
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module experiments1_functions

use utils
use prolates
use prolates10
use prolates_phase
use prolates_exps

contains


subroutine test_accuracy1(exp3,gamma,sigma1,sigma2,errmax)
implicit double precision (a-h,o-z)
!
!  Check the error for one value of gamma
!
type(exp3_data)               :: exp3
real*10                       :: chi0,gamma0
real*10, allocatable          :: coefsps(:)

pi      = acos(-1.0d0)
errmax  = 0
nn      = 100
inormal = 0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,dd,sigma,n,gamma0,chi0,nn1,nn2,coefsps,chi)
!$OMP DO SCHEDULE(DYNAMIC)
do i=1,nn
call random_number(dd)
sigma  = sigma1 + dd*(sigma2-sigma1)
n      = sigma*gamma
sigma  = n/gamma
gamma0 = gamma

call prolates10_integer(inormal,gamma0,n,chi0,nn1,nn2,coefsps)
call prolates_exp3_eval0(exp3,gamma,sigma,chi)

errrel = abs(chi-chi0)/abs(chi0)
errmax = max(errrel,errmax)

end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine


subroutine test_accuracy2(exp3,gamma1,gamma2,sigma1,sigma2,errmax)
implicit double precision (a-h,o-z)
!
!  Check the error for a range of gammas ... sample 100,000 points
!
double precision, allocatable :: coefsps(:)
type(exp3_data)                :: exp3

ngamma = 100
errmax = 0 

do i=1,ngamma
dd    = (i-1.0d0)/(ngamma-1.0d0)
gamma = gamma1 + (gamma2-gamma1)*dd
call test_accuracy1(exp3,gamma,sigma1,sigma2,derr)
errmax = max(errmax,derr)
end do

end subroutine


subroutine test_time1(exp3,gamma1,gamma2,sigma1,sigma2,tphase)
implicit double precision (a-h,o-z)
!
!  Check the time
!
type(exp3_data)               :: exp3

pi      = acos(-1.0d0)
errmax  = 0
ngamma  = 100
nsigma  = 100
ntimes  = 100
inormal = 0

call elapsed(t1)
do j=1,ngamma
gamma = gamma1 + (gamma2-gamma1)*dd1
do i=1,nsigma

call random_number(dd1)
call random_number(dd2)

sigma = sigma1 + (sigma2-sigma1)*dd2
do l=1,ntimes
call prolates_exp3_eval0(exp3,gamma,sigma,chi)
end do

end do
end do

call elapsed(t2)
tphase = (t2-t1)/(ngamma*nsigma*ntimes+0.0d0)

end subroutine




subroutine test_time2(exp3,gamma1,gamma2,sigma1,sigma2,tphase,trokhlin)
implicit double precision (a-h,o-z)
!
!  Check the time
!
type(exp3_data)               :: exp3
double precision, allocatable :: coefsps(:)


pi      = acos(-1.0d0)
errmax  = 0
ngamma  = 100
nsigma  = 100
ntimes  = 100
inormal = 0

call elapsed(t1)
do j=1,ngamma
call random_number(dd1)
gamma = gamma1 + (gamma2-gamma1)*dd1
do i=1,nsigma

call random_number(dd2)
sigma = sigma1 + (sigma2-sigma1)*dd2
n     = gamma*sigma
sigma = n /gamma

do l=1,ntimes
call prolates_exp3_eval0(exp3,gamma,sigma,chi)
end do

end do
end do

call elapsed(t2)
tphase = (t2-t1)/(ngamma*nsigma*ntimes+0.0d0)


call elapsed(t1)
do j=1,ngamma
call random_number(dd1)
gamma = gamma1 + (gamma2-gamma1)*dd1
do i=1,nsigma

call random_number(dd2)

sigma = sigma1 + (sigma2-sigma1)*dd2
n     = gamma*sigma

!do l=1,ntimes
call prolates_integer(inormal,gamma,n,chi,nn1,nn2,coefsps)
!end do

end do
end do

call elapsed(t2)
trokhlin = (t2-t1)/(ngamma*nsigma+0.0d0)

end subroutine



end module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Main program
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program experiments1

use utils
use prolates
use prolates_phase
use prolates_exps
use experiments1_functions

implicit double precision (a-h,o-z)
! type(chebexps_data)            :: chebdata
! type(prolates_phase_data)      :: phase
! type(exp3_data)                :: exp3
! double precision, allocatable  :: coefsps(:), coefs(:)
! double precision, allocatable  :: xs(:), ys(:,:), gammas(:), sigmas(:)


type(chebexps_data)            :: chebdata
type(prolates_phase_data)      :: phase
type(exp1_data)                :: exp1
type(exp3_data)                :: exp3
double precision, allocatable  :: coefsps(:), coefs(:)
double precision, allocatable  :: xs(:), ys(:,:), gammas(:), sigmas(:)
double precision, allocatable  :: derrs(:,:),times1(:,:), times2(:,:)

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

call elapsed(t1)
do j=1,ntimes
call prolates_integer(inormal,gamma,n,chi0,n1,n2,coefsps)
end do
call elapsed(t2)
time1 = (t2-t1)/(ntimes+0.0d0)

nexps = exp3%nints
ncost = k + 2*nexps

do i=1,nexps
  do j=1,k
    nints = exp3%exps(i)%exps(j)%nints
    kk     = exp3%exps(i)%exps(j)%k
    print *,nexps,nints,kk
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
!    print *,nn,kk
    ncost = ncost + nn*nints
    deallocate(coefs)
  end do
end do
dd = ncost*8.0d0/(1024.0d0*1024.0d0)
call prini("ncost = ",ncost)
call prin2("cost of chi expansion (in MB) = ",dd)
stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Table giving time, comparison with OXR
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call elapsed(t3)
nsigma = 4
ngamma = 7

allocate(sigmas(nsigma+1),gammas(ngamma+1),range_labels(7))

sigmas(1)  = 0.00d0
sigmas(2)  = 0.25d0
sigmas(3)  = 0.50d0
sigmas(4)  = 0.75d0
sigmas(5)  = 1.00d0

gammas(1)  = 4**3
gammas(2)  = 4**4
gammas(3)  = 4**5
gammas(4)  = 4**6
gammas(5)  = 4**7
gammas(6)  = 4**8
gammas(7)  = 4**9
gammas(8)  = 4**10

! label1     = "100 to 500"
! label2     = "500 to 1,000"
! label3     = "1,000 to 5,000"
! label4     = "5,000 to 10,000"
! label5     = "10,000 to 50,000"
! label6     = "50,000 to 100,000"
! label7     = "100,000 to 500,000"
! label8     = "500,000 to 1,000,000"

label1     = "$4^3$ to $4^4$"
label2     = "$4^4$ to $4^5$"
label3     = "$4^5$ to $4^6$"
label4     = "$4^6$ to $4^7$"
label5     = "$4^7$ to $4^8$"
label6     = "$4^8$ to $4^9$"
label7     = "$4^9$ to $4^{10}$"

range_labels(1) = label1
range_labels(2) = label2
range_labels(3) = label3
range_labels(4) = label4
range_labels(5) = label5
range_labels(6) = label6
range_labels(7) = label7
!range_labels(8) = label8


iw = 1001
open(iw,FILE="table3.tex")
write(iw,"(A)") "\begin{tabular}{ccccc}"
write(iw,"(A)") "\toprule"
write(iw,"(A)",advance="no") "Range of $\gamma$ & Range of $\sigma$ & Average Time   & Average time    \\"
write(iw,"(A)") "                               &                   &  expansion     & Rokhlin, et. al.\\"
write(iw,"(A)") "\midrule"

allocate(times1(ngamma,nsigma),times2(ngamma,nsigma))

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,sigma1,sigma2,gamma1,gamma2,tphase,rokhlin)
!$OMP DO SCHEDULE(DYNAMIC)
do i=1,ngamma
do j=1,nsigma

sigma1 = sigmas(j)
sigma2 = sigmas(j+1)

gamma1 = gammas(i)
gamma2 = gammas(i+1)

call test_time2(exp3,gamma1,gamma2,sigma1,sigma2,tphase,trokhlin)
print *,"   ",gamma1,gamma2,sigma1,sigma2,tphase,trokhlin


times1(i,j) = tphase
times2(i,j) = trokhlin
end do
end do
!$OMP END DO
!$OMP END PARALLEL


do i=1,ngamma
do j=1,nsigma

tphase   = times1(i,j)
trokhlin = times2(i,j)

! Write to the table
if (j .eq. 1) then
write(iw,"(A,A)",advance = 'no') range_labels(i)," &"
else
write(iw,"(A)") " & "
endif

write(iw,"(A,F4.2,A,F4.2,A)", advance='no')  "$",sigmas(j)," - ",sigmas(j+1), " $ & "

write(iw,"(A)",advance='no') "$"
call write_table_double(iw,tphase)
write(iw,"(A)",advance='no') "$ &"

write(iw,"(A)",advance='no') "$"
call write_table_double(iw,trokhlin)
write(iw,"(A)",advance='no') "$ \\ "


write(iw,"(A)")  " \addlinespace[.125em]"

end do
write(iw,"(A)")  " \addlinespace[.25em]"
end do

write(iw,"(A)") "\bottomrule"
write(iw,"(A)") "\end{tabular}"

close(iw)

deallocate(gammas, sigmas, range_labels)

call elapsed(t4)

print *,"table103:  = ",t4-t3





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Create a table for accuracy
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nsigma = 4
ngamma = 7

allocate(sigmas(nsigma+1),gammas(ngamma+1),range_labels(8))

sigmas(1)  = 0.00d0
sigmas(2)  = 0.25d0
sigmas(3)  = 0.50d0
sigmas(4)  = 0.75d0
sigmas(5)  = 1.00d0

gammas(1)  = 4**3
gammas(2)  = 4**4
gammas(3)  = 4**5
gammas(4)  = 4**6
gammas(5)  = 4**7
gammas(6)  = 4**8
gammas(7)  = 4**9
gammas(8)  = 4**10

! label1     = "100 to 500"
! label2     = "500 to 1,000"
! label3     = "1,000 to 5,000"
! label4     = "5,000 to 10,000"
! label5     = "10,000 to 50,000"
! label6     = "50,000 to 100,000"
! label7     = "100,000 to 500,000"
! label8     = "500,000 to 1,000,000"

label1     = "$4^3$ to $4^4$"
label2     = "$4^4$ to $4^5$"
label3     = "$4^5$ to $4^6$"
label4     = "$4^6$ to $4^7$"
label5     = "$4^7$ to $4^8$"
label6     = "$4^8$ to $4^9$"
label7     = "$4^9$ to $4^{10}$"

range_labels(1) = label1
range_labels(2) = label2
range_labels(3) = label3
range_labels(4) = label4
range_labels(5) = label5
range_labels(6) = label6
range_labels(7) = label7
!range_labels(8) = label8



allocate(derrs(ncount,ncount))

iw = 1001
open(iw,FILE="table1.tex")
write(iw,"(A)") "\begin{tabular}{ccr@{\hspace{2em}}ccr}"
write(iw,"(A)") "\toprule"
write(iw,"(A)",advance="no") "Range of $\gamma$ & Range of $\sigma$ & Max relative error  &"
write(iw,"(A)",advance="no") "Range of $\gamma$ & Range of $\sigma$ & Max relative error  \\"
write(iw,"(A)") "\midrule"

do i=1,ceiling(ngamma/2.0d0)
do j=1,nsigma

sigma1 = sigmas(j)
sigma2 = sigmas(j+1)

! Find the first error
gamma1 = gammas(i)
gamma2 = gammas(i+1)


call test_accuracy2(exp3,gamma1,gamma2,sigma1,sigma2,errmax)
print *,gamma1,gamma2,sigma1,sigma2,errmax

if (j .eq. 1) then
write(iw,"(A,A)",advance = 'no') range_labels(i)," &"
else
write(iw,"(A)") " & "
endif


write(iw,"(A,F4.2,A,F4.2,A)", advance='no')  "$",sigmas(j)," - ",sigmas(j+1), " $ & "
call write_table_double(iw,errmax)
write(iw,"(A)", advance='no')  " &"


i2 = ceiling(ngamma/2.0d0)+i
if (i2 .le. ngamma) then
! Find the second error
gamma1 = gammas(i2)
gamma2 = gammas(i2+1)
call test_accuracy2(exp3,gamma1,gamma2,sigma1,sigma2,errmax)
print *,gamma1,gamma2,sigma1,sigma2,errmax


if (j .eq. 1) then
nn1 = gamma1
nn2 = gamma2
write(iw,"(A,A)",advance = 'no') range_labels(i2)," &"
else
write(iw,"(A)") " & "
endif

write(iw,"(A,F4.2,A,F4.2,A)", advance='no')  "$",sigmas(j)," - ",sigmas(j+1), " $ & "
call write_table_double(iw,errmax)


endif

write(iw,"(A)", advance='no')  "   \\"
!write(iw,"(A)")  "\addlinespace[.125em]"


end do
write(iw,"(A)")  "\addlinespace[.25em]"
end do


write(iw,"(A)") "\bottomrule"
write(iw,"(A)") "\end{tabular}"

close(iw)

deallocate(gammas, sigmas, derrs, range_labels)

stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Make a plots of chi as a function of sigma for gamma = 2 , gamma = 100, gamma = 10000, 
!  gamma = 10^6
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nn     = 100
nfuns  = 1
gamma  = 2.0d0
sigma1 = 0
sigma2 = 1.0d0

call prolates_exp1(chebdata,sigma1,sigma2,gamma,exp1)


allocate(xs(nn),ys(nfuns,nn))
do i=1,nn
dd    = (i-1.0d0)/(nn-1.0d0)
sigma = sigma1 + (sigma2-sigma1)*dd

call prolates_exp1_eval(exp1,sigma,chi,apval,apppval)

xs(i)   = sigma
ys(1,i) = chi

end do

call plot_functions("plotchi1.py","plotchi1.pdf",nfuns,"$\sigma$","$\chi$",0,0,   &
   sigma1,sigma2,0.0d0,10.0d00,nn,xs,ys,"",&
   "**","b-*r--*g-.*k:*mx*")


gamma  = 100.0d0
sigma1 = 0
sigma2 = 1.0d0

call prolates_exp1(chebdata,sigma1,sigma2,gamma,exp1)


do i=1,nn
dd    = (i-1.0d0)/(nn-1.0d0)
sigma = sigma1 + (sigma2-sigma1)*dd

call prolates_exp1_eval(exp1,sigma,chi,apval,apppval)

xs(i)   = sigma
ys(1,i) = chi

!print *,gamma,sigma,chi
end do

call plot_functions("plotchi2.py","plotchi2.pdf",nfuns,"$\sigma$","$\chi$",0,0,   &
   sigma1,sigma2,0.0d0,16000.0d0,nn,xs,ys,"",&
   "**","b-*r--*g-.*k:*mx*")

gamma  = 10000.0d0
sigma1 = 0
sigma2 = 1.0d0

call prolates_exp1(chebdata,sigma1,sigma2,gamma,exp1)


do i=1,nn
dd    = (i-1.0d0)/(nn-1.0d0)
sigma = sigma1 + (sigma2-sigma1)*dd

call prolates_exp1_eval(exp1,sigma,chi,apval,apppval)

xs(i)   = sigma
ys(1,i) = chi

!print *,gamma,sigma,chi
end do

call plot_functions("plotchi3.py","plotchi3.pdf",nfuns,"$\sigma$","$\chi$",0,0,   &
   sigma1,sigma2,0.0d0,1.6d8,nn,xs,ys,"",&
   "**","b-*r--*g-.*k:*mx*")


gamma  = 1.0d6
sigma1 = 0
sigma2 = 1.0d0

call prolates_exp1(chebdata,sigma1,sigma2,gamma,exp1)


do i=1,nn
dd    = (i-1.0d0)/(nn-1.0d0)
sigma = sigma1 + (sigma2-sigma1)*dd

call prolates_exp1_eval(exp1,sigma,chi,apval,apppval)

xs(i)   = sigma
ys(1,i) = chi

!print *,gamma,sigma,chi
end do

call plot_functions("plotchi4.py","plotchi4.pdf",nfuns,"$\sigma$","$\chi$",0,0,   &
   sigma1,sigma2,0.0d0,1.6d12,nn,xs,ys,"",&
   "**","b-*r--*g-.*k:*mx*")

deallocate(xs,ys)


! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !
! ! !  Plot time required to evaluate chi_\sigma(\gamma) as a function of sigma for several
! ! !  values of gamma
! ! !
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! ! nfuns = 4
! ! ! nn    = 200
! ! ! allocate(xs(nn), ys(nfuns,nn), gammas(nfuns) )

! ! ! gammas(1) = 100.0d0
! ! ! gammas(2) = 500.0d0
! ! ! gammas(3) = 1000.0d0
! ! ! gammas(4) = 5000.0d0

! ! ! sigma1 = 0.0d0
! ! ! sigma2 = 1.0d0
! ! ! ntimes = 500          ! number of times to average

! ! ! do i=1,nn


! ! ! dd    = (i-1.0d0)/(nn-1.0d0)
! ! ! sigma = sigma1 + (sigma2-sigma1)*dd

! ! ! do j=1,nfuns
! ! ! gamma = gammas(j)
! ! ! call elapsed(t1)
! ! ! do ii=1,ntimes
! ! ! call prolates_exp3_eval0(exp3,gamma,sigma,chi)
! ! ! end do
! ! ! call elapsed(t2)
! ! ! dtime   = (t2-t1)/(ntimes+0.0d0)
! ! ! dtime   = dtime*1000000
! ! ! ys(j,i) = dtime
! ! ! end do

! ! ! xs(i)   = sigma

! ! ! end do

! ! ! call plot_functions("plot1.py","plot1.pdf",nfuns,"$\sigma$","Time (in microseconds)",0,0,   &
! ! !    0.0d0,1.0d0,0.0d0,3.0d0,nn,xs,ys,"best",&
! ! !    "$\gamma=100$*$\gamma=500$*$\gamma=1,000$*$\gamma=5,000$*", &
! ! !    "b-*r--*g-.*k:*")

! ! ! deallocate(xs,ys,gammas)


! ! ! nfuns = 4
! ! ! nn    = 200
! ! ! allocate(xs(nn), ys(nfuns,nn), gammas(nfuns) )

! ! ! gammas(1) = 10000.0d0
! ! ! gammas(2) = 100000.0d0
! ! ! gammas(3) = 500000.0d0
! ! ! gammas(4) = 1000000.0d0

! ! ! sigma1 = 0.0d0
! ! ! sigma2 = 1.0d0
! ! ! ntimes = 500          ! number of times to average

! ! ! do i=1,nn


! ! ! dd    = (i-1.0d0)/(nn-1.0d0)
! ! ! sigma = sigma1 + (sigma2-sigma1)*dd

! ! ! do j=1,nfuns
! ! ! gamma = gammas(j)
! ! ! call elapsed(t1)
! ! ! do ii=1,ntimes
! ! ! call prolates_exp3_eval0(exp3,gamma,sigma,chi)
! ! ! end do
! ! ! call elapsed(t2)
! ! ! dtime   = (t2-t1)/(ntimes+0.0d0)
! ! ! dtime   = dtime*1000000
! ! ! ys(j,i) = dtime
! ! ! end do

! ! ! xs(i)   = sigma

! ! ! end do

! ! ! call plot_functions("plot3.py","plot3.pdf",nfuns,"$\sigma$","Time (in microseconds)",0,0,   &
! ! !    0.0d0,1.0d0,0.0d0,3.0d0,nn,xs,ys,"best",&
! ! !    "$\gamma=10^4$*$\gamma=10^5$*$\gamma = 5 \\times 10^5$*$\gamma=10^6$*", &
! ! !    "b-*r--*g-.*k:*")

! ! ! deallocate(xs,ys,gammas)


! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !
! ! !  Plot time required to evaluate chi_n(\gamma) as a function of gamma for several
! ! !  values of sigma
! ! !
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! ! nfuns = 4
! ! ! nn    = 200
! ! ! allocate(xs(nn), ys(nfuns,nn), sigmas(nfuns) )

! ! ! gamma1   = 100.0d0
! ! ! gamma2   = 1000000.0d0
! ! ! ntimes   = 500

! ! ! sigmas(1) = 0.00d0
! ! ! sigmas(2) = 0.10d0
! ! ! sigmas(3) = 0.250
! ! ! sigmas(4) = 0.33d0


! ! ! do i=1,nn

! ! ! dd    = (i-1.0d0)/(nn-1.0d0)
! ! ! gamma = gamma1  + (gamma2-gamma1)*dd

! ! ! do j=1,nfuns
! ! ! sigma = sigmas(j)
! ! ! call elapsed(t1)
! ! ! do ii=1,ntimes
! ! ! call prolates_exp3_eval0(exp3,gamma,sigma,chi)
! ! ! end do
! ! ! call elapsed(t2)
! ! ! dtime   = (t2-t1)/(ntimes+0.0d0)
! ! ! dtime   = dtime*1000000
! ! ! ys(j,i) = dtime
! ! ! end do

! ! ! xs(i)   = gamma

! ! ! end do

! ! ! call plot_functions("plot2.py","plot2.pdf",nfuns,"$\gamma$","Time (in microseconds)",1,0,   &
! ! !    100.0d0,1000000.0d0,0.0d0,3.0d0,nn,xs,ys,"best",&
! ! !    "$\sigma=0.00$*$\sigma=0.10$*$\sigma=0.25$*$\sigma=0.33$*", &
! ! !    "b-*r--*g-.*k:*")

! ! ! deallocate(xs,ys,sigmas)


! ! ! nfuns = 4
! ! ! nn    = 200
! ! ! allocate(xs(nn), ys(nfuns,nn), sigmas(nfuns) )

! ! ! gamma1   = 100.0d0
! ! ! gamma2   = 1000000.0d0
! ! ! ntimes   = 500

! ! ! sigmas(1) = 0.50d0
! ! ! sigmas(2) = 0.66d0
! ! ! sigmas(3) = 0.75d0
! ! ! sigmas(4) = 1.0d0


! ! ! do i=1,nn

! ! ! dd    = (i-1.0d0)/(nn-1.0d0)
! ! ! gamma = gamma1  + (gamma2-gamma1)*dd

! ! ! do j=1,nfuns
! ! ! sigma = sigmas(j)
! ! ! call elapsed(t1)
! ! ! do ii=1,ntimes
! ! ! call prolates_exp3_eval0(exp3,gamma,sigma,chi)
! ! ! end do
! ! ! call elapsed(t2)
! ! ! dtime   = (t2-t1)/(ntimes+0.0d0)
! ! ! dtime   = dtime*1000000
! ! ! ys(j,i) = dtime
! ! ! end do

! ! ! xs(i)   = gamma

! ! ! end do

! ! ! call plot_functions("plot4.py","plot4.pdf",nfuns,"$\gamma$","Time (in microseconds)",1,0,   &
! ! !    100.0d0,1000000.0d0,0.0d0,3.0d0,nn,xs,ys,"upper left",&
! ! !    "$\sigma=0.50$*$\sigma=0.66$*$\sigma=0.75$*$\sigma=1.00$*", &
! ! !    "b-*r--*g-.*k:*")

! ! ! deallocate(xs,ys,sigmas)




! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !
! ! !  Plots of accuracy
! ! !
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end program

