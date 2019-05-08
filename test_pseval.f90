!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This file contains code for testing the accuracy and speed with which the
!  code in prolates_fast.f90 evaluates the functions ps_n(x,c) via comparison with
!  the Xiao-Rokhlin algorithm.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module test_pseval_files

use utils
use chebyshev
use tensor
use prolates
use prolates_window
use prolates_fast

type        phase_data
integer                       :: nints
double precision, allocatable :: ab(:,:),apcoefs(:,:),acoefs(:,:),appcoefs(:,:)
double precision              :: cc1,chi
end type    phase_data

type        xr_data
integer                       :: ncoefs
double precision, allocatable :: coefsps(:),coefsqs(:)
double precision              :: cc2,cc1
double complex                :: clambda
end type    xr_data

contains


subroutine pseval(c1,c2,ncs,netas,nts,errabs,dtime1,dtime2,dtime3,dtime4)
implicit double precision (a-h,o-z)

type(chebexps_data)                :: chebdata

double precision, allocatable      :: coefsps(:),coefsqs(:)
double precision, allocatable      :: ab(:,:),acoefs(:,:),apcoefs(:,:),appcoefs(:,:)
double complex                     :: clambda,ima

double precision, allocatable      :: vals(:,:,:),vals0(:,:,:),roots(:),ts(:,:,:),cs(:)
double precision, allocatable      :: dtimes1(:,:),dtimes2(:,:),dtimes3(:,:),dtimes4(:,:)
double precision, allocatable      :: etas(:,:)
double precision, allocatable      :: derrs(:)

type(phase_data), allocatable      :: phases(:,:)
type(xr_data), allocatable         :: xr(:,:)


call elapsed(totaltime1)

allocate(dtimes1(netas,ncs),dtimes2(netas,ncs))
allocate(dtimes3(netas,ncs),dtimes4(netas,ncs))

allocate(cs(ncs),etas(netas,ncs),ts(nts,netas,ncs))
allocate(vals(nts,netas,ncs),vals0(nts,netas,ncs))

allocate(phases(netas,ncs))
allocate(xr(netas,ncs))

write (*,"(2(D16.8,2X))",advance="no") c1,c2


!
!  Form the list of points to evaluate the function at
!

ima = (0.0d0,1.0d0)
pi  = acos(-1.0d0)
k   = 30
call chebexps(k,chebdata)

do ic=1,ncs
c    = c1 + (c2-c1) *(ic-1.0d0)/(ncs-1.0d0)

cs(ic) = c
eta1   = 200
eta2   = c

do ieta=1,netas
call random_number(dd)
eta = eta1 + (eta2-eta1) * dd
n   = eta
eta = n

etas(ieta,ic) = eta

call prolexps(c,n,chi,apval,appval,apppval)


a      =  0.0d0
b      =  1.0d0
call qprolates_find_roots(0,c,chi,a,b,nroots,roots)
if (nroots .gt. 0) then
b      = roots(1)
endif
b      = min(.99d0,b)


do it=1,nts
call random_number(dd)
ts(it,ieta,ic) =  a + (b-a) * dd
end do


end do
end do


!
!  Construct the phase functions; find the average time
!

do ic=1,ncs
do ieta=1,netas
c   = cs(ic)
eta = etas(ieta,ic)
n   = eta

call elapsed(t1)
call prolates_fast_phase(chebdata,c,n,phases(ieta,ic)%chi,phases(ieta,ic)%nints,&
  phases(ieta,ic)%ab,phases(ieta,ic)%acoefs,phases(ieta,ic)%apcoefs,phases(ieta,ic)%appcoefs)
call elapsed(t2)
dtimes1(ieta,ic) = t2-t1

call prolates_fast_legendre0(n,valp,derp,valq,derq)
call prolates_fast_evalder(k,phases(ieta,ic)%nints,phases(ieta,ic)%ab,&
  phases(ieta,ic)%acoefs,phases(ieta,ic)%apcoefs,phases(ieta,ic)%appcoefs, &
  0.0d0,valps0,derps0,valqs0,derqs0)

if ( mod(n,2) .eq. 0) then
cc1 = valp/valps0
else
cc1 = derp/derps0
endif

phases(ieta,ic)%cc1 = cc1
end do
end do

dtime1 = sum(dtimes1)/(ncs*netas)


! evaluate using the phase functions

do ic=1,ncs
do ieta=1,netas

c   = cs(ic)
eta = etas(ieta,ic)
n   = eta

call elapsed(t1)
do it=1,nts
call prolates_fast_eval0(k,phases(ieta,ic)%nints,&
  phases(ieta,ic)%ab,phases(ieta,ic)%acoefs,phases(ieta,ic)%apcoefs,ts(it,ieta,ic),vals(it,ieta,ic))
vals(it,ieta,ic) = vals(it,ieta,ic)*phases(ieta,ic)%cc1

end do
call elapsed(t2)
dtimes2(ieta,ic) = t2-t1
end do
end do

dtime2 = sum(dtimes2)/(nts*ncs*netas)
write (*,"(2(D16.8,2X))",advance="no") dtime1,dtime2



!
!  Build the coefficient expansions
!

do ic=1,ncs
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ieta,c,eta,n,valps0,derps0,valqs0,derqs0,t1,t2,&
!$OMP   valp,derp,valq,derq,chi) 
!$OMP DO SCHEDULE(DYNAMIC)
do ieta=1,netas

c   = cs(ic)
eta = etas(ieta,ic)
n   = eta

call elapsed(t1)
call prolates_integer(c,n,chi,xr(ieta,ic)%clambda,xr(ieta,ic)%ncoefs,xr(ieta,ic)%coefsps,  &
  xr(ieta,ic)%coefsqs)
call elapsed(t2)
dtimes3(ieta,ic) = t2-t1

call prolates_integer_evalder(n,xr(ieta,ic)%ncoefs,xr(ieta,ic)%coefsps,xr(ieta,ic)%coefsqs,&
  0.0d0,valps0,derps0,valqs0,derqs0)
call prolates_fast_legendre0(n,valp,derp,valq,derq)

if ( mod(n,2) .eq. 0) then
xr(ieta,ic)%cc1 = valp/valps0
else
xr(ieta,ic)%cc1 = derp/derps0
endif


end do
!$OMP END DO
!$OMP END PARALLEL

end do

dtime3 = sum(dtimes3)/(ncs*netas)



do ic=1,ncs

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ieta,c,eta,n,t1,it,t2) 
!$OMP DO SCHEDULE(DYNAMIC)
do ieta=1,netas
c   = cs(ic)
eta = etas(ieta,ic)
n   = eta

call elapsed(t1)
do it=1,nts
call prolates_integer_eval0(n,xr(ieta,ic)%ncoefs,xr(ieta,ic)%coefsps,ts(it,ieta,ic),vals0(it,ieta,ic))
vals0(it,ieta,ic) = vals0(it,ieta,ic)*xr(ieta,ic)%cc1
end do
call elapsed(t2)

dtimes4(ieta,ic) = t2-t1
end do
!$OMP END DO
!$OMP END PARALLEL

end do

dtime4 = sum(dtimes4)/(ncs*netas*nts)
errabs = maxval(abs(vals-vals0))


call elapsed(totaltime2)

write (*,"(4(D16.8,2X))") dtime3,dtime4,errabs,totaltime2-totaltime1



end subroutine



end module



program test_pseval
use utils
use chebyshev
use prolates
use prolates_window
use prolates_exps
use test_pseval_files
implicit double precision (a-h,o-z)

double precision, allocatable :: dtimes1(:),dtimes2(:),davgs1(:),davgs2(:),derrs(:)
double precision, allocatable :: cs(:,:)

! 100 * 100 * 100 = 1,000,000 evaluations per row of table

ncs    = 100
netas  = 100
nts    = 100

call prolexps_ncds(nintscd)

allocate(dtimes1(nintscd),dtimes2(nintscd))
allocate(davgs1(nintscd),davgs2(nintscd))
allocate(derrs(nintscd))


do i=1,nintscd
call prolexps_intcd(i,c1,c2)
n1 = c1
n2 = c2


call pseval(c1,c2,ncs,netas,nts,errabs,dtime1,davg1,dtime2,davg2)

dtimes1(i) = dtime1
dtimes2(i) = dtime2
davgs1(i) = davg1
davgs2(i) = davg2
derrs(i)  = errabs

end do

!
!  Write the table to the disk
!


iw = 20
open(iw,FILE='pstable1.tex')


write(iw,'(A)') "\begin{tabular}{cccc}"
write(iw,'(A)') "\toprule"
write(iw,'(A)', advance='no') "        &"
write(iw,'(A)', advance='no') "Maximum absolute         &"
write(iw,'(A)', advance='no') "Average       &"
write(iw,'(A)'              ) "Average       \\"

write(iw,'(A)', advance='no') "Range of $\gamma$                         &"
write(iw,'(A)', advance='no') "error                    &"
write(iw,'(A)', advance='no') "evaluation time                     &"
write(iw,'(A)')               "evaluation time                    \\"

write(iw,'(A)', advance='no') "                         &"
write(iw,'(A)', advance='no') "                         &"
write(iw,'(A)', advance='no') "phase algorithm            &"
write(iw,'(A)')               "Xiao-Rokhlin             \\"

write(iw,'(A)')               "\midrule"

do i=1,nintscd

call prolexps_intcd(i,c1,c2)

n1 = c1
n2 = c2

errabs = derrs(i)
davg1  = davgs1(i)
davg2  = davgs2(i)
dtime1 = dtimes1(i)
dtime2 = dtimes2(i)

call write_table_integer_range(iw,n1,n2)
call write_table_next(iw)
call write_table_double(iw,errabs)
call write_table_next(iw)
call write_table_double(iw,davg1)
call write_table_next(iw)
call write_table_double(iw,davg2)

call write_table_nextline(iw)

end do

write (iw,'(A)') "\bottomrule"
write (iw,'(A)') "\end{tabular}"
close(iw)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11111111111


iw = 20
open(iw,FILE='pstable2.tex')


write(iw,'(A)') "\begin{tabular}{ccc}"
write(iw,'(A)') "\toprule"
write(iw,'(A)', advance='no') "                         &"
write(iw,'(A)', advance='no') "Average   &"
write(iw,'(A)'              ) "Average    \\"

write(iw,'(A)', advance='no') "Range of $\gamma$                             &"
write(iw,'(A)', advance='no') "precomp time                     &"
write(iw,'(A)')               "precomp time                    \\"

write(iw,'(A)', advance='no') "                         &"
write(iw,'(A)', advance='no') "phase algorithm          &"
write(iw,'(A)')               "Xiao-Rokhlin             \\"

write(iw,'(A)')               "\midrule"

do i=1,nintscd
call prolexps_intcd(i,c1,c2)
n1 = c1
n2 = c2

errabs = derrs(i)
davg1  = davgs1(i)
davg2  = davgs2(i)
dtime1 = dtimes1(i)
dtime2 = dtimes2(i)

call write_table_integer_range(iw,n1,n2)
call write_table_next(iw)
call write_table_double(iw,dtime1)
call write_table_next(iw)
call write_table_double(iw,dtime2)

call write_table_nextline(iw)

end do

write (iw,'(A)') "\bottomrule"
write (iw,'(A)') "\end{tabular}"
close(iw)



end program
