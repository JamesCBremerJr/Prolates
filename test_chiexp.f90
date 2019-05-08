!
!  Run extensive tests to check that the expansion of chi is correct via
!  comparison with the Xiao-Rokhlin algorithm
!

module test_chiexps_functions
use utils
use prolates
contains


subroutine check_chiexp(ncs,netas,c1,c2,dmax,dtime1,dtime2)
implicit double precision (a-h,p-z)
double precision              :: c1,c2
integer                       :: ncs,netas

!
!  Test the expansion of chi for a specified range of c's -- report the largest
!  relative error encountered and the maximum time per evaluation of t
!

double precision, allocatable :: chis(:,:),chis0(:,:),etas(:,:),cs(:)
double precision, allocatable :: times0(:,:)
double precision, allocatable :: coefsps(:),coefsqs(:)
double complex                :: clambda

allocate(chis(netas,ncs),chis0(netas,ncs))
allocate(times0(netas,ncs),etas(netas,ncs))
allocate(cs(ncs))


do ic=1,ncs
dd     = (ic-1.0d0)/(ncs-1.0d0)
c      = c1 + dd * (c2-c1)
cs(ic) = c

eta1   = 200
eta2   = c

do ieta=1,netas
dd            = (ieta-1.0d0)/(netas-1.0d0)
call random_number(dd)
eta           = eta1+(eta2-eta1)*dd
etas(ieta,ic) = eta
end do
end do


call elapsed(time1)


!
!  Compute all of the chi's using the expansion ... these c
!

call elapsed(t1)
do ic   = 1,ncs
do ieta = 1,netas
c   = cs(ic)
eta = etas(ieta,ic)
n   = eta
call prolexps(c,n,chis(ieta,ic),apval,appval,apppval)
end do
end do
call elapsed(t2)
dtime1 = (t2-t1)     / (ncs*netas)

write (*,"(3(D16.7,2X))",advance="no") c1,c2,dtime1

!
!  Do the same using the Xiao-Rokhlin algorithm ... record the times
!  so we can parallelize this because otherwise it is too slow
!

do ic   = 1,ncs
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c,ieta,eta,n,t1,t2,dd)
!$OMP DO SCHEDULE(DYNAMIC)
do ieta = 1,netas
c   = cs(ic)
eta = etas(ieta,ic)
n   = eta
call elapsed(t1)
call prolates_integer_eignvalue(c,n,chis0(ieta,ic))
call elapsed(t2)
times0(ieta,ic) = t2-t1
end do
!$OMP END DO
!$OMP END PARALLEL


end do






dtime2 = sum(times0) / (ncs*netas)
dmax   = maxval(abs(chis-chis0)/abs(chis0))


call elapsed(time2)

write (*,"(3(D16.7,2X))") dtime2,dmax,time2-time1



end subroutine


end module


program test_chiexps
use utils
use prolates
use prolates_window
use test_chiexps_functions

implicit double precision (a-h,o-z)

double precision, allocatable :: coefsps(:),coefsqs(:)
double complex                :: clambda

double precision, allocatable :: coefs(:),w(:)



! c =      525338.67735470936d0
! n =         100.00000000000000    
! call prolexps(c,n,chi,apval,appval,apppval)
! call prolates_integer_eignvalue(c,n,chi0)

! !chi0 = 105588023.156673170612d0

! print *,chi
! print *,chi0
! print *,(chi-chi0)/chi
! stop



ncs   = 500   ! number of values of c to sample in each range
netas = 500   ! number of values of n to sample in each range

!
!  Test the expansion of chi for various ranges of c and output a latex
!  showing the results
!

iw = 20
open(iw,FILE='chitable.tex')


write(iw,'(A)') "\begin{tabular}{cccc}"
write(iw,'(A)') "\toprule"
write(iw,'(A)', advance='no') "Range of $\gamma$        &"
write(iw,'(A)', advance='no') "Maximum relative         &"
write(iw,'(A)', advance='no') "Average time             &"
write(iw,'(A)'              ) "Average time             \\"

write(iw,'(A)', advance='no') "                         &"
write(iw,'(A)', advance='no') "difference               &"
write(iw,'(A)', advance='no') "expansion                &"
write(iw,'(A)')               "Xiao-Roklin             \\"

write(iw,'(A)')               "\midrule"


call prolexps_ncds(nintscd)


do intcd=1,nintscd
call prolexps_intcd(intcd,c1,c2)
n1 = c1
n2 = c2

call check_chiexp(ncs,netas,c1,c2, dmax,dtime1,dtime2)
call write_table_integer_range(iw,n1,n2)
call write_table_next(iw)
call write_table_double(iw,dmax)
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
