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

double precision, allocatable :: chis(:,:),chis0(:,:),etas(:,:)
double precision, allocatable :: times0(:,:)
double precision, allocatable :: coefsps(:),coefsqs(:)
double complex                :: clambda

allocate(chis(ncs,netas),chis0(ncs,netas))
allocate(times0(ncs,netas),etas(netas,ncs))

call elapsed(time1)
!
!  Compute all of the chi's using the expansion ... these c
!

call elapsed(t1)
do ic   = 1,ncs
c   = c1 + (c2-c1) * (ic-1.0d0) / (ncs-1.0d0)


!call prolexps_eta_range(c,eta1,eta2,eta3,eta4)
eta1 = 100
eta2 = c

do ieta = 1,netas
eta = eta1 + (eta2-eta1) * (ieta-1.0d0)/(netas-1.0d0)
call random_number(dd)
eta           = eta1 + (eta2-eta1) * dd
etas(ieta,ic) = eta
n             = eta
call prolexps(c,n,chis(ic,ieta),apval,appval,apppval)
end do

! do ieta = 1,netas
! !eta = eta3 + (eta4-eta3) * (ieta-1.0d0)/(netas-1.0d0)
! call random_number(dd)
! eta                 = eta3 + (eta4-eta3) * dd
! etas(ieta+netas,ic) = eta
! n                   = eta
! call prolexps(c,n,chis(ic,ieta+netas),apval,appval,apppval)
! end do

end do
call elapsed(t2)
dtime1 = (t2-t1)     / (ncs*netas)
print *,""

write (*,"(3(D16.7,2X))",advance="no") c1,c2,dtime1

!
!  Do the same using the Xiao-Rokhlin algorithm ... record the times
!  so we can parallelize this because otherwise it is too slow
!

do ic   = 1,ncs
c   = c1 + (c2-c1) * (ic-1.0d0) / (ncs-1.0d0)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ieta,eta,n,t1,t2)
!$OMP DO
do ieta = 1,netas
eta = etas(ieta,ic)
n   = eta
call elapsed(t1)
call prolates_integer_eignvalue(c,n,chis0(ic,ieta))
call elapsed(t2)


times0(ic,ieta) = t2-t1
end do
!$OMP END DO
!$OMP END PARALLEL


! !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ieta,eta,n,t1,t2)
! !$OMP DO
! do ieta = 1,netas
! !eta = eta3 + (eta4-eta3) * (ieta-1.0d0)/(netas-1.0d0)
! eta = etas(ieta+netas,ic)
! n   = eta
! call elapsed(t1)
! call prolates_integer_eignvalue(c,n,chis0(ic,ieta+netas))
! call elapsed(t2)
! times0(ic,ieta+netas) = t2-t1
! end do
! !$OMP END DO
! !$OMP END PARALLEL

end do


dtime2 = sum(times0) / (ncs*netas)
dmax   = maxval(abs(chis-chis0)/abs(chis0))


call elapsed(time2)

write (*,"(3(D16.7,2X))",advance="no") dtime2,dmax,time2-time1

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


ncs   = 100    ! number of values of c to sample in each range
netas = 100   ! number of values of n to sample in each range

!
!  Test the expansion of chi for various ranges of c and output a latex
!  showing the results
!

iw = 20
open(iw,FILE='chitable.tex')


write(iw,'(A)') "\begin{tabular}{rccc}"
write(iw,'(A)') "\toprule"
write(iw,'(A)', advance='no') "Range of $\gamma$        &"
write(iw,'(A)', advance='no') "Maximum relative         &"
write(iw,'(A)', advance='no') "Average time             &"
write(iw,'(A)'              ) "Average time             \\"

write(iw,'(A)', advance='no') "                         &"
write(iw,'(A)', advance='no') "error                    &"
write(iw,'(A)', advance='no') "our algorithm            &"
write(iw,'(A)')               "Xiao-Roklin             \\"

write(iw,'(A)')               "\midrule"


n1 = 500
n2 = 1000

do i=1,11
c1 = n1
c2 = n2
call check_chiexp(ncs,netas,c1,c2, dmax,dtime1,dtime2)

call write_table_integer_range(iw,n1,n2)
call write_table_next(iw)
call write_table_double(iw,dmax)
call write_table_next(iw)
call write_table_double(iw,dtime1)
call write_table_next(iw)
call write_table_double(iw,dtime2)
call write_table_nextline(iw)

n1 = n1*2
n2 = n2*2
end do


write (iw,'(A)') "\bottomrule"
write (iw,'(A)') "\end{tabular}"
close(iw)

!

end program
