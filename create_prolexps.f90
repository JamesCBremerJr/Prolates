module create_prolexps_functions
use utils

contains

subroutine write_integer(iw,varname,int)
implicit double precision (a-h,o-z)
integer                       :: int
character(len=*)              :: varname
write(iw,'("integer, parameter :: ",A12," = ",I8.8)') varname,int
end subroutine

subroutine write_array2(iw,arrayname,a)
implicit double precision (a-h,o-z)
double precision              :: a(:,:)
character(len=*)              :: arrayname

n1 = size(a,1)
n2 = size(a,2)

write(iw,'("double precision :: ",A12,"(",I3.3,",",I3.3,")")') arrayname,n1,n2
write(iw,'("data ",A12," /                    &")')  arrayname
write (iw,'(D24.16,",  &")') a(1:n1,1:n2-1)
write (iw,'(D24.16,",  &")') a(1:n1-1,n2)
write (iw,'(D24.16,"   /")') a(n1,n2)
write (iw,*)

end subroutine


subroutine write_array3(iw,arrayname,a)
implicit double precision (a-h,o-z)
double precision             :: a(:,:,:)
character(len=*)             :: arrayname

n1 = size(a,1)
n2 = size(a,2)
n3 = size(a,3)

write(iw,'("double precision :: ",A12,"(",I3.3,",",I3.3,",",I3.3,")")') arrayname,n1,n2,n3
write(iw,'("data ",A12," /                    &")')  arrayname

if (n3 .eq. 1) then
write (iw,'(D24.16,",   &")') a(1:n1,1:n2-1,n3)
write (iw,'(D24.16,",   &")') a(1:n1-1,n2,n3)
write (iw,'(D24.16,"    /")') a(n1,n2,n3)
else
write (iw,'(D24.16,",   &")') a(1:n1,1:n2,1:n3-1)
write (iw,'(D24.16,",   &")') a(1:n1,1:n2-1,n3)
write (iw,'(D24.16,",   &")') a(1:n1-1,n2,n3)
write (iw,'(D24.16,"    /")') a(n1,n2,n3)
endif

write (iw,*)

end subroutine


end module


program create_prolexps
use utils
use chebyshev
use tensor
use prolates_window
use prolates_exps
use create_prolexps_functions

implicit double precision (a-h,o-z)

type(chebexps_data)                 :: chebdata

double precision, allocatable       :: ab(:,:),cd(:,:),chicoefs(:,:,:),eta1coefs(:,:)
double precision, allocatable       :: eta2coefs(:,:),apcoefs(:,:,:),appcoefs(:,:,:)
double precision, allocatable       :: apppcoefs(:,:,:)

double precision, allocatable       :: ab2(:,:),cd2(:,:),chicoefs2(:,:,:),eta1coefs2(:,:)
double precision, allocatable       :: eta2coefs2(:,:),apcoefs2(:,:,:),appcoefs2(:,:,:)
double precision, allocatable       :: apppcoefs2(:,:,:)

eps0 = epsilon(0.0d0)

if (eps0 .lt. 1.d0-20) then
eps0 = 1.0d-20
elseif (eps0 .lt. 1.0d-17) then
eps = 1.0d-15
else
eps = 1.0d-13
endif

ic1     = 0
ic2     = 12

isubst  = 0
k       = 40
call chebexps(k,chebdata)

! 
!  Setup the header
!
iw = 20
open(iw,FILE='prolexps.f90')

write(iw,'(A)') "subroutine prolexps(c,n,chi,apval,appval,apppval)"
write(iw,'(A)') "implicit double precision (a-h,o-z)"
write(iw,'(A)') "double precision, intent(in)           :: c"
write(iw,'(A)') "integer, intent(in)                    :: n"
write(iw,'(A)') "double precision, intent(out)          :: chi, apval, appval, apppval"

!
!  Build the expansions
!

irange = 1
call prolates_exps_construct(eps,irange,ic1,ic2,chebdata,isubst,nintsab,ab,nintscd,cd, &
  chicoefs,apcoefs,appcoefs,apppcoefs,eta1coefs,eta2coefs)

irange = 2
call prolates_exps_construct(eps,irange,ic1,ic2,chebdata,isubst,nintsab2,ab2,nintscd2,cd2, &
  chicoefs2,apcoefs2,appcoefs2,apppcoefs2,eta1coefs2,eta2coefs2)

!
!  report on the sizes
!

dd1 = 8.0d0/(1024.0d0 * 1024.0d0)
dd1 = dd1 * k*(k+1)/2
dd1 = dd1 * nintsab * nintscd

dd2 = 8.0d0/(1024.0d0 * 1024.0d0)
dd2 = dd2 * k*(k+1)/2
dd2 = dd2 * nintsab2 * nintscd2

dd  = dd1*4 + dd2*4

call prind("size of each range 1 expansion (in MB)      = ",dd1)
call prind("size of each range 2 expansion (in MB)      = ",dd2)
call prind("total size (in MB)                          = ",dd)

call write_integer(iw,"k",k)
call write_integer(iw,"nintsab",nintsab)
call write_integer(iw,"nintsab2",nintsab2)
call write_integer(iw,"nintscd",nintscd)
call write_integer(iw,"nintscd2",nintscd2)


call write_array2(iw,"ab",ab)
call write_array2(iw,"cd",cd)
call write_array2(iw,"eta1coefs",eta1coefs)
call write_array2(iw,"eta2coefs",eta2coefs)
call write_array3(iw,"chicoefs",chicoefs)
call write_array3(iw,"apcoefs",apcoefs)
call write_array3(iw,"appcoefs",appcoefs)
call write_array3(iw,"apppcoefs",apppcoefs)


call write_array2(iw,"ab2",ab2)
call write_array2(iw,"cd2",cd2)
call write_array2(iw,"eta1coefs2",eta1coefs2)
call write_array2(iw,"eta2coefs2",eta2coefs2)
call write_array3(iw,"chicoefs2",chicoefs2)
call write_array3(iw,"apcoefs2",apcoefs2)
call write_array3(iw,"appcoefs2",appcoefs2)
call write_array3(iw,"apppcoefs2",apppcoefs2)

!
!  Write the code 
!

write(iw,'(A)') "call chebpw_eval2(nintscd,cd,k,eta1coefs,eta2coefs,c,eta1,eta2)"

write(iw,'(A)') "if( n .gt. eta2) then"
write(iw,'(A)') "call chebpw_eval2(nintscd2,cd2,k,eta1coefs2,eta2coefs2,c,eta1,eta2)"
write(iw,'(A)') "u = (n-eta1)/(eta2-eta1)"
write(iw,'(A,A)') "call tensorpw_eval4(k,nintsab2,ab2,nintscd2,cd2,chicoefs2,apcoefs2,", &
                  "appcoefs2,apppcoefs2,u,c,chi,apval,appval,apppval)"
write(iw,'(A)') "else"
write(iw,'(A)') "u = (n-eta1)/(eta2-eta1)"
write(iw,'(A,A)') "call tensorpw_eval4(k,nintsab,ab,nintscd,cd,chicoefs,apcoefs,", &
                  "appcoefs,apppcoefs,u,c,chi,apval,appval,apppval)"
write(iw,'(A)') "endif"
write(iw,'(A)') "end subroutine"
write(iw,*)     ""

!
!  Write the auxillary information subroutine which is used by the test code
!

write(iw,'(A)') "subroutine prolexps_eta_range(c,eta1,eta2,eta3,eta4)"
write(iw,'(A)') "implicit double precision (a-h,o-z)"
write(iw,'(A)') "double precision, intent(in)           :: c"
write(iw,'(A)') "double precision, intent(out)          :: eta1,eta2,eta3,eta4"

call write_integer(iw,"k",k)
call write_integer(iw,"nintscd",nintscd)
call write_integer(iw,"nintscd2",nintscd2)
call write_array2(iw,"cd",cd)
call write_array2(iw,"eta1coefs",eta1coefs)
call write_array2(iw,"eta2coefs",eta2coefs)
call write_array2(iw,"cd2",cd2)
call write_array2(iw,"eta1coefs2",eta1coefs2)
call write_array2(iw,"eta2coefs2",eta2coefs2)

write(iw,'(A)') "call chebpw_eval2(nintscd,cd,k,eta1coefs,eta2coefs,c,eta1,eta2)"
write(iw,'(A)') "call chebpw_eval2(nintscd2,cd2,k,eta1coefs2,eta2coefs2,c,eta3,eta4)"

write(iw,'(A)') "end subroutine"
write(iw,*)     ""


write(iw,'(A)') "subroutine prolexps_ncds(nintscd)"
write(iw,'(A)') "implicit double precision (a-h,o-z)"
call write_integer(iw,"nintscd0",nintscd)
write(iw,'(A)') "nintscd = nintscd0"
write(iw,'(A)') "end subroutine"
write(iw,*)     ""

write(iw,'(A)') "subroutine prolexps_intcd(intcd,c1,c2)"
write(iw,'(A)') "implicit double precision (a-h,o-z)"
call write_array2(iw,"cd",cd)
write(iw,'(A)') "c1 = cd(1,intcd)"
write(iw,'(A)') "c2 = cd(2,intcd)"
write(iw,'(A)') "end subroutine"
write(iw,*)     ""


close(iw)



!
!  Copy the auxillary routines into the file
!

call system("cat prolexps_aux.f90 >> prolexps.f90")


end program
