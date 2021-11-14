!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains miscellaneous utility routines which can be divided into the 
!  following categories:
!
!  Printing routines:
!  -----------------
!
!  The folllowing routines output the contents of arrays.  They have two principal
!  advantages over simply using "print": (1) the output of these routines is also
!  written into the file fort.13 for later consideration, (2) the output is formatted
!  in a consistent fashion.
!
!    prin2 - output an array of doubles with 7 digits displayed 
!    prind - output an array of doubles with 15 digits displayed
!    princ - output an array of complex numbers with 7 digits displayed
!    prinz - output an array of complex numbers with 15 digits displayed
!    prini - output an array of integers
!    prinl - output an array of long integers
!    prina - output a string
!
!  Plotting routines:
!  -----------------
!
!    plot_function - this is a quick utility routine which plots a function of
!     one variable using python and the matplotlib package
!
!    plot_2functions - this is a quick utility routine which plots 2 functions of
!     one variable on the same graph using python and the matplotlib package
!
!    plot_3functions - this is a quick utility routine which plots 3 functions of
!     one variable on the same graph using python and the matplotlib package
!
!    plot_function3d - this is a quick utility routine for producing a crude
!     3D plot of a functoin of two variables using python and matplotlib
!
!    plot_rectangles - produce a PDF plot of a collection of rectangles
!     using python and pyplot
!
!    plot_colored_rectangles - produce a PDF plot of a collection of rectangles
!     of varying colors rectangles using python and matplotlib
!
!    plot_matrix - produce a PDF plot of a matrix (in the style of MATLAB's
!     "imagesc" command) using python and matplotlib
!
!    plot_functions - produces publication quality PDF plots using
!      python and matplotlib but is a bit cumbsersome to use
!
!
!  Table generatation:
!  ------------------
!
!   This module contains the following extremely primitive routines for generating
!   laTeX tables:
!
!     latex_table_init - write the header for a laTeX table to a specified file
!     latex_table_double - write a double precision value in the table
!     latex_table_integer - write an integer value to the table
!     latex_table_10ton - write "10^n" into a table entry
!
!     latex_table_
!
!     latex_table_nextline - advance to the next line in a table
!     latex_table_finish - write out the trailing code for the laTeX table
!
!  Other routines:
!  ---------------
!
!    elapsed - return the wall clock time in seconds which has elapsed since some
!     arbitrary point in the past 
!    insort - sort an array of real numbers
!    insorti  - sort an array of integers
!    quicksort - sort an array of real numbers using the quicksort algorithm
!    quicksorti - sort an array of integers using the quicksort algorithm
!    iexclude - remove from a list of integers all entries which appear in
!      a seccond list
!    eye - return the (k,k) identity matrix
!    qrsolv - solve a system of linear equations with a real coefficient matrix
!      via a QR decomposition
!    cqrsolv - solve a system of linear equations with a complex coefficient matrix
!      via a QR decomposition
!    randperm - return a random permutation in S_n
!    equispaced_intervals - return an array representing a decomposition of the
!      interval [a,b] into equispaced subintervals
!    bisected_intervals - return an array representing the decomposition of the
!      interval [a,b] into bisected intervals
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module utils

interface prin2
   module procedure prin2_0
   module procedure prin2_1
   module procedure prin2_2
   module procedure prin2_3
end interface prin2

interface prind
   module procedure prind_0
   module procedure prind_1
   module procedure prind_2
end interface prind

interface princ
   module procedure princ_0
   module procedure princ_1
   module procedure princ_2
end interface princ

interface prinz
   module procedure prinz_0
   module procedure prinz_1
   module procedure prinz_2
end interface prinz

interface prini
   module procedure prini_0
   module procedure prini_1
   module procedure prini_2
end interface prini

interface prinl
   module procedure prinl_0
   module procedure prinl_1
   module procedure prinl_2
end interface prinl

contains

subroutine prin2_0(str,a)
implicit double precision (a-h,o-z)

double precision a
character(len=*), intent(in) :: str

print *,str
print "(6(2x,e15.7))",a

write (13,*) str
write (13,"(6(2x,e15.7))") a

end subroutine

subroutine prin2_1(str,a)
implicit double precision (a-h,o-z)

double precision, intent (in) :: a(:)
character(len=*), intent(in) :: str

print *,str
print "(6(2x,e15.7))",a

write (13,*) str
write (13,"(6(2x,e15.7))") a

end subroutine

subroutine prin2_2(str,a)
implicit double precision (a-h,o-z)

double precision, intent (in) :: a(:,:)
character(len=*), intent(in) :: str

print *,str
print "(6(2x,e15.7))",a

write (13,*) str
write (13,"(8(2x,e15.7))") a

end subroutine


subroutine prin2_3(str,a)
implicit double precision (a-h,o-z)

double precision, intent (in) :: a(:,:,:)
character(len=*), intent(in) :: str

print *,str
print "(6(2x,e15.7))",a

write (13,*) str
write (13,"(6(2x,e15.7))") a

end subroutine

subroutine prind_0(str,a)
implicit double precision (a-h,o-z)

double precision :: a
character(len=*), intent(in) :: str

print *,str
print "(2(2x,e24.16))",a

write (13,*) str
write (13,"(2(2x,e24.16))") a

end subroutine

subroutine prind_1(str,a)
implicit double precision (a-h,o-z)

double precision,intent(in) :: a(:)
character(len=*), intent(in) :: str

print *,str
print "(2(2x,e24.16))",a

write (13,*) str
write (13,"(2(2x,e24.16))") a

end subroutine


subroutine prind_2(str,a)
implicit double precision (a-h,o-z)

double precision,intent(in) :: a(:,:)
character(len=*), intent(in) :: str

print *,str
print "(2(2x,e24.16))",a

write (13,*) str
write (13,"(2(2x,e24.16))") a

end subroutine



subroutine princ_0(str,a)
implicit double precision (a-h,o-z)
double complex               :: a
character(len=*), intent(in) :: str

print *,str
print "(3(d15.7,',',d15.7,2X))",a

write (13,*) str
write (13,"(3(d15.7,',',d15.7,2X))") a

end subroutine


subroutine princ_1(str,a)
implicit double precision (a-h,o-z)

double complex,intent(in)    :: a(:)
character(len=*), intent(in) :: str

print *,str
print "(3(d15.7,',',d15.7,2X))",a

write (13,*) str
write (13,"(3(d15.7,',',d15.7,2X))") a

end subroutine


subroutine princ_2(str,a)
implicit double precision (a-h,o-z)

double complex,intent(in)    :: a(:,:)
character(len=*), intent(in) :: str

print *,str
print "(3(d15.7,',',d15.7,2X))",a

write (13,*) str
write (13,"(3(d15.7,',',d15.7,2X))") a

end subroutine


subroutine prinz_0(str,a)
implicit double precision (a-h,o-z)
double complex               :: a
character(len=*), intent(in) :: str

print *,str
print "(2(d24.15,',',d24.15,2X))",a

write (13,*) str
write (13,"(2(d24.15,',',d24.15,2X))") a

end subroutine


subroutine prinz_1(str,a)
implicit double precision (a-h,o-z)

double complex,intent(in)    :: a(:)
character(len=*), intent(in) :: str

print *,str
print "(2(d24.15,',',d24.15,2X))",a

write (13,*) str
write (13,"(2(d24.15,',',d24.15,2X))") a

end subroutine


subroutine prinz_2(str,a)
implicit double precision (a-h,o-z)

double complex,intent(in)    :: a(:,:)
character(len=*), intent(in) :: str

print *,str
print "(2(d24.15,',',d24.15,2X))",a

write (13,*) str
write (13,"(2(d24.15,',',d24.15,2X))") a

end subroutine


subroutine prini_0(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a
character(len=*), intent(in) :: str

print *,str
print "(8(2x,I9))",a

write (13,*) str
write (13,"(8(2x,I9))") a

end subroutine

subroutine prini_1(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a(:)
character(len=*), intent(in) :: str

print *,str
print "(8(2x,I8))",a

write (13,*) str
write (13,"(8(2x,I8))") a

end subroutine



subroutine prini_2(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a(:,:)
character(len=*), intent(in) :: str

print *,str
print "(8(2x,I8))",a

write (13,*) str
write (13,"(8(2x,I8))") a

end subroutine

subroutine prina(str)
implicit double precision (a-h,o-z)

character(len=*), intent(in) :: str

print *,str
write (13,*) str

end subroutine


subroutine prinl_0(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a
character(len=*), intent(in) :: str

print *,str
print "(6(2x,I13))",a

write (13,*) str
write (13,"(6(2x,I13))") a

end subroutine

subroutine prinl_1(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a(:)
character(len=*), intent(in) :: str

print *,str
print "(6(2x,I13))",a

write (13,*) str
write (13,"(6(2x,I13))") a

end subroutine

subroutine prinl_2(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a(:,:)
character(len=*), intent(in) :: str

print *,str
print "(6(2x,I13))",a

write (13,*) str
write (13,"(6(2x,I13))") a

end subroutine



subroutine elapsed(t)
implicit double precision (a-h,o-z)
integer*8 i,irate
real t1
call system_clock(i,irate)

dd = i
dd = dd/irate
t = dd
return
end subroutine




subroutine insort(k,a)
implicit double precision (a-h,o-z)
integer, intent(in)              :: k
double precision, intent (inout) :: a(k)

if (k .le. 1) return

do i=2,k
val=a(i)
j=i-1
do while (j .ge. 1 .AND. a(j) .gt. val) 
a(j+1)=a(j)
j=j-1
end do
a(j+1)=val
end do
end subroutine


subroutine insorti(k,ia)
implicit double precision (a-h,o-z)
integer, intent(in)              :: k
integer, intent (inout)          :: ia(k)

if (k .le. 1) return

do i=2,k
ival=ia(i)
j=i-1
do while (j .ge. 1 .AND. ia(j) .gt. ival) 
ia(j+1)=ia(j)
j=j-1
end do
ia(j+1)=ival
end do
end subroutine


subroutine insorti2(k,ia,idxs)
implicit double precision (a-h,o-z)
integer, intent(in)              :: k
integer, intent (inout)          :: ia(k),idxs(k)

if (k .le. 1) return

do i=2,k
ival   = ia(i)
idxval = idxs(i)
j=i-1
do while (j .ge. 1 .AND. ia(j) .gt. ival) 
ia(j+1)   = ia(j)
idxs(j+1) = idxs(j)
j=j-1
end do
ia(j+1)   = ival
idxs(j+1) = idxval
end do

end subroutine

subroutine insorti3(k,ia,ia2,ia3)
implicit double precision (a-h,o-z)
integer, intent(in)              :: k
integer, intent (inout)          :: ia(k),ia2(k),ia3(k)

if (k .le. 1) return

do i=2,k
ival   = ia(i)
ival2  = ia2(i)
ival3  = ia3(i)

j=i-1
do while (j .ge. 1 .AND. ia(j) .gt. ival) 
ia(j+1)   = ia(j)
ia2(j+1)  = ia2(j)
ia3(j+1)  = ia3(j)
j=j-1
end do
ia(j+1)   = ival
ia2(j+1)  = ival2
ia3(j+1)  = ival3
end do

end subroutine


subroutine insort3(k,ia,vals)
implicit double precision (a-h,o-z)
integer, intent(in)              :: k
integer, intent (inout)          :: ia(k)
double complex                   :: vals(k)
double complex                   :: val

if (k .le. 1) return

do i=2,k
ival   = ia(i)
val    = vals(i)

j=i-1
do while (j .ge. 1 .AND. ia(j) .gt. ival) 
ia(j+1)   = ia(j)
vals(j+1) = vals(j)
j=j-1
end do
ia(j+1)   = ival
vals(j+1) = val
end do
end subroutine


subroutine quicksort(n,vals)
implicit double precision (a-h,o-z)
dimension istack(2,20000)
dimension vals(1),idxs(1)
!
!       Sort a list of double precision numbers.
!
k        = 60

if (n .lt. k) then
call insort(n,vals)
return
endif

maxstack = 10000

m = 1
istack(1,1) = 1
istack(2,1) = n
!
 1000 continue
if (m .eq. 0) goto 1100
i1 = istack(1,m)
i2 = istack(2,m)
m=m-1
!
l = i2-i1+1
if (l .le. k) then
call insort(l,vals(i1))
goto 1000
endif
!
! Otherwise perform quicksort.
!
call quicksort01(vals,i1,i2,i3)
!
!       This should never happen, but just in case ...
!
! if (m+2 .ge. maxstack) then
! print *,"quicksort out of memory"
! stop
! endif
!
!  Make sure the smaller half is processed first to reduce storage
!  to O(logn).
!             
n1 = i3-i1+1
n2 = i2-i3
!
if (n2 .lt. n1) then
!
m = m+1
istack(1,m) = i1
istack(2,m) = i3
!
m = m+1
istack(1,m) = i3+1
istack(2,m) = i2
!
else
!
m = m+1
istack(1,m) = i3+1
istack(2,m) = i2
!
m = m+1
istack(1,m) = i1
istack(2,m) = i3
!
endif
!
goto 1000
 1100 continue 
end subroutine


        subroutine quicksort01(vals,i1,i2,i3)
        implicit double precision (a-h,o-z)
        dimension vals(1)
!
!       Randomly choose a pivot index.
!
!        call corrand3(1,r)
        call random_number(r)
        ipiv = i1+floor((i2-i1)*r)
!
!        ipiv = i1+(i2-i1)/2
!
        val  = vals(ipiv)
!
!       Swap the pivot element and the last element.
!
        vals(ipiv) = vals(i2)
        vals(i2)   = val
!
       i3 = i1
!
        do 1000 i=i1,i2-1
!
        if( vals(i) .lt. val) then
        d  = vals(i)
!
        vals(i)  = vals(i3)
        vals(i3) = d       
!
        i3=i3+1
        endif
!
 1000 continue
!
        dd = vals(i3)
        vals(i3) = vals(i2)
        vals(i2) = dd
!
        end subroutine



        subroutine quicksorti(n,ivals)
        implicit double precision (a-h,o-z)
        dimension istack(2,20000),ivals(1)
!
!       Sort a list of integers.
!
        if (n .lt. 60) then
        call insorti(n,ivals)
        return
        endif
!
        maxstack = 10000
        k        = 60
!
        nstack      = 1
        istack(1,1) = 1
        istack(2,1) = n
!
 1000 continue
        if (nstack .eq. 0) goto 1100
        i1 = istack(1,nstack)
        i2 = istack(2,nstack)
        nstack=nstack-1
!
!
        l = i2-i1+1
        if (l .le. k) then
        call insorti(l,ivals(i1))
        goto 1000
        endif
!
!       Otherwise perform quicksort.
!
        call quicksorti0(ivals,i1,i2,i3)
!
!       This should never happen, but just in case ...
!
        if (nstack+2 .ge. maxstack) then
        print *,"quicksort out of memory"
        stop
        endif
!
!       Make sure the smaller half is processed first to reduce storage
!       to O(logn).
!             
        n1 = i3-i1+1
        n2 = i2-(i3+1)+1
!
        if (n2 .lt. n1) then
!
        nstack = nstack+1
        istack(1,nstack) = i1
        istack(2,nstack) = i3
!
        nstack = nstack+1
        istack(1,nstack) = i3+1
        istack(2,nstack) = i2
!
        else
!
        nstack=nstack+1
        istack(1,nstack) = i3+1
        istack(2,nstack) = i2
!
        nstack=nstack+1
        istack(1,nstack) = i1
        istack(2,nstack) = i3
!
        endif
!
        goto 1000
 1100 continue
        end subroutine


        subroutine quicksorti0(ivals,i1,i2,i3)
        implicit double precision (a-h,o-z)
        dimension ivals(1)
!
!       Randomly choose a pivot index.
!
!        call corrand3(1,r)
        call random_number(r)
        ipiv = i1+floor((i2-i1)*r)
!
!        ipiv = i1+(i2-i1)/2
!
        ival  = ivals(ipiv)
!
!       Swap the pivot element and the last element.
!
        ivals(ipiv) = ivals(i2)
        ivals(i2)   = ival
!
        i3 = i1
!
        do 1000 i=i1,i2-1
        if( ivals(i) .lt. ival) then
        id  = ivals(i)
!
        ivals(i)  = ivals(i3)
        ivals(i3) = id
!
        i3=i3+1
        endif
 1000 continue
!
        id        = ivals(i3)
        ivals(i3) = ivals(i2)
        ivals(i2) = id
!
        end subroutine



        subroutine quicksorti2(n,ivals,ivals2)
        implicit double precision (a-h,o-z)
        dimension istack(2,20000),ivals(1),ivals2(1)
!
!       Sort a list of integers.
!
        if (n .lt. 60) then
        call insorti2(n,ivals,ivals2)
        return
        endif
!
        maxstack = 10000
        k        = 60
!
        nstack      = 1
        istack(1,1) = 1
        istack(2,1) = n
!
 1000 continue
        if (nstack .eq. 0) goto 1100
        i1 = istack(1,nstack)
        i2 = istack(2,nstack)
        nstack=nstack-1
!
!
        l = i2-i1+1
        if (l .le. k) then
        call insorti2(l,ivals(i1),ivals2(i1))
        goto 1000
        endif
!
!       Otherwise perform quicksort.
!
        call quicksorti20(ivals,ivals2,i1,i2,i3)
!
!       This should never happen, but just in case ...
!
        if (nstack+2 .ge. maxstack) then
        print *,"quicksort out of memory"
        stop
        endif
!
!       Make sure the smaller half is processed first to reduce storage
!       to O(logn).
!             
        n1 = i3-i1+1
        n2 = i2-(i3+1)+1
!
        if (n2 .lt. n1) then
!
        nstack = nstack+1
        istack(1,nstack) = i1
        istack(2,nstack) = i3
!
        nstack = nstack+1
        istack(1,nstack) = i3+1
        istack(2,nstack) = i2
!
        else
!
        nstack=nstack+1
        istack(1,nstack) = i3+1
        istack(2,nstack) = i2
!
        nstack=nstack+1
        istack(1,nstack) = i1
        istack(2,nstack) = i3
!
        endif
!
        goto 1000
 1100 continue
        end subroutine


        subroutine quicksorti20(ivals,ivals2,i1,i2,i3)
        implicit double precision (a-h,o-z)
        dimension ivals(1),ivals2(1)
!
!       Randomly choose a pivot index.
!
!        call corrand3(1,r)
        call random_number(r)
        ipiv = i1+floor((i2-i1)*r)
!
!        ipiv = i1+(i2-i1)/2
!
        ival  = ivals(ipiv)
        ival2 = ivals2(ipiv)
!
!       Swap the pivot element and the last element.
!
        ivals(ipiv)  = ivals(i2)
        ivals2(ipiv) = ivals2(i2)
        ivals(i2)    = ival
        ivals2(i2)   = ival2
!
        i3 = i1
!
        do 1000 i=i1,i2-1
        if( ivals(i) .lt. ival) then
        id  = ivals(i)
        id2 = ivals2(i)
!
        ivals(i)   = ivals(i3)
        ivals2(i)  = ivals2(i3)
        ivals(i3)  = id
        ivals2(i3) = id2
!
        i3=i3+1
        endif
 1000 continue
!
        id         = ivals(i3)
        id2        = ivals2(i3)
        ivals(i3)  = ivals(i2)
        ivals2(i3) = ivals2(i2)
        ivals(i2)  = id
        ivals2(i2) = id2
        
!
        end subroutine



        subroutine quicksorti3(n,ivals,ivals2,ivals3)
        implicit double precision (a-h,o-z)
        dimension istack(2,20000),ivals(1),ivals2(1),ivals3(1)
!
!       Sort a list of integers.
!
        if (n .lt. 60) then
        call insorti3(n,ivals,ivals2,ivals3)
        return
        endif
!
        maxstack = 10000
        k        = 60
!
        nstack      = 1
        istack(1,1) = 1
        istack(2,1) = n
!
 1000 continue
        if (nstack .eq. 0) goto 1100
        i1 = istack(1,nstack)
        i2 = istack(2,nstack)
        nstack=nstack-1
!
!
        l = i2-i1+1
        if (l .le. k) then
        call insorti3(l,ivals(i1),ivals2(i1),ivals3(i1))
        goto 1000
        endif
!
!       Otherwise perform quicksort.
!
        call quicksorti30(ivals,ivals2,ivals3,i1,i2,i3)
!
!       This should never happen, but just in case ...
!
        if (nstack+2 .ge. maxstack) then
        print *,"quicksort out of memory"
        stop
        endif
!
!       Make sure the smaller half is processed first to reduce storage
!       to O(logn).
!             
        n1 = i3-i1+1
        n2 = i2-(i3+1)+1
!
        if (n2 .lt. n1) then
!
        nstack = nstack+1
        istack(1,nstack) = i1
        istack(2,nstack) = i3
!
        nstack = nstack+1
        istack(1,nstack) = i3+1
        istack(2,nstack) = i2
!
        else
!
        nstack=nstack+1
        istack(1,nstack) = i3+1
        istack(2,nstack) = i2
!
        nstack=nstack+1
        istack(1,nstack) = i1
        istack(2,nstack) = i3
!
        endif
!
        goto 1000
 1100 continue
        end subroutine


        subroutine quicksorti30(ivals,ivals2,ivals3,i1,i2,i3)
        implicit double precision (a-h,o-z)
        dimension ivals(1),ivals2(1),ivals3(1)
!
!       Randomly choose a pivot index.
!
!        call corrand3(1,r)
        call random_number(r)
        ipiv = i1+floor((i2-i1)*r)
!
!        ipiv = i1+(i2-i1)/2
!
        ival  = ivals(ipiv)
        ival2 = ivals2(ipiv)
        ival3 = ivals3(ipiv)
!
!       Swap the pivot element and the last element.
!
        ivals(ipiv)  = ivals(i2)
        ivals2(ipiv) = ivals2(i2)
        ivals3(ipiv) = ivals3(i2)
        ivals(i2)    = ival
        ivals2(i2)   = ival2
        ivals3(i2)   = ival3
!
        i3 = i1
!
        do 1000 i=i1,i2-1
        if( ivals(i) .lt. ival) then
        id  = ivals(i)
        id2 = ivals2(i)
        id3 = ivals3(i)
!
        ivals(i)   = ivals(i3)
        ivals2(i)  = ivals2(i3)
        ivals3(i)  = ivals3(i3)
        ivals(i3)  = id
        ivals2(i3) = id2
        ivals3(i3) = id3
!
        i3=i3+1
        endif
 1000 continue
!
        id         = ivals(i3)
        id2        = ivals2(i3)
        id3        = ivals3(i3)
        ivals(i3)  = ivals(i2)
        ivals2(i3) = ivals2(i2)
        ivals3(i3) = ivals3(i2)
        ivals(i2)  = id
        ivals2(i2) = id2
        ivals3(i2) = i3
        
!
        end subroutine



        subroutine iremove(n,ia,m,ib)
        implicit double precision (a-h,o-z)
        dimension ia(n),ib(m)
!
!       Remove from the list ia of length n all integers appearing in 
!       the list ib of length m.  Both the list ia and the list ib
!       must be sorted before this call is made.  The results will
!       also be sorted.
!

!        call quicksorti(n,ia)
!        call quicksorti(m,ib)

        isrc = 1
        itar = 1
        ii   = 1
 1000 continue

        if (ii .gt. m)   goto 2000
        if (isrc .gt. n) goto 3000
        if (ia(isrc) .gt. ib(ii)) then
        ii=ii+1
        goto 1000
        endif

        if (ia(isrc) .lt. ib(ii)) then         
        ia(itar) = ia(isrc)
        itar=itar+1
        isrc=isrc+1
        goto 1000
        endif
        isrc=isrc+1
        goto 1000

 2000 continue
        if (isrc .gt. n) goto 3000
        ia(itar) = ia(isrc)
        itar=itar+1
        isrc=isrc+1
        goto 2000

 3000 continue
        n = itar-1
        return        
        
        end


subroutine qrsolv(a,n,rhs)
implicit double precision (a-h,o-z)

integer            :: n
double precision   :: a(n,n),rhs(n)

!
!  This subroutine uses a version of QR-decomposition to solve the equation
!  A x = b.  Both the input matrix a and the right-hand side are destroyed
!  b ythe routine.
!
!  Input parameters:
!    a - the (n,n) matrix of coefficients
!    n - an integer specifying the size of the system of 
!    rhs - a vector of length n speciying the rhs of the system
!
!  Output parameters:
!
!   rhs - upon return, the solution of the linear system


double precision :: aa(2),u(2,2)

! 
! transpose the input matrix a 
!

size22=0
do i=1,n
do j=1,i
d=a(j,i)
a(j,i)=a(i,j)
a(i,j)=d
size22=size22+a(j,i)**2
size22=size22+a(i,j)**2
end do
end do

!
!  Reduce to upper triangular 
!
do i=1,n-1
do j=n,i+1,-1
aa(1)=a(i,j-1)
aa(2)=a(i,j)

u22=-aa(1)
u12=aa(2)
d=u22**2+u12**2
if(d .lt. size22*1.0d-66) then
u(2,2)=1
u(1,2)=0
u(1,1)=1
u(2,1)=0
else
d=sqrt(d)
u(2,2)=u22/d
u(1,2)=u12/d
u(1,1)=-u(2,2)
u(2,1)=u(1,2)
endif

do ii=i,n
d1=u(1,1)*a(ii,j-1)+u(1,2)*a(ii,j)
d2=u(2,1)*a(ii,j-1)+u(2,2)*a(ii,j) 
a(ii,j-1)=d1
a(ii,j)=d2
end do

d1=u(1,1)*rhs(j-1)+u(1,2)*rhs(j)
d2=u(2,1)*rhs(j-1)+u(2,2)*rhs(j)
rhs(j-1)=d1
rhs(j)=d2
end do
end do

!
!  Apply the inverse of the triangular matrix
! 

rhs(n)=rhs(n)/a(n,n)
do i=n-1,1,-1
d=0
do j=n,i+1,-1
d=d+a(j,i)*rhs(j)
end do
rhs(i)=(rhs(i)-d)/a(i,i)
end do

return
end subroutine



subroutine cqrsolv(n,a,rhs)
implicit complex *16 (a-h,o-z)
double complex         :: a(n,n),u(2,2),rhs(n)
double precision       :: rcond
! 
!  This subroutine uses a version of QR-decomposition to solve
!  the user-supplied system of linear algebraic equations
!  Ax = y.   The matrix  a is destroyed in the process.
! 
!  Input parameters:
!    n - the dimensionality of the system being solved
!    a - the (n,n) complex-valued coefficient matrix WHICH WILL
!      BE DESTROYED BY THIS ROUUTINE
!    rhs - the right-hand side of the system to be solved, which
!      is overwritten by this subroutine
! 
!  Output parameters:
!    rhs - the solution of the linear system
!

double precision       ::  dmax,dmin,size22,dd,eps0
double complex         ::  aa(2)
integer, allocatable   ::  ipiv(:)

eps0 = epsilon(0.0d0)

if (eps0 .gt. 1.0d-17) then
allocate(ipiv(n))
!call zgesv(n,1,a,n,ipiv,rhs,n,info)
if (info .ne. 0) then
call prini("after zgesv, info = ",info)
stop
endif
return
endif

!
! Transpose the input matrix a
! 

size22=0
do i=1,n
do j=1,i
d=a(j,i)
a(j,i)=a(i,j)
a(i,j)=d
size22=size22+a(j,i)*conjg(a(j,i))
size22=size22+a(i,j)*conjg(a(i,j))
end do
end do

! 
! Eliminate the upper right triangle
! 
do i=1,n-1
! 
do j=n,i+1,-1
! 
aa(1)=a(i,j-1)
aa(2)=a(i,j)
u21=-aa(2)
u22=aa(1)
dd=u22*conjg(u22)+u21*conjg(u21)
if(dd .lt. size22*1.0d-66) then
u(2,2)=1
u(1,2)=0
u(1,1)=1
u(2,1)=0
else
dd=sqrt(dd)
u(2,2)=u22/dd
u(2,1)=u21/dd
u(1,1)=-conjg(u(2,2))
u(1,2)=conjg(u(2,1))
endif


do ii=i,n
d1=u(1,1)*a(ii,j-1)+u(1,2)*a(ii,j)
d2=u(2,1)*a(ii,j-1)+u(2,2)*a(ii,j)
a(ii,j-1) = d1
a(ii,j)   = d2
end do

d1=u(1,1)*rhs(j-1)+u(1,2)*rhs(j)
d2=u(2,1)*rhs(j-1)+u(2,2)*rhs(j)
rhs(j-1)=d1
rhs(j)=d2
end do
end do
! 
! Estimate the condition number
! 
dmax=-1
dmin=abs(a(1,1))
do i=1,n
if(dmax .lt. abs(a(i,i)) ) dmax=abs(a(i,i))
if(dmin .gt. abs(a(i,i)) ) dmin=abs(a(i,i))
end do
if(dmin .eq. 0) then
rcond=-1
return
endif
rcond=dmax/dmin

!
!  Apply the inverse of the triangular matrix to rhs
!
!call cqrtrin(a,rhs,n)

rhs(n)=rhs(n)/a(n,n)
do i=n-1,1,-1
d=0
do j=n,i+1,-1
d=d+a(j,i)*rhs(j)
end do
rhs(i)=(rhs(i)-d)/a(i,i)
end do

end subroutine


! function eye(n) result(a)
! implicit double precision (a-h,o-z)

! double precision, allocatable :: a(:,:)
! integer n

! !
! !  Return an identity matrix which is dimensioned (n,n).
! !

! allocate(a(n,n))
! a = 0
! do i=1,n
! a(i,i)=1.0d0
! end do

! end function


subroutine plot_rectangles(ifshow,filename,nrects,rects)
implicit double precision (a-h,o-z)
character(len=*)             :: filename
double precision             :: rects(4,nrects)

!
!  Produce a PDF plot of a collection of rectangles in the plane.
!
!  Input parameters:
!     filename - the name of the PDF file to produce
!     nrects - the number of rectangles to plot
!     rects - a (4,nrects) array which gives the coordinates of the
!        lower left corner of the rectangle and the upper-right corner
!        of the rectangle
!
!  Output parameters:
!     NONE
!

character(len=12)          :: scriptname
character(len=21)          :: command

call random_number(dd)
ii = dd * (10000)

write(scriptname,"(A,I5.5,A)") "plot",ii,".py"
write(command,"(A,A,I5.5,A)")  "python ","plot",ii,".py &"

iw = 20
open(iw,FILE=scriptname)
write(iw,"(A)") "import numpy as np"

write(iw,"(A)") "import matplotlib as mpl"
write(iw,"(A)") "mpl.use('Agg')"
write(iw,"(A)") "import matplotlib.pyplot as plt"
write(iw,"(A)") "import matplotlib.patches as patches"

!
!  Find the maximum extents
!
xmax = -1d300
xmin =  1d300
ymax = -1d300
ymin =  1d300

do i=1,nrects
x1 = rects(1,i)
y1 = rects(2,i)
x2 = rects(3,i)
y2 = rects(4,i)
xmax = max(xmax,x1)
xmin = min(xmin,x1)
xmax = max(xmax,x2)
xmin = min(xmin,x2)
ymax = max(xmax,y1)
ymin = min(xmin,y1)
ymax = max(xmax,y2)
ymin = min(xmin,y2)
end do


write(iw,"(A)") "fig, ax = plt.subplots()"

do i=1,nrects
x1 = rects(1,i)
y1 = rects(2,i)
x2 = rects(3,i)
y2 = rects(4,i)
dx = x2-x1
dy = y2-y1

write(iw,"('ax.add_patch(patches.Rectangle((',F24.15,',',F24.15,'),',F24.15,',',F24.15,',',A,'))')") &
  x1,y1,dx,dy,"facecolor = 'white', edgecolor = 'black'"

end do

write(iw,"(A,F24.16,',',F24.16,A)") "plt.xlim(",xmin,xmax,")"
write(iw,"(A,F24.16,',',F24.16,A)") "plt.ylim(",ymin,ymax,")"

write(iw,"(A,A,A)") 'fig.savefig("',filename,'")'



if (ifshow .eq. 1) then
write(iw,"(A,A,A)") 'fig.canvas.set_window_title("',filename,'")'
write(iw,"(A)") 'plt.show()'
endif

close(iw)

 call system(command)

! close(iw)

! call system("python plotscript.py")
! call system("rm -f plotscript.py")


end subroutine




subroutine plot_matrix(ifshow,filename,amatr)
implicit double precision (a-h,o-z)
character(len=*)             :: filename
double precision             :: amatr(:,:)

!
!  Produce a PDF plot of a matrix ala MATLAB's imagesc command.
!
!  Input parameters:
!     filename - the name of the PDF file to produce
!     amatr - an (n,m) real-valued matrix
!
!  Output parameters:
!     NONE
!

character(len=12)          :: scriptname
character(len=21)          :: command

n  = size(amatr,1)
m  = size(amatr,2)

call random_number(dd)
ii = dd * (10000)

write(scriptname,"(A,I5.5,A)") "plot",ii,".py"
write(command,"(A,A,I5.5,A)")  "python ","plot",ii,".py &"


iw = 20
open(iw,FILE=scriptname)
write(iw,"(A)") "import numpy as np"
write(iw,"(A)") "import matplotlib.pyplot as plt"

!
!  Find the maximum extents
!
xmax = -1d300
xmin =  1d300

do i=1,n
do j=1,m
xmax = max(xmax,amatr(i,j))
xmin = min(xmin,amatr(i,j))
end do
end do

write(iw,"(A)") "fig, ax = plt.subplots()"





write(iw,"(A,I5,',',I5,A)") "vals = np.zeros(",n,m,")"

do i=1,n
do j=1,m
write(iw,"(A,I5,',',I5,A,E24.16)") "vals[",i-1,j-1,"] = ",amatr(i,j)
end do
end do

write(iw,"(A)") "ax.imshow(grid)"

! do i=1,nrects
! x1 = rects(1,i)
! y1 = rects(2,i)
! x2 = rects(3,i)
! y2 = rects(4,i)
! dx = x2-x1
! dy = y2-y1

! write(iw,"('ax.add_patch(patches.Rectangle((',F24.15,',',F24.15,'),',F24.15,',',F24.15,',',A,'))')") &
!   x1,y1,dx,dy,"facecolor = 'white', edgecolor = 'black'"

! end do

! write(iw,"(A,F24.16,',',F24.16,A)") "plt.xlim(",xmin,xmax,")"
! write(iw,"(A,F24.16,',',F24.16,A)") "plt.ylim(",ymin,ymax,")"

! write(iw,"(A,A,A)") 'fig.savefig("',filename,'")'


! write(iw,"(A,A,A)") 'fig.canvas.set_window_title("',filename,'")'

if (ifshow .eq. 1) then
write(iw,"(A)") 'plt.show()'
endif

close(iw)

call system(command)


end subroutine


subroutine plot_functions(scriptname,filename,nfuns,xlabel,ylabel,iflogx,iflogy,x1_,x2_,y1_,y2_, &
   n,xs,ys,legend_loc,legends,styles)
implicit double precision (a-h,o-z)

character(len=*)               :: scriptname,filename,legends,legend_loc,styles
character(len=*)               :: xlabel,ylabel
double precision               :: xs(n),ys(nfuns,n)
character(len=:), allocatable  :: command, legend, style

!
!  Produce a python script which generates a PDF file containing a plot of one or more
!  functions of one variable suitable for inclusion in a paper.
!
!  Input parameters:
!    scriptname - the name of the python script
!    filename - the name of the PDF file to produce
!
!    xlabel - a label for the x-axis (no label will be present if the string is empty)
!    ylabel - a label for the y-axis (no label will be present if the string is empty)
!
!    iflogx - an integer parameter specifying whether the x-axis should be logarithm
!      scale or not
!    iflogy - an integer parameter specifying whether the y-axis should be logarithm
!      scale or not
!
!    (x1,x2,y1,y2) - extents of the axes ---- if x2 <= x1 then these will be set
!      automatically; likewise if y2 <= y1
!
!    legend_loc - a string specifying the location for the legend --- this can
!       be "upper right", "upper left", "lower right", etc  OR a blank string if 
!       no legend is to be included
!
!    legends - a string specifying the legend labels for each function ...
!       each should be terminated by an asterik "*" so "1*2*3*" specifies
!       the label 1 for the first function, 2 for the second, and 3 for the
!       third
!
!       this is ignored if legend_loc is blank
!
!    styles - a string specifying styles for each function ... separated as in
!       the legends string ... ignored if empty
!     
!    n - the number of point in the graph of the function to be specified
!    xs - the x-coordinates of the points to plot
!    ys - an (nfuns,n) array whose jth column gives the y-coordinates on
!        the graph of the jth function
!
!  Output parameters:
!    N/A
!


x1 = x1_
x2 = x2_
y1 = y1_
y2 = y2_

if (x2 .le. x1) then
x1 =  1d300
x2 = -1d300
do i=1,n
x1 = min(x1,xs(i))
x2 = max(x2,xs(i))
end do

endif

if (y2 .le. y1) then
y1 =  1d300
y2 = -1d300

do i=1,n
do j=1,nfuns
y1 = min(y1,ys(j,i))
y2 = max(y2,ys(j,i))
end do
end do

if (iflogy .eq. 0) then
y1 = y1*0.98d0
y2 = y2*1.02d0
endif

endif



iw = 1001
open(iw,FILE=scriptname)
write(iw,"(A)") "import numpy as np"
write(iw,"(A)") "import matplotlib as mpl"
write(iw,"(A)") "mpl.use('Agg')"
write(iw,"(A)") "import matplotlib.pyplot as plt"

write(iw,"(A,I5,A)") "xs  = np.zeros(",n,")"
write(iw,"(A,I5,A,I5,A)") "ys = np.zeros((",nfuns,",",n,"))"


do i=1,n
write(iw,"(A,I5,A,ES30.18E3)") "xs[",i-1,"]  = ",xs(i)
end do

do i=1,n
do j=1,nfuns
write(iw,"(A,I5,A,I5,A,ES30.18E3)") "ys[",j-1,",",i-1,"] = ",ys(j,i)
end do
end do

write(iw,"(A)") "fig, ax = plt.subplots()"

if (len(xlabel) .gt. 0) then
write (iw,"(A,A,A)") 'ax.set(xlabel="',xlabel,'")'
endif

if (len(ylabel) .gt. 0) then
write (iw,"(A,A,A)") 'ax.set(ylabel="',ylabel,'")'
endif



!allocate(legend(0),style(0))

idx1 = 1
idx2 = 1

idx3 = 1 
idx4 = 1

do j=1,nfuns


! find the legend string
if (len(legend_loc) .gt. 0) then
do while (legends(idx2:idx2) .ne. '*') 
idx2 = idx2+1
end do
ll = idx2-idx1
allocate(character(ll) :: legend )
legend(1:ll) = legends(idx1:idx2-1)
idx1 = idx2+1
idx2 = idx1
else
allocate(character(0) :: legend )
endif



! find the style string
if (len(styles) .gt. 0) then
do while (styles(idx4:idx4) .ne. '*') 
idx4 = idx4+1
end do
ll = idx4-idx3
allocate(character(ll) :: style )
style(1:ll) = styles(idx3:idx4-1)
idx3 = idx4+1
idx4 = idx3
else
allocate(character(0) :: style  )
endif

!print *,j,style," ",legend

write(iw,"(A,I5,A,A,A,A,A)") 'ax.plot(xs,ys[',j-1,',:],"',style,'",label="',legend,'")'

deallocate(style)
deallocate(legend)

end do


if (len(legend_loc) .gt. 0) then
write(iw,"(A,A,A)") 'plt.legend(loc="',legend_loc,'")'
endif

if (iflogx .eq. 1) then
write(iw,"(A)") 'plt.xscale("log")'
endif

if (iflogy .eq. 1) then
write(iw,"(A)") 'plt.yscale("log")'
endif

write(iw,"(A,E24.15,A,E25.15,A)") "plt.xlim([",x1,",",x2,"])"
write(iw,"(A,E24.15,A,E25.15,A)") "plt.ylim([",y1,",",y2,"])"
write(iw,"(A)") 'ax.grid()'
write(iw,"(A,A,A)") 'fig.savefig("',filename,'")'

close(iw)

allocate(character(7+len(scriptname)) :: command )
write(command,"(A,A)") "python ",scriptname
call system(command)

end subroutine



subroutine plot_colored_rectangles(filename,nrects,rects,icolors)
implicit double precision (a-h,o-z)
character(len=*)             :: filename
double precision             :: rects(4,nrects)
integer                      :: icolors(nrects)

!
!  Produce a PDF plot of a collection of rectangles in the plane.
!
!  Input parameters:
!     filename - the name of the PDF file to produce
!     nrects - the number of rectangles to plot
!     rects - a (4,nrects) array which gives the coordinates of the
!        lower left corner of the rectangle and the upper-right corner
!        of the rectangle
!
!  Output parameters:
!     NONE
!

iw = 20
open(iw,FILE="plotscript.py")
write(iw,"(A)") "import numpy as np"
write(iw,"(A)") "import matplotlib.pyplot as plt"
write(iw,"(A)") "import matplotlib.patches as patches"

!
!  Find the maximum extents
!
xmax = -1d300
xmin =  1d300
ymax = -1d300
ymin =  1d300

do i=1,nrects
x1 = rects(1,i)
y1 = rects(2,i)
x2 = rects(3,i)
y2 = rects(4,i)
xmax = max(xmax,x1)
xmin = min(xmin,x1)
xmax = max(xmax,x2)
xmin = min(xmin,x2)
ymax = max(xmax,y1)
ymin = min(xmin,y1)
ymax = max(xmax,y2)
ymin = min(xmin,y2)
end do

! write(iw,"(A,I5,A)") "xs = np.zeros(",n,")"
! write(iw,"(A,I5,A)") "ys = np.zeros(",n,")"

! do i=1,n
! write(iw,"(A,I5,A,E24.16)") "xs[",i-1,"] = ",xs(i)
! write(iw,"(A,I5,A,E24.16)") "ys[",i-1,"] = ",ys(i)
! end do

write(iw,"(A)") "fig, ax = plt.subplots()"

do i=1,nrects
x1 = rects(1,i)
y1 = rects(2,i)
x2 = rects(3,i)
y2 = rects(4,i)
dx = x2-x1
dy = y2-y1
icolor = icolors(i)

if (icolor .eq. 1) then
write(iw,"('ax.add_patch(patches.Rectangle((',F24.15,',',F24.15,'),',F24.15,',',F24.15,',',A,'))')") &
  x1,y1,dx,dy,"facecolor = 'white', edgecolor = 'black'"
else if(icolor .eq. 2) then
write(iw,"('ax.add_patch(patches.Rectangle((',F24.15,',',F24.15,'),',F24.15,',',F24.15,',',A,'))')") &
  x1,y1,dx,dy,"facecolor = 'blue', edgecolor = 'black'"
else if(icolor .eq. 3) then
write(iw,"('ax.add_patch(patches.Rectangle((',F24.15,',',F24.15,'),',F24.15,',',F24.15,',',A,'))')") &
  x1,y1,dx,dy,"facecolor = 'red', edgecolor = 'black'"
else if(icolor .eq. 4) then
write(iw,"('ax.add_patch(patches.Rectangle((',F24.15,',',F24.15,'),',F24.15,',',F24.15,',',A,'))')") &
  x1,y1,dx,dy,"facecolor = 'green', edgecolor = 'black'"
else if(icolor .eq. 5) then
write(iw,"('ax.add_patch(patches.Rectangle((',F24.15,',',F24.15,'),',F24.15,',',F24.15,',',A,'))')") &
  x1,y1,dx,dy,"facecolor = 'purple', edgecolor = 'black'"
endif

end do

write(iw,"(A,F24.16,',',F24.16,A)") "plt.xlim(",xmin,xmax,")"
write(iw,"(A,F24.16,',',F24.16,A)") "plt.ylim(",ymin,ymax,")"

write(iw,"(A,A,A)") 'fig.savefig("',filename,'")'

!write(iw,"(A,A,A)") 'fig.canvas.set_window_title("',filename,'")'
!write(iw,"(A)") 'plt.show()'

close(iw)

call system("python plotscript.py")
call system("rm -f plotscript.py")


end subroutine



subroutine plot_function(filename,xlabel,ylabel,n,xs,ys)
implicit double precision (a-h,o-z)

character(len=*)             :: filename,xlabel,ylabel
double precision             :: xs(n),ys(n)

!
!  Use Python and the pyplot package to produce a plot a function
!  of one variable.  The final output is a PDF file; a python
!  script which is executed using the system command is an intermediate
!  
!
!  Input parameters:
!    filename - the name of the PDF file to produce
!    xlabel - a label for the x-axis
!    ylabel - a label for the y-axis
!    n - the number of point in the graph of the function to be specified
!    xs - a vector of length n speciying the x coordinate of each point to plot
!    ys - a vector of length n speciying the y coordinate of each point to plot
!
!  Output parameters:
!    N/A
!
!


iw = 20
open(iw,FILE="plotscript.py")
write(iw,"(A)") "import numpy as np"
write(iw,"(A)") "import matplotlib as mpl"
write(iw,"(A)") "mpl.use('Agg')"
write(iw,"(A)") "import matplotlib.pyplot as plt"


write(iw,"(A,I5,A)") "xs = np.zeros(",n,")"
write(iw,"(A,I5,A)") "ys = np.zeros(",n,")"

do i=1,n
write(iw,"(A,I5,A,ES30.18E3)") "xs[",i-1,"] = ",xs(i)
write(iw,"(A,I5,A,ES30.18E3)") "ys[",i-1,"] = ",ys(i)
end do

write(iw,"(A)") "fig, ax = plt.subplots()"

if (len(xlabel) .gt. 0) then
write (iw,"(A,A,A)") 'ax.set(xlabel="',xlabel,'")'
endif

if (len(ylabel) .gt. 0) then
write (iw,"(A,A,A)") 'ax.set(ylabel="',ylabel,'")'
endif

write(iw,"(A)") "ax.plot(xs,ys)"
!write (iw,"(A)") "plt.ylim([-36,1])"
!write (iw,"(A)") "plt.xlim([0,1z])"

write(iw,"(A)") 'ax.grid()'
write(iw,"(A,A,A)") 'fig.savefig("',filename,'")'

! write(iw,"(A,A,A)") 'fig.canvas.set_window_title("',filename,'")'
! write(iw,"(A)") 'plt.show()'

close(iw)

call system("python plotscript.py")
! call system("rm -f plotscript.py")

end subroutine


subroutine plot_2functions(filename,xlabel,ylabel,legend1,legend2,n,xs,ys1,ys2)
implicit double precision (a-h,o-z)

character(len=*)             :: filename,xlabel,ylabel, legend1, legend2
double precision             :: xs(n),ys1(n),ys2(n)

!
!  Use Python and the pyplot package to produce a plot a function
!  of one variable.  The final output is a PDF file; a python
!  script which is executed using the system command is an intermediate
!  step.
!
!  Input parameters:
!    filename - the name of the PDF file to produce
!    xlabel - a label for the x-axis
!    ylabel - a label for the y-axis
!    legend1 - a label for the first functionm
!    legend2 - a label for the second function
!    n - the number of point in the graph of the function to be specified
!    xs - a vector of length n speciying the x coordinate of each point to plot
!    ys - a vector of length n speciying the y coordinate of each point to plot
!
!  Output parameters:
!    N/A
!
!


iw = 20
open(iw,FILE="plotscript.py")
write(iw,"(A)") "import numpy as np"
write(iw,"(A)") "import matplotlib as mpl"
write(iw,"(A)") "mpl.use('Agg')"
write(iw,"(A)") "import matplotlib.pyplot as plt"

write(iw,"(A,I5,A)") "xs  = np.zeros(",n,")"
write(iw,"(A,I5,A)") "ys1 = np.zeros(",n,")"
write(iw,"(A,I5,A)") "ys2 = np.zeros(",n,")"

do i=1,n
write(iw,"(A,I5,A,ES30.18E3)") "xs[",i-1,"]  = ",xs(i)
write(iw,"(A,I5,A,ES30.18E3)") "ys1[",i-1,"] = ",ys1(i)
write(iw,"(A,I5,A,ES30.18E3)") "ys2[",i-1,"] = ",ys2(i)

end do

write(iw,"(A)") "fig, ax = plt.subplots()"

if (len(xlabel) .gt. 0) then
write (iw,"(A,A,A)") 'ax.set(xlabel="',xlabel,'")'
endif

if (len(ylabel) .gt. 0) then
write (iw,"(A,A,A)") 'ax.set(ylabel="',ylabel,'")'
endif



if (len(legend1) .gt. 0) then
write(iw,"(A,A,A)") 'ax.plot(xs,ys1,label="',legend1,'")'
else
write(iw,"(A)") 'ax.plot(xs,ys1)'
endif

if (len(legend2) .gt. 0) then
write(iw,"(A,A,A)") 'ax.plot(xs,ys2,label="',legend2,'")'
else
write(iw,"(A)") "ax.plot(xs,ys2)"
endif

write(iw,"(A)") 'plt.legend(loc="upper left")'
write(iw,"(A)") 'ax.grid()'
write(iw,"(A,A,A)") 'fig.savefig("',filename,'")'


close(iw)

call system("python plotscript.py")
! call system("rm -f plotscript.py")

end subroutine


subroutine plot_3functions(filename,xlabel,ylabel,legend1,legend2,legend3,n,xs,ys1,ys2,ys3)
implicit double precision (a-h,o-z)

character(len=*)             :: filename,xlabel,ylabel,legend1,legend2,legend3
double precision             :: xs(n),ys1(n),ys2(n),ys3(n)

!
!  Use Python and the pyplot package to produce a plot a function
!  of one variable.  The final output is a PDF file; a python
!  script which is executed using the system command is an intermediate
!  step.
!
!  Input parameters:
!    filename - the name of the PDF file to produce
!    xlabel - a label for the x-axis
!    ylabel - a label for the y-axis
!    n - the number of point in the graph of the function to be specified
!    xs - a vector of length n speciying the x coordinate of each point to plot
!    ys - a vector of length n speciying the y coordinate of each point to plot
!
!  Output parameters:
!    N/A
!
!


iw = 20
open(iw,FILE="plotscript.py")
write(iw,"(A)") "import numpy as np"
write(iw,"(A)") "import matplotlib as mpl"
write(iw,"(A)") "mpl.use('Agg')"
write(iw,"(A)") "import matplotlib.pyplot as plt"

write(iw,"(A,I5,A)") "xs  = np.zeros(",n,")"
write(iw,"(A,I5,A)") "ys1 = np.zeros(",n,")"
write(iw,"(A,I5,A)") "ys2 = np.zeros(",n,")"
write(iw,"(A,I5,A)") "ys3 = np.zeros(",n,")"

do i=1,n
write(iw,"(A,I5,A,ES30.18E3)") "xs[",i-1,"]  = ",xs(i)
write(iw,"(A,I5,A,ES30.18E3)") "ys1[",i-1,"] = ",ys1(i)
write(iw,"(A,I5,A,ES30.18E3)") "ys2[",i-1,"] = ",ys2(i)
write(iw,"(A,I5,A,ES30.18E3)") "ys3[",i-1,"] = ",ys3(i)
end do

write(iw,"(A)") "fig, ax = plt.subplots()"

if (len(xlabel) .gt. 0) then
write (iw,"(A,A,A)") 'ax.set(xlabel="',xlabel,'")'
endif

if (len(ylabel) .gt. 0) then
write (iw,"(A,A,A)") 'ax.set(ylabel="',ylabel,'")'
endif


if (len(legend1) .gt. 0) then
write(iw,"(A,A,A)") 'ax.plot(xs,ys1,linestyle="solid",color="#1E90FF",label="',legend1,'")'
else
write(iw,"(A)") 'ax.plot(xs,ys1)'
endif

if (len(legend2) .gt. 0) then
write(iw,"(A,A,A)") 'ax.plot(xs,ys2,linestyle="dashed",color="orange",label="',legend2,'")'
else
write(iw,"(A)") "ax.plot(xs,ys2)"
endif

if (len(legend3) .gt. 0) then
write(iw,"(A,A,A)") 'ax.plot(xs,ys3,linestyle="dotted",color="red",label="',legend3,'")'
else
write(iw,"(A)") "ax.plot(xs,ys3)"
endif

write(iw,"(A)") 'plt.legend(loc="upper left")'
write(iw,"(A)") 'ax.grid()'
write(iw,"(A,A,A)") 'fig.savefig("',filename,'")'

close(iw)

call system("python plotscript.py")

end subroutine


subroutine plot_4functions(filename,xlabel,ylabel,legend1,legend2,legend3,legend4,n,xs,ys1,ys2,ys3,ys4)
implicit double precision (a-h,o-z)

character(len=*)             :: filename,xlabel,ylabel,legend1,legend2,legend3,legend4
double precision             :: xs(n),ys1(n),ys2(n),ys3(n),ys4(n)

!
!  Use Python and the pyplot package to produce a plot a function
!  of one variable.  The final output is a PDF file; a python
!  script which is executed using the system command is an intermediate
!  step.
!
!  Input parameters:
!    filename - the name of the PDF file to produce
!    xlabel - a label for the x-axis
!    ylabel - a label for the y-axis
!    n - the number of point in the graph of the function to be specified
!    xs - a vector of length n speciying the x coordinate of each point to plot
!    ys - a vector of length n speciying the y coordinate of each point to plot
!
!  Output parameters:
!    N/A
!
!


iw = 20
open(iw,FILE="plotscript.py")
write(iw,"(A)") "import numpy as np"
write(iw,"(A)") "import matplotlib as mpl"
write(iw,"(A)") "mpl.use('Agg')"
write(iw,"(A)") "import matplotlib.pyplot as plt"

write(iw,"(A,I5,A)") "xs  = np.zeros(",n,")"
write(iw,"(A,I5,A)") "ys1 = np.zeros(",n,")"
write(iw,"(A,I5,A)") "ys2 = np.zeros(",n,")"
write(iw,"(A,I5,A)") "ys3 = np.zeros(",n,")"
write(iw,"(A,I5,A)") "ys4 = np.zeros(",n,")"

do i=1,n
write(iw,"(A,I5,A,ES30.18E3)") "xs[",i-1,"]  = ",xs(i)
write(iw,"(A,I5,A,ES30.18E3)") "ys1[",i-1,"] = ",ys1(i)
write(iw,"(A,I5,A,ES30.18E3)") "ys2[",i-1,"] = ",ys2(i)
write(iw,"(A,I5,A,ES30.18E3)") "ys3[",i-1,"] = ",ys3(i)
write(iw,"(A,I5,A,ES30.18E3)") "ys4[",i-1,"] = ",ys4(i)
end do

write(iw,"(A)") "fig, ax = plt.subplots()"

if (len(xlabel) .gt. 0) then
write (iw,"(A,A,A)") 'ax.set(xlabel="',xlabel,'")'
endif

if (len(ylabel) .gt. 0) then
write (iw,"(A,A,A)") 'ax.set(ylabel="',ylabel,'")'
endif



if (len(legend1) .gt. 0) then
write(iw,"(A,A,A)") 'ax.plot(xs,ys1,label="',legend1,'")'
else
write(iw,"(A)") 'ax.plot(xs,ys1)'
endif

if (len(legend2) .gt. 0) then
write(iw,"(A,A,A)") 'ax.plot(xs,ys2,label="',legend2,'")'
else
write(iw,"(A)") "ax.plot(xs,ys2)"
endif

if (len(legend3) .gt. 0) then
write(iw,"(A,A,A)") 'ax.plot(xs,ys3,label="',legend3,'")'
else
write(iw,"(A)") "ax.plot(xs,ys3)"
endif

if (len(legend4) .gt. 0) then
write(iw,"(A,A,A)") 'ax.plot(xs,ys4,label="',legend4,'")'
else
write(iw,"(A)") "ax.plot(xs,ys4)"
endif

write(iw,"(A)") 'plt.legend(loc="upper left")'
!write (iw,"(A)") "plt.ylim([-36,1])"
!write (iw,"(A)") "plt.xlim([0,1z])"

write(iw,"(A)") 'ax.grid()'
write(iw,"(A,A,A)") 'fig.savefig("',filename,'")'

! write(iw,"(A,A,A)") 'fig.canvas.set_window_title("',filename,'")'
! write(iw,"(A)") 'plt.show()'

close(iw)

call system("python plotscript.py")
! call system("rm -f plotscript.py")

end subroutine


subroutine plot_points(ifshow,filename,n,xs,ys)
implicit double precision (a-h,o-z)

character(len=*)             :: filename
double precision             :: xs(n),ys(n)

!
!  Use Python and the pyplot package to produce a scatter plot.
!
!  Input parameters:
!    filename - the name of the PDF file to produce
!    n - the number of point in the graph of the function to be specified
!    xs - a vector of length n speciying the x coordinate of each point to plot
!    ys - a vector of length n speciying the y coordinate of each point to plot
!
!  Output parameters:
!    N/A
!
!

character(len=12)          :: scriptname
character(len=21)          :: command

call random_number(dd)
ii = dd * (10000)

write(scriptname,"(A,I5.5,A)") "plot",ii,".py"
write(command,"(A,A,I5.5,A)")  "python ","plot",ii,".py &"


iw = 20
open(iw,FILE=scriptname)
write(iw,"(A)") "import numpy as np"
write(iw,"(A)") "import matplotlib.pyplot as plt"


write(iw,"(A,I5,A)") "xs = np.zeros(",n,")"
write(iw,"(A,I5,A)") "ys = np.zeros(",n,")"

do i=1,n
write(iw,"(A,I5,A,E24.16)") "xs[",i-1,"] = ",xs(i)
write(iw,"(A,I5,A,E24.16)") "ys[",i-1,"] = ",ys(i)
end do

write(iw,"(A)") "fig, ax = plt.subplots()"
write(iw,"(A)") "ax.scatter(xs,ys,s=10)"
write(iw,"(A)") 'ax.grid()'
write(iw,"(A,A,A)") 'fig.savefig("',filename,'")'


if (ifshow .eq. 1) then
write(iw,"(A,A,A)") 'fig.canvas.set_window_title("',filename,'")'
write(iw,"(A)") 'plt.show()'
endif

close(iw)

 call system(command)

!call system("python plotscript.py")
!call system("rm -f plotscript.py")

end subroutine



subroutine plot_function3d(filename,n,xs,ys,zs)
implicit double precision (a-h,o-z)

character(len=*)             :: filename
double precision             :: xs(n),ys(n),zs(n)

!
!  Use Python and the pyplot package to produce a plot a function
!  of two variables.  The final output is a PDF file; a python
!  script which is executed using the system command is an intermediate
!
!  Input parameters:
!    filename - the name of the PDF file to produce
!    xlabel - a label for the x-axis
!    ylabel - a label for the y-axis
!    n - the number of point in the graph of the function to be specified
!    xs - a vector of length n speciying the x coordinate of each point to plot
!    ys - a vector of length n speciying the y coordinate of each point to plot
!
!  Output parameters:
!    N/A
!
!

character(len=12)          :: scriptname
character(len=21)          :: command

call random_number(dd)
ii = dd * (10000)

write(scriptname,"(A,I5.5,A)") "plot",ii,".py"
write(command,"(A,A,I5.5,A)")  "python ","plot",ii,".py &"


iw = 20
open(iw,FILE=scriptname)
write(iw,"(A)") "import numpy as np"

! write(iw,"(A)") "import matplotlib as mpl"
! write(iw,"(A)") "mpl.use('Agg')"

write(iw,"(A)") "import matplotlib.pyplot as plt"
write(iw,"(A)") "from mpl_toolkits.mplot3d import Axes3D"

! fig = plt.figure()
! ax = plt.axes(projection='3d')
! ax.plot_wireframe(X, Y, Z, color='black')
! ax.set_title('wireframe');

write(iw,"(A,I5,A)") "xs = np.zeros(",n,")"
write(iw,"(A,I5,A)") "ys = np.zeros(",n,")"
write(iw,"(A,I5,A)") "zs = np.zeros(",n,")"

do i=1,n
write(iw,"(A,I5,A,E24.16)") "xs[",i-1,"] = ",xs(i)
write(iw,"(A,I5,A,E24.16)") "ys[",i-1,"] = ",ys(i)
write(iw,"(A,I5,A,E24.16)") "zs[",i-1,"] = ",zs(i)
end do

!write(iw,"(A)") "X, Y = np.meshgrid(xs,ys)"

write(iw,"(A)") "fig = plt.figure()"
write(iw,"(A)") "ax  = plt.axes(projection='3d')"


write(iw,"(A)") "ax.scatter(xs,ys,zs)"

write(iw,"(A)")     'ax.grid()'
write(iw,"(A,A,A)") 'fig.savefig("',filename,'")'

write(iw,"(A,A,A)") 'fig.canvas.set_window_title("',filename,'")'

write(iw,"(A)") 'plt.show()'

close(iw)

 call system(command)
! call system("rm -f plotscript.py")

end subroutine

subroutine randperm(n,iperm)
implicit double precision (a-h,o-z)
integer     :: iperm(n)

do i=1,n
iperm(i) = i
end do

do i=1,n-1
call random_number(dd)
j = i+1 + (n-i-1)*dd
ival     = iperm(j)
iperm(j) = iperm(i)
iperm(i) = ival
end do

end subroutine






        subroutine quicksort2(n,vals,idxs)
        implicit double precision (a-h,o-z)
        dimension istack(2,20000)
        dimension vals(1),idxs(1)
!
!       Sorts the list of double precision numbers in vals and keep track of
!       indices.
!
        if (n .lt. 100) then
        call insort2(n,vals,idxs)
        return
        endif
!
        maxstack = 10000
        k        = 60
!
        m = 1
        istack(1,1) = 1
        istack(2,1) = n
!
 1000 continue
        if (m .eq. 0) goto 1100
        i1 = istack(1,m)
        i2 = istack(2,m)
        m=m-1
!
        l = i2-i1+1
        if (l .le. k) then
        call insort2(l,vals(i1),idxs(i1))
        goto 1000
        endif
!
!       Otherwise perform quicksort.
!
        call quicksort21(vals,idxs,i1,i2,i3)
!
!       This should never happen, but just in case ...
!
        if (m+2 .ge. maxstack) then
        print *,"quicksort out of memory"
        stop
        endif
!
!       Make sure the smaller half is processed first to ensure storage
!       requirements are O(log(n))
!             
        n1 = i3-i1+1
        n2 = i2-i3
!
        if (n2 .lt. n1) then
!
        m = m+1
        istack(1,m) = i1
        istack(2,m) = i3
!
        m = m+1
        istack(1,m) = i3+1
        istack(2,m) = i2
!
        else
!
        m = m+1
        istack(1,m) = i3+1
        istack(2,m) = i2
!
        m = m+1
        istack(1,m) = i1
        istack(2,m) = i3
!
        endif
!
        goto 1000
 1100 continue
        end


        subroutine quicksort21(vals,idxs,i1,i2,i3)
        implicit double precision (a-h,o-z)
        dimension vals(1),idxs(1)
!
!       Randomly choose a pivot index.
!
!        call corrand3(1,r)
        call random_number(r)
        ipiv = i1+floor((i2-i1)*r)
!
        ipiv = i1+(i2-i1)/2
!
        val  = vals(ipiv)
        ival = idxs(ipiv)
!
!       Swap the pivot element and the last element.
!
        vals(ipiv) = vals(i2)
        vals(i2)   = val
!
        idxs(ipiv) = idxs(i2)
        idxs(i2)   = ival
!
        i3 = i1
!
        do 1000 i=i1,i2-1
!
        if( vals(i) .lt. val) then
        d  = vals(i)
        id = idxs(i)
!
        vals(i)  = vals(i3)
        vals(i3) = d       
!
        idxs(i)  = idxs(i3)
        idxs(i3) = id
        i3=i3+1
        endif
 1000 continue
!
        dd = vals(i3)
        vals(i3) = vals(i2)
        vals(i2) = dd
!
        idd = idxs(i3)
        idxs(i3) = idxs(i2)
        idxs(i2) = idd
!
        end


        subroutine insort2(k,a,ib)
        implicit double precision (a-h,o-z)
        dimension a(1),ib(1)
        if (k .le. 1) return
        do 1000 i=2,k
        val=a(i)
        ival=ib(i)
        j=i-1
        do 1100 while (j .ge. 1 .AND. a(j) .gt. val) 
        a(j+1)=a(j)
        ib(j+1)=ib(j)
        j=j-1
 1100 continue
        a(j+1)=val
        ib(j+1)=ival
 1000 continue
        end




subroutine write_table_double(iw,xx)
implicit double precision (a-h,o-z)

if (xx .eq. 0) then
i1 = 0
i2 = 0
i3 = 0
nexp = 0

write (iw,'(I1,".",I1,I1,"\e{+",I2.2,"}")',advance="no") i1,i2,i3,nexp

return
endif

nexp = 0
x    = xx

if (x .gt. 1) then

do while (x .ge. 10.0d0) 
x    = x / 10.0d0
nexp = nexp + 1
end do


i1 = floor(x)
i2 = floor(x*10-10*i1)
i3 = floor(x*100-100*i1-10*i2)

write (iw,'(I1,".",I1,I1,"\e{+",I2.2,"} ")',advance="no") i1,i2,i3,nexp

else

do while (x .le. 1.0d0) 
x    = x * 10.0d0
nexp = nexp + 1
end do


i1 = floor(x)
i2 = floor(x*10-10*i1)
i3 = floor(x*100-100*i1-10*i2)

write (iw,'(I1,".",I1,I1,"\e{-",I2.2,"}")',advance="no") i1,i2,i3,nexp

endif

end subroutine





subroutine write_table_integer_range(iw,nu1,nu2)
implicit double precision (a-h,o-z)

call write_table_integer_sep(iw,nu1)
write(iw,"(A)",advance="no") " - "
call write_table_integer_sep(iw,nu2)

end subroutine



subroutine write_table_integer_sep(iw,n)
implicit double precision (a-h,o-z)

integer :: idigits(1000)

m      = n
ncount = 0
isep   = 0

if (n .eq. 0) then
write(iw,"(A)",advance="no") "0"
return
endif

do while (m .ge. 1)

if (isep .eq. 3) then
ncount          = ncount+1
idigits(ncount) = -1
isep            = 0
endif

ncount          = ncount+1
idigits(ncount) = mod(m,10)
m               = m /10
isep            = isep+1

end do

do i=1,ncount
id = idigits(ncount-i+1)
if (id .ge. 0) then
write(iw,"(I1)",advance="no") id
else
write(iw,"(A)",advance="no") "\sep,"
endif

end do

end subroutine


subroutine write_table_integer_power2(iw,n)
implicit double precision (a-h,o-z)


dd = log(n+0.0d0)/log(2.0d0)
if(dd .ge. 10) then
write(iw,"('$2^{',I2,'}$')",advance="no") int(dd)
else
write(iw,"('$2^{',I1,'}$')",advance="no") int(dd)
endif

end subroutine


subroutine write_table_nextline(iw)
implicit double precision (a-h,o-z)
write (iw,*) " \\"
end subroutine


subroutine write_table_next(iw)
implicit double precision (a-h,o-z)
write (iw,"(' & ')",advance="no") 
end subroutine






subroutine pyplot_begin(iw,istatus)
implicit double precision (a-h,o-z)

integer           :: iw,istatus

!
!  Begin the process of constructing a python script for producing a plot
!  of a function or collection of functions.
!  
!  Input parameters:
!    iw - the fortran unit number of the file to write output to; this must
!      be <= 255
!
!  Output parameters:
!    istatus - a status word used to keep track of certain parameters
!

write(iw,"(A)") "import numpy as np"
write(iw,"(A)") "import matplotlib.pyplot as plt"
write(iw,"(A)") "fig, ax = plt.subplots()"


write(iw,"(A)") "from matplotlib import rc"
write(iw,"(A)") "rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':16})"
write(iw,"(A)") "rc('text', usetex=True)"


istatus = iw


end subroutine



! subroutine pyplot_xlogscale(istatus)
! implicit double precision (a-h,o-z)

! integer               :: istatus
! integer, parameter    :: mask1 = Z'000000FF'
! integer, parameter    :: mask2 = Z'0000FF00'
! integer, parameter    :: mask3 = Z'00010000'
! integer, parameter    :: mask4 = Z'00020000'
! integer, parameter    :: mask5 = Z'00040000'


! iw      = iand(istatus,mask1)
! idx     = ishft(iand(istatus,mask2),-8)
! ilabel  = ishft(iand(istatus,mask3),-16)
! ilogx   = ishft(iand(istatus,mask4),-17)
! ilogy   = ishft(iand(istatus,mask5),-18)

! ilogx = 1

! write (iw,"(A)") "ax.set_xscale('log')"
! istatus = iw + ishft(idx,8) + ishft(ilabel,16) + ishft(ilogx,17) + ishft(ilogy,18)

! end subroutine


! subroutine pyplot_ylogscale(istatus)
! implicit double precision (a-h,o-z)

! integer               :: istatus
! integer, parameter    :: mask1 = Z'000000FF'
! integer, parameter    :: mask2 = Z'0000FF00'
! integer, parameter    :: mask3 = Z'00010000'
! integer, parameter    :: mask4 = Z'00020000'
! integer, parameter    :: mask5 = Z'00040000'


! iw      = iand(istatus,mask1)
! idx     = ishft(iand(istatus,mask2),-8)
! ilabel  = ishft(iand(istatus,mask3),-16)
! ilogx   = ishft(iand(istatus,mask4),-17)
! ilogy   = ishft(iand(istatus,mask5),-18)

! ilogy = 1

! write (iw,"(A)") "ax.set_yscale('log')"
! istatus = iw + ishft(idx,8) + ishft(ilabel,16) + ishft(ilogx,17) + ishft(ilogy,18)

! end subroutine


! subroutine pyplot_xlabel(istatus,xlabel)
! implicit double precision (a-h,o-z)

! !
! !
! !

! integer               :: istatus
! character(len=*)      :: xlabel

! integer, parameter    :: mask1 = Z'000000FF'
! integer, parameter    :: mask2 = Z'0000FF00'
! integer, parameter    :: mask3 = Z'00010000'
! integer, parameter    :: mask4 = Z'00020000'
! integer, parameter    :: mask5 = Z'00040000'


! iw      = iand(istatus,mask1)
! idx     = ishft(iand(istatus,mask2),-8)
! ilabel  = ishft(iand(istatus,mask3),-16)
! ilogx   = ishft(iand(istatus,mask4),-17)
! ilogy   = ishft(iand(istatus,mask5),-18)

! write (iw,"(A,A,A)") 'ax.set(xlabel=r"',xlabel,'")'

! istatus = iw + ishft(idx,8) + ishft(ilabel,16) + ishft(ilogx,17) + ishft(ilogy,18)

! end subroutine


! subroutine pyplot_ylabel(istatus,ylabel)
! implicit double precision (a-h,o-z)

! !
! !
! !

! integer               :: istatus
! character(len=*)      :: ylabel

! integer, parameter    :: mask1 = Z'000000FF'
! integer, parameter    :: mask2 = Z'0000FF00'
! integer, parameter    :: mask3 = Z'00010000'
! integer, parameter    :: mask4 = Z'00020000'
! integer, parameter    :: mask5 = Z'00040000'


! iw      = iand(istatus,mask1)
! idx     = ishft(iand(istatus,mask2),-8)
! ilabel  = ishft(iand(istatus,mask3),-16)
! ilogx   = ishft(iand(istatus,mask4),-17)
! ilogy   = ishft(iand(istatus,mask5),-18)

! write (iw,"(A,A,A)") 'ax.set(ylabel=r"',ylabel,'")'

! istatus = iw + ishft(idx,8) + ishft(ilabel,16) + ishft(ilogx,17) + ishft(ilogy,18)

! end subroutine


! subroutine pyplot_add_function(istatus,istyle,label,n,xs,ys)
! implicit double precision (a-h,o-z)

! integer                  :: istyle, idx
! character(len=*)         :: label
! double precision         :: xs(n),ys(n)
! character(len=4)         :: xname,yname
! character(:),allocatable :: style

! integer               :: istatus
! integer, parameter    :: mask1 = Z'000000FF'
! integer, parameter    :: mask2 = Z'0000FF00'
! integer, parameter    :: mask3 = Z'00010000'
! integer, parameter    :: mask4 = Z'00020000'
! integer, parameter    :: mask5 = Z'00040000'


! iw      = iand(istatus,mask1)
! idx     = ishft(iand(istatus,mask2),-8)
! ilabel  = ishft(iand(istatus,mask3),-16)
! ilogx   = ishft(iand(istatus,mask4),-17)
! ilogy   = ishft(iand(istatus,mask5),-18)


! idx     = idx+1

! write (xname,"(A,I2.2)") "xs",idx
! write (yname,"(A,I2.2)") "ys",idx


! write(iw,"(A,A,I5,A)") xname," = np.zeros(",n,")"
! write(iw,"(A,A,I5,A)") yname," = np.zeros(",n,")"

! do i=1,n
! write(iw,"(A,A,I5,A,E24.16)") xname,"[",i-1,"] = ",xs(i)
! write(iw,"(A,A,I5,A,E24.16)") yname,"[",i-1,"] = ",ys(i)
! end do


! if (istyle .eq. 1) then
! allocate(character(2) :: style)
! style='b-'
! elseif (istyle .eq. 2) then
! allocate(character(2) :: style)
! style = "b:"
! elseif (istyle .eq. 3) then
! allocate(character(3) :: style)
! style = "b--"
! elseif (istyle .eq. 4) then
! allocate(character(3) :: style)
! style = "b-."
! elseif (istyle .eq. 5) then
! allocate(character(2) :: style)
! style='r-'
! elseif (istyle .eq. 6) then
! allocate(character(2) :: style)
! style = "r:"
! elseif (istyle .eq. 7) then
! allocate(character(3) :: style)
! style = "r--"
! elseif (istyle .eq. 8) then
! allocate(character(3) :: style)
! style = "r-."
! elseif (istyle .eq. 9) then
! allocate(character(2) :: style)
! style='k-'
! elseif (istyle .eq. 10) then
! allocate(character(2) :: style)
! style = "k:"
! elseif (istyle .eq. 11) then
! allocate(character(3) :: style)
! style = "k--"
! elseif (istyle .eq. 12) then
! allocate(character(3) :: style)
! style = "k-."
! endif


! if (len(label) .gt. 0) then
! write(iw,"(A,A,A,A,A,A,A,A,A,A)") "ax.plot(",xname,",",yname,",'",&
!   style,"'",',label = r"',label,'",linewidth=2)'
! ilabel = 1
! else
! write(iw,"(A,A,A,A,A,A,A)") "ax.plot(linewidth=.5,",xname,",",yname,",'",style,&
!   ",linewidth=2')"
! endif


! istatus = iw + ishft(idx,8) + ishft(ilabel,16) + ishft(ilogx,17) + ishft(ilogy,18)


! end subroutine


! subroutine pyplot_end(istatus,out_name)
! implicit double precision (a-h,o-z)

! character(len=*)    :: out_name

! !
! !  Input parameters:
! !


! integer               :: istatus
! integer, parameter    :: mask1 = Z'000000FF'
! integer, parameter    :: mask2 = Z'0000FF00'
! integer, parameter    :: mask3 = Z'00010000'
! integer, parameter    :: mask4 = Z'00020000'
! integer, parameter    :: mask5 = Z'00040000'


! iw      = iand(istatus,mask1)
! idx     = ishft(iand(istatus,mask2),-8)+1
! ilabel  = ishft(iand(istatus,mask3),-16)
! ilogx   = ishft(iand(istatus,mask4),-17)
! ilogy   = ishft(iand(istatus,mask5),-18)


! if (ilabel .eq. 1) then
! write(iw,"(A)")     "legend=ax.legend(loc='upper left', fontsize='small')"
! endif

! write(iw,"(A)") 'ax.grid()'
! write(iw,"(A,A,A)") 'fig.savefig("',out_name,'")'
! write(iw,"(A)") 'plt.show()'

! end subroutine



subroutine equispaced_intervals(nints,a,b,ab)
implicit double precision (a-h,o-z)

integer, intent(in)            :: nints
double precision, intent(in)   :: a,b
double precision, allocatable, intent(out) :: ab(:,:)

allocate(ab(2,nints))
do int=1,nints
ab(1,int) = a + (b-a) * (int-1.0d0)/(nints)
ab(2,int) = a + (b-a) * (int+0.0d0)/(nints)
end do

end subroutine


subroutine bisected_intervals(nints,a,b,ab)
implicit double precision (a-h,o-z)

integer, intent(in)            :: nints
double precision, intent(in)   :: a,b
double precision, allocatable, intent(out) :: ab(:,:)

allocate(ab(2,nints))

do int=1,nints
ab(1,int) = 2.0d0**(-nints+int-1) 
ab(2,int) = 2.0d0**(-nints+int)
end do

ab = a + (b-a)*ab
end subroutine




end module
