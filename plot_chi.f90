!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This program produces various plots of the Sturm-Liouville eigenvalue chi as
!  a function of characteristic exponent nu and as a function of the value of 
!  a multiple of the nonoscillatory phase at zero, which we will call eta.
!
!  The plots are written to the PDF files chiplot?.pdf.  The plots are
!  generated via python scripts which call matplotlib and which are produced by this 
!  Fortran code.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program plot_chi

use utils
use chebyshev
use prolates
use prolates_window
use prolates_exps

implicit double precision (a-h,o-z)

double precision, allocatable       :: coefs(:,:),chis(:)
double precision, allocatable       :: xs(:,:),ys(:,:),xs2(:),ys2(:)
double precision, allocatable       :: xs3(:),ys3(:),xs4(:,:),ys4(:,:)


type(chebexps_data)                 :: chebdata
double precision, allocatable       :: ab0(:,:),chivals0(:,:),apvals0(:,:),appvals0(:,:)
double precision, allocatable       :: apppvals0(:,:)

c    =  500.0d0

!
!  Evaluate chi as a function of the characteristic exponent
! 

ndnus = 20
n1    = 100
n2    = 150

n3    = 150
n4    = 160

nn    = ndnus*n

allocate(xs(ndnus,n1:n4),ys(ndnus,n1:n4))

a    = -0.5d0
b    =  0.5d0

do i=1,ndnus

dnu0 = a + (b-a)/(ndnus + 1.0d0) * (i)

call elapsed(t1)
call prolates_noninteger_all(c,dnu0,n4,m,chis,coefs)
call elapsed(t2)
print *,"calculting chi_nu(c) ...   "," i = ",i,"(",t2-t1," seconds)"


do n=n1,n4
xs(i,n) = dnu0+n
ys(i,n) = chis(n)
end do

end do


!
!  Evaluate chi as a function of eta = -2/pi (alpha+1)
!


eps    = 1.0d-12
irange = 1
isubst = 0
k      = 30


call chebexps(k,chebdata)

call prolates_exps_one(eps,irange,chebdata,c,isubst,nints0,ab0,&
  chivals0,apvals0,appvals0,apppvals0,eta1,eta2)


a = ab0(1,1)
b = ab0(2,nints0)

nn2      = 1000
allocate(xs2(nn2),ys2(nn2))

do n=1,nn2

u  = (n-1.0d0)/(nn2-1.0d0)
dn = eta1 + (eta2-eta1) *u

call chebpw_eval(nints0,ab0,k,chebdata%xs,chivals0,u,chi)
print *,dn,chi

xs2(n)  = dn
ys2(n)  = chi

end do




!
!  Construct the first plot -- chi as a function of characeristic exponent
!

iw = 20
open(iw,FILE="chiplot1.py")
write(iw,"(A)") "import numpy as np"
write(iw,"(A)") "import matplotlib as mpl"
write(iw,"(A)") "mpl.use('Agg')"
write(iw,"(A)") "import matplotlib.pyplot as plt"

write(iw,"(A)") "fig, ax = plt.subplots()"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iw,"(A,I5,A)") "xs = np.zeros(",ndnus,")"
write(iw,"(A,I5,A)") "ys = np.zeros(",ndnus,")"

do n=n1,n2
do i=1,ndnus
write(iw,"(A,I5,A,E24.16)") "xs[",i-1,"] = ",xs(i,n)
write(iw,"(A,I5,A,E24.16)") "ys[",i-1,"] = ",ys(i,n)
end do
write(iw,"(A)") "ax.plot(xs,ys,'b')"
end do


write(iw,"(A)") 'ax.grid()'
write(iw,"(A,A,A)") 'fig.savefig("chiplot1.pdf")'

close(iw)
call system("python chiplot1.py")
!call system("rm -f chiplot1.py")

!
!  Construct the second plot --- chi as a function of eta
!

call plot_function("chiplot2.pdf","","",nn2,xs2,ys2)



!
!  Now build a combined plot with both functions
!



! a   = n3
! b   = n4
! nn  = n4-n3+1

! allocate(xs3(nn),ys3(nn))

! do i=1,nn

! dd = (i-1.0d0)/(nn-1.0d0)
! n  = n3 + (n4-n3)*dd
! u  = (n-eta1)/(eta2-eta1)

! call chebpw_eval(nints0,ab0,k,chebdata%xs,chivals0,u,chi)
! xs3(i) = n
! ys3(i) = chi
! end do




! !
! !  Construct the combined plot
! !

! iw = 20
! open(iw,FILE="chiplot3.py")
! write(iw,"(A)") "import numpy as np"
! write(iw,"(A)") "import matplotlib as mpl"
! write(iw,"(A)") "mpl.use('Agg')"
! write(iw,"(A)") "import matplotlib.pyplot as plt"

! write(iw,"(A)") "fig, ax = plt.subplots()"

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write(iw,"(A,I5,A)") "xs = np.zeros(",ndnus,")"
! write(iw,"(A,I5,A)") "ys = np.zeros(",ndnus,")"

! do n=n3,n4
! do i=1,ndnus
! write(iw,"(A,I5,A,E24.16)") "xs[",i-1,"] = ",xs(i,n)
! write(iw,"(A,I5,A,E24.16)") "ys[",i-1,"] = ",ys(i,n)
! end do
! write(iw,"(A)") "ax.plot(xs,ys,'b')"
! end do

! write(iw,"(A,I5,A)") "xs3 = np.zeros(",nn,")"
! write(iw,"(A,I5,A)") "ys3 = np.zeros(",nn,")"
! do i=1,nn
! write(iw,"(A,I5,A,E24.16)") "xs3[",i-1,"] = ",xs3(i)
! write(iw,"(A,I5,A,E24.16)") "ys3[",i-1,"] = ",ys3(i)
! end do
! write(iw,"(A)") "ax.plot(xs3,ys3,'k')"


! write(iw,"(A)") 'ax.grid()'
! write(iw,"(A,A,A)") 'fig.savefig("chiplot3.pdf")'

! close(iw)
! call system("python chiplot3.py")

end program
