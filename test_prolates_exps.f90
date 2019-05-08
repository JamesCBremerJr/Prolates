program test_prolates_exps
use utils
use chebyshev
use tensor
use prolates
use prolates_window
use prolates_exps

implicit double precision (a-h,o-z)

type(chebexps_data)                 :: chebdata

double precision, allocatable       :: ab0(:,:),chivals0(:,:),apvals0(:,:),appvals0(:,:)
double precision, allocatable       :: apppvals0(:,:)

double precision, allocatable       :: ab(:,:),cd(:,:),chicoefs(:,:,:),eta1coefs(:,:)
double precision, allocatable       :: eta2coefs(:,:),apcoefs(:,:,:),appcoefs(:,:,:)
double precision, allocatable       :: apppcoefs(:,:,:),w(:),coefs(:)

double complex                      :: clambda
double precision, allocatable       :: coefsps(:),coefsqs(:)


ic1    = 0
ic2    = 1

eps0   = epsilon(0.0d0)
eps    = 1.0d-12

if (eps0 .lt. 1.0d-17) then
eps    = 1.0d-14
endif

pi     = acos(-1.0d0)
k      = 30
call chebexps(k,chebdata)

!
!  Test expansion built for fixed c 
!

call prina("---------------[ testing prolates_exp_one ]------------------------")
call prina("")

irange = 1
isubst = 0

c      = 5000.0d0
call prolates_exps_one(eps,irange,chebdata,c,isubst,nints0,ab0,&
  chivals0,apvals0,appvals0,apppvals0,eta1,eta2)
eta = ceiling((3*eta1+eta2)/4)
u   = (eta-eta1)/(eta2-eta1)
n   = eta
call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)
call prolates_window_avals(c,chi,isubst,aval,apval,appval,apppval)

call chebpw_eval(nints0,ab0,k,chebdata%xs,chivals0,u,chi0)
call chebpw_eval(nints0,ab0,k,chebdata%xs,apvals0,u,apval0)
call chebpw_eval(nints0,ab0,k,chebdata%xs,appvals0,u,appval0)
call chebpw_eval(nints0,ab0,k,chebdata%xs,apppvals0,u,apppval0)
errrel = abs(chi-chi0)/abs(chi)
call prind("relative error in chi                   = ",errrel)
errrel = abs(apval-apval0)/abs(apval)
call prind("relative error in alpha'(0)             = ",errrel)
call prina("")




c      = 50000.0d0
call prolates_exps_one(eps,irange,chebdata,c,isubst,nints0,ab0,&
  chivals0,apvals0,appvals0,apppvals0,eta1,eta2)
eta = ceiling((3*eta1+eta2)/4)
u   = (eta-eta1)/(eta2-eta1)
n   = eta
call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)
call prolates_window_avals(c,chi,isubst,aval,apval,appval,apppval)
call chebpw_eval(nints0,ab0,k,chebdata%xs,chivals0,u,chi0)
call chebpw_eval(nints0,ab0,k,chebdata%xs,apvals0,u,apval0)
errrel = abs(chi-chi0)/abs(chi)
call prind("relative error in chi                   = ",errrel)
errrel = abs(apval-apval0)/abs(apval)
call prind("relative error in alpha'(0)             = ",errrel)
call prina("")

c      = 500000.0d0
call prolates_exps_one(eps,irange,chebdata,c,isubst,nints0,ab0,&
  chivals0,apvals0,appvals0,apppvals0,eta1,eta2)
eta = ceiling((3*eta1+eta2)/4)
u   = (eta-eta1)/(eta2-eta1)
n   = eta
call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)
call prolates_window_avals(c,chi,isubst,aval,apval,appval,apppval)
call chebpw_eval(nints0,ab0,k,chebdata%xs,chivals0,u,chi0)
call chebpw_eval(nints0,ab0,k,chebdata%xs,apvals0,u,apval0)
errrel = abs(chi-chi0)/abs(chi)
call prind("relative error in chi                   = ",errrel)
errrel = abs(apval-apval0)/abs(apval)
call prind("relative error in alpha'(0)             = ",errrel)
call prina("")

c      = 1000000.0d0
call prolates_exps_one(eps,irange,chebdata,c,isubst,nints0,ab0,&
  chivals0,apvals0,appvals0,apppvals0,eta1,eta2)
eta = ceiling((3*eta1+eta2)/4)
u   = (eta-eta1)/(eta2-eta1)
n   = eta
call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)
call prolates_window_avals(c,chi,isubst,aval,apval,appval,apppval)
call chebpw_eval(nints0,ab0,k,chebdata%xs,chivals0,u,chi0)
call chebpw_eval(nints0,ab0,k,chebdata%xs,apvals0,u,apval0)
errrel = abs(chi-chi0)/abs(chi)
call prind("relative error in chi                   = ",errrel)
errrel = abs(apval-apval0)/abs(apval)
call prind("relative error in alpha'(0)             = ",errrel)
call prina("")

irange = 2
isubst = 0

c      = 5000.0d0
call prolates_exps_one(eps,irange,chebdata,c,isubst,nints0,ab0,&
  chivals0,apvals0,appvals0,apppvals0,eta1,eta2)
eta = ceiling((3*eta1+eta2)/4)
u   = (eta-eta1)/(eta2-eta1)
n   = eta
call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)
call prind("chi = ",chi)

call prolates_window_avals(c,chi,isubst,aval,apval,appval,apppval)
call chebpw_eval(nints0,ab0,k,chebdata%xs,chivals0,u,chi0)
call chebpw_eval(nints0,ab0,k,chebdata%xs,apvals0,u,apval0)
errrel = abs(chi-chi0)/abs(chi)
call prind("relative error in chi                   = ",errrel)
errrel = abs(apval-apval0)/abs(apval)
call prind("relative error in alpha'(0)             = ",errrel)
call prina("")


c      = 50000.0d0
call prolates_exps_one(eps,irange,chebdata,c,isubst,nints0,ab0,&
  chivals0,apvals0,appvals0,apppvals0,eta1,eta2)
eta = ceiling((3*eta1+eta2)/4)
u   = (eta-eta1)/(eta2-eta1)
n   = eta
call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)
call prolates_window_avals(c,chi,isubst,aval,apval,appval,apppval)
call chebpw_eval(nints0,ab0,k,chebdata%xs,chivals0,u,chi0)
call chebpw_eval(nints0,ab0,k,chebdata%xs,apvals0,u,apval0)
errrel = abs(chi-chi0)/abs(chi)
call prind("relative error in chi                   = ",errrel)
errrel = abs(apval-apval0)/abs(apval)
call prind("relative error in alpha'(0)             = ",errrel)
call prina("")

c      = 500000.0d0
call prolates_exps_one(eps,irange,chebdata,c,isubst,nints0,ab0,&
  chivals0,apvals0,appvals0,apppvals0,eta1,eta2)
eta = ceiling((3*eta1+eta2)/4)
u   = (eta-eta1)/(eta2-eta1)
n   = eta
call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)
call prolates_window_avals(c,chi,isubst,aval,apval,appval,apppval)
call chebpw_eval(nints0,ab0,k,chebdata%xs,chivals0,u,chi0)
call chebpw_eval(nints0,ab0,k,chebdata%xs,apvals0,u,apval0)
errrel = abs(chi-chi0)/abs(chi)
call prind("relative error in chi                   = ",errrel)
errrel = abs(apval-apval0)/abs(apval)
call prind("relative error in alpha'(0)             = ",errrel)
call prina("")

c      = 1000000.0d0
call prolates_exps_one(eps,irange,chebdata,c,isubst,nints0,ab0,&
  chivals0,apvals0,appvals0,apppvals0,eta1,eta2)
eta = ceiling((3*eta1+eta2)/4)
u   = (eta-eta1)/(eta2-eta1)
n   = eta
call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)
call prolates_window_avals(c,chi,isubst,aval,apval,appval,apppval)
call chebpw_eval(nints0,ab0,k,chebdata%xs,chivals0,u,chi0)
call chebpw_eval(nints0,ab0,k,chebdata%xs,apvals0,u,apval0)
errrel = abs(chi-chi0)/abs(chi)
call prind("relative error in chi                   = ",errrel)
errrel = abs(apval-apval0)/abs(apval)
call prind("relative error in alpha'(0)             = ",errrel)
call prina("")



!
!  Test expansion of c and eta
!

call prina("---------------[ testing prolates_exps_construct ]----------------------------")



irange = 1
isubst = 0
call prolates_exps_construct(eps,irange,ic1,ic2,chebdata,isubst,nintsab,ab,nintscd,cd, &
  chicoefs,apcoefs,appcoefs,apppcoefs,eta1coefs,eta2coefs)

dsize = 8.0d0/(1024.0d0*1024.0d0)
dsize = dsize *k * (k+1)/2
dsize = dsize * nintscd * nintsab

call prini("nintsab = ",nintsab)
call prin2("ab = ",ab)
call prini("nintscd = ",nintscd)
call prin2("cd = ",cd)

call prin2("expansion size (in MB) = ",dsize)


do intcd=1,nintscd

c1     = cd(1,intcd)
c2     = cd(2,intcd)
c      = (3*c1+c2)/4

call chebpw_eval2(nintscd,cd,k,eta1coefs,c,eta1)
call chebpw_eval2(nintscd,cd,k,eta2coefs,c,eta2)

n     = (eta1+eta2)/2
u     = (n-eta1)/(eta2-eta1)

call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)
call prolates_window_avals(c,chi,isubst,aval,apval,appval,apppval)

call tensorpw_eval(k,nintsab,ab,nintscd,cd,chicoefs,u,c,chi0)
call tensorpw_eval(k,nintsab,ab,nintscd,cd,apcoefs,u,c,apval0)
call tensorpw_eval(k,nintsab,ab,nintscd,cd,appcoefs,u,c,appval0)
call tensorpw_eval(k,nintsab,ab,nintscd,cd,apppcoefs,u,c,apppval0)

print *,apval,apval0
print *,appval,appval0
print *,apppval,apppval0

stop

call prind("chi as computed via 2d expansion        = ",chi0)
call prind("chi as computed via Rokhlin's algorithm = ",chi)
errrel = abs(chi-chi0)/abs(chi)
call prind("relative error in chi                   = ",errrel)


errrel = abs(apval-apval0)/abs(apval)
call prind("relative error in alpha'(0)             = ",errrel)
end do


stop


irange = 2
call prolates_exps_construct(eps,irange,ic1,ic2,chebdata,isubst,nintsab,ab,nintscd,cd, &
  chicoefs,apcoefs,appcoefs,apppcoefs,eta1coefs,eta2coefs)

dsize = 8.0d0/(1024.0d0*1024.0d0)
dsize = dsize *k * (k+1)/2
dsize = dsize * nintscd * nintsab

call prini("nintsab = ",nintsab)
call prin2("ab = ",ab)
call prini("nintscd = ",nintscd)
call prin2("cd = ",cd)

call prin2("expansion size (in MB) = ",dsize)


do intcd=1,nintscd

c1     = cd(1,intcd)
c2     = cd(2,intcd)
c      = (3*c1+c2)/4

call chebpw_eval2(nintscd,cd,k,eta1coefs,c,eta1)
call chebpw_eval2(nintscd,cd,k,eta2coefs,c,eta2)

n     = (eta1+eta2)/2
u     = (n-eta1)/(eta2-eta1)

call prolates_integer(c,n,chi,clambda,ncoefs,coefsps,coefsqs)
call prolates_window_avals(c,chi,isubst,aval,apval,appval,apppval)

call tensorpw_eval(k,nintsab,ab,nintscd,cd,chicoefs,u,c,chi0)
call tensorpw_eval(k,nintsab,ab,nintscd,cd,apcoefs,u,c,apval0)

call prind("chi as computed via 2d expansion        = ",chi0)
call prind("chi as computed via Rokhlin's algorithm = ",chi)
errrel = abs(chi-chi0)/abs(chi)
call prind("relative error in chi                   = ",errrel)


errrel = abs(apval-apval0)/abs(apval)
call prind("relative error in alpha'(0)             = ",errrel)
end do


end program
