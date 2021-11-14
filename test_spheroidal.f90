program test_spheroidal
use spheroidal

implicit double precision (a-h,o-z)
double precision, allocatable :: coefsps(:)
double precision, allocatable :: vals(:)

pi    = acos(-1.0d0)
gamma =   10
m     =   10
n     =   21
x     =   0.9d0

! reference values
chi0   = 501.101795657320319701294212789709019d0
valps0 = 0.551920158973020332856865217900369385d0

call spheroidal_integer(gamma,n,m,chi,n1,n2,coefsps)

call prini("m = ",m)
call prini("n = ",n)
call prind("chi = ",chi)
call prini("n1 = ",2*n1+n)
call prini("n2 = ",2*n2+n)
call prin2("coefsps = ",coefsps)

call spheroidal_integer_eval(n,m,n1,n2,coefsps,x,valps)
err1 = abs(chi-chi0)/abs(chi0)
err2 = abs(valps-valps0)/abs(valps0)

call prin2("relative error in chi = ",err1)
call prin2("relative error in valps = ",err2)


end program
