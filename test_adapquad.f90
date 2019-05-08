module test_adapquad_funs

contains

subroutine testfun2(n,xs,vals,userptr)
use    iso_c_binding
implicit double precision (a-h,o-z)
double precision         :: xs(n)
type(c_ptr)              :: userptr
double precision         :: vals(n)
vals = log(xs**2)
end subroutine



end module


program test_adapquad
use utils
use adapquad
use iso_c_binding
use test_adapquad_funs

implicit double precision (a-h,o-z)
type(c_ptr)              :: userptr
double precision         :: val,val0

eps  = 1.0d-13
x1   = 0.0d0
x2   = 1.0d0
call adapint(ier,eps,x1,x2,testfun2,userptr,val,ntotal)

val0 = -2.0d0

call prind("integral value = ",val)
call prind("relative error in value = ",(val-val0)/val0)
call prini("nodes in quadrature = ",ntotal)

end program
