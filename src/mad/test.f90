program test
use esvar
use estatic, only: setmadvar, initmadelung

implicit none

 ldip = 2;
 call setmadvar()
 call initmadelung()


 write(*,*)'test program ended...'
end program test
