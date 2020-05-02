
! Ahsan 02-05-2020
! This module contains ELK-6.2.8 files that calculate gaunt coefficients
! gaunt coeff will be used to calculate electron-electron interaction tensors Vee
! for each TM species using Slater-integrals F(0), F(2),F(4) [for d-orbitals].
! everything borrowed from ELK assumes complex spherical harmonics,
! so we need to change the basis from our real spherical harmnics to these
! complex in our atomic density matrices, calc Vee in complex Ylm basis, 
! and then transform the basis of Vee back to our real Ylm. 
! (our real: our sign convention, standard, not the convention by Siesta/Fernandos used for spin-orbit calculations)
module modgaunt
implicit none

!real(8) :: gaunt, wigner3j,factnm,factr

!public :: mkdgaunt ! our driver routine
!private:: gaunt, wigner3j,factnm,factr ! all from ELK code

contains

!======================================================================
! our driver routine
! calculate the matrix g(k,m4,m3,m2,m1) that can be used to calculate Vee
! for any TM species
subroutine mkdgaunt(l,g)
implicit none
integer, intent(in) :: l
double precision, dimension(0:l,-l:l,-l:l,-l:l,-l:l), intent(out) :: g
! local
integer :: m1,m2,m3,m4,k,q
double precision :: sum1, sum2, t1
! external
!real(8) gaunt
!external gaunt

do m1=-l,l
  do m2=-l,l
    do m3=-l,l
      do m4=-l,l
        sum1=0.d0
        do k=0,2*l,2
          sum2=0.d0
          do q=-k,k
            t1=gaunt(l,k,l,m1,q,m2)*gaunt(l,k,l,m3,-q,m4)
            if (mod(q,2).eq.0) then
              sum2=sum2+t1
            else
              sum2=sum2-t1
            end if
          end do
          g(k/2,m4,m3,m2,m1) = sum2/dble(2*k+1)          
          !sum1=sum1+f(k)*sum2/dble(2*k+1)
        end do
        !vee(m1,m3,m2,m4)=fourpi*sum1
      end do
    end do
  end do
end do

return
end subroutine mkdgaunt

!======================================================================
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: gaunt
! !INTERFACE:
real(8) function gaunt(l1,l2,l3,m1,m2,m3)
! !INPUT/OUTPUT PARAMETERS:
!   l1, l2, l3 : angular momentum quantum numbers (in,integer)
!   m1, m2, m3 : magnetic quantum numbers (in,integer)
! !DESCRIPTION:
!   Returns the Gaunt coefficient given by
!   $$  \langle Y^{l_1}_{m_1}|Y^{l_2}_{m_2}|Y^{l_3}_{m_3} \rangle
!    = (-1)^{m_1}\left[\frac{(2l_1+1)(2l_2+1)(2l_3+1)}{4\pi} \right]
!    ^{\frac{1}{2}}
!    \begin{pmatrix} l_1 & l_2 & l_3 \\  0   & 0   & 0   \end{pmatrix}
!    \begin{pmatrix} l_1 & l_2 & l_3 \\ -m_1 & m_2 & m_3 \end{pmatrix}. $$
!   Suitable for $l_i$ less than 50.
!
! !REVISION HISTORY:
!   Created November 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: l1,l2,l3
integer, intent(in) :: m1,m2,m3
! local variables
integer j,j1,j2,j3,jh
real(8) t1
! real constant 1/sqrt(4*pi)
real(8), parameter :: c1=0.28209479177387814347d0
! external functions
!real(8) wigner3j,factnm,factr
!external wigner3j,factnm,factr
if ((l1.lt.0).or.(l2.lt.0).or.(l3.lt.0).or.(abs(m1).gt.l1).or.(abs(m2).gt.l2) &
 .or.(abs(m3).gt.l3)) then
  write(*,*)
  write(*,'("Error(gaunt): non-physical arguments :")')
  write(*,'("l1 = ",I8," l2 = ",I8," l3 = ",I8)') l1,l2,l3
  write(*,'("m1 = ",I8," m2 = ",I8," m3 = ",I8)') m1,m2,m3
  write(*,*)
  stop
end if
if ((l1.gt.50).or.(l2.gt.50).or.(l3.gt.50)) then
  write(*,*)
  write(*,'("Error(gaunt): angular momenta out of range : ",3I8)') l1,l2,l3
  write(*,*)
  stop
end if
if (m1-m2-m3.ne.0) then
  gaunt=0.d0
  return
end if
j1=l2-l1+l3
j2=l1-l2+l3
j3=l1+l2-l3
if ((j1.lt.0).or.(j2.lt.0).or.(j3.lt.0)) then
  gaunt=0.d0
  return
end if
j=l1+l2+l3
if (mod(j,2).ne.0) then
  gaunt=0.d0
  return
end if
jh=j/2
t1=sqrt(dble((2*l1+1)*(2*l2+1)*(2*l3+1))*factr(j1,j+1)*factnm(j2,1) &
 *factnm(j3,1))
t1=t1*factr(jh,jh-l1)/(factnm(jh-l2,1)*factnm(jh-l3,1))
gaunt=t1*c1*wigner3j(l1,l2,l3,-m1,m2,m3)
if (mod(m1+jh,2).ne.0) gaunt=-gaunt
return
end function
!EOC

!======================================================================
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: wigner3j
! !INTERFACE:
real(8) function wigner3j(j1,j2,j3,m1,m2,m3)
! !INPUT/OUTPUT PARAMETERS:
!   j1, j2, j3 : angular momentum quantum numbers (in,integer)
!   m1, m2, m3 : magnetic quantum numbers (in,integer)
! !DESCRIPTION:
!   Returns the Wigner $3j$-symbol. There are many equivalent formulae for
!   the $3j$-symbols, the following provides high accuracy for $j\le 50$
!   \begin{align*}
!    &\begin{pmatrix} j_1 & j_2 & j_3 \\ m_1 & m_2 & m_3 \end{pmatrix}= \\
!    &(-1)^{j1+j2+m3}\sqrt{\frac{(j_1+m_1)!\,(j_2+m_2)!\,(j_3+m_3)!\,
!    (j_3-m_3)!\,(j_1-m_1)!\,(j_2-m_2)!}{(j_2-j_1+j_3)!\,(j_1-j_2+j_3)!\,
!    (j_1+j_2-j_3)!\,(1+j_1+j_2+j_3)!}}\,\sum_k(-1)^k \\
!    &\frac{(j_2-j_1+j_3)!\,(j_1-j_2+j_3)!\,(j_1+j_2-j_3)!}{(j_3-j_1-m_2+k)!\,
!    (j_3-j_2+m_1+k)!\,(j_1+j_2-j_3-k)!\,k!\,(j_1-m_1-k)!\,(j_2+m_2-k)!},
!   \end{align*}
!   where the sum is over all integers $k$ for which the factorials in the
!   summand are non-negative.
!
! !REVISION HISTORY:
!   Created November 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: j1,j2,j3
integer, intent(in) :: m1,m2,m3
! local variables
integer k,k1,k2,l1,l2,l3,n1,n2
real(8) sgn,sum,t1
! external functions
!real(8) factnm,factr
!external factnm,factr
! check input variables
if ((j1.lt.0).or.(j2.lt.0).or.(j3.lt.0).or.(abs(m1).gt.j1).or.(abs(m2).gt.j2) &
 .or.(abs(m3).gt.j3)) then
  write(*,*)
  write(*,'("Error(wigner3j): invalid arguments :")')
  write(*,'("j1 = ",I8," j2 = ",I8," j3 = ",I8)') j1,j2,j3
  write(*,'("m1 = ",I8," m2 = ",I8," m3 = ",I8)') m1,m2,m3
  write(*,*)
  stop
end if
if ((j1.eq.0).and.(j2.eq.0).and.(j3.eq.0)) then
  wigner3j=1.d0
  return
end if
if ((j1.gt.50).or.(j2.gt.50).or.(j3.gt.50)) then
  write(*,*)
  write(*,'("Error(wigner3j): angular momenta out of range : ",3I8)') j1,j2,j3
  write(*,*)
  stop
end if
l1=j2-j1+j3
l2=j1-j2+j3
l3=j1+j2-j3
if ((m1+m2+m3.ne.0).or.(l1.lt.0).or.(l2.lt.0).or.(l3.lt.0)) then
  wigner3j=0.d0
  return
end if
n1=j1-m1
n2=j2+m2
k1=max(0,n1-l2,n2-l1)
k2=min(l3,n1,n2)
if (mod(k1-j1+j2+m3,2).ne.0) then
  sgn=-1.d0
else
  sgn=1.d0
end if
sum=0.d0
do k=k1,k2
  t1=sgn*factr(l1,l1-n2+k)*factr(l2,l2-n1+k)*factr(l3,l3-k)
  sum=sum+t1/(factnm(k,1)*factnm(n1-k,1)*factnm(n2-k,1))
  sgn=-sgn
end do
t1=factr(j1+m1,l1)*factr(j2+m2,l2)*factr(j3+m3,l3)
t1=t1*factr(j3-m3,1+j1+j2+j3)*factnm(j1-m1,1)*factnm(j2-m2,1)
wigner3j=sum*sqrt(t1)
return
end function
!EOC
!======================================================================
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: factr
! !INTERFACE:
real(8) function factr(n,d)
! !INPUT/OUTPUT PARAMETERS:
!   n : numerator (in,integer)
!   d : denominator (in,integer)
! !DESCRIPTION:
!   Returns the ratio $n!/d!$ for $n,d\ge 0$. Performs no under- or overflow
!   checking.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n,d
! local variables
integer i
! external functions
!real(8) factnm
!external factnm
if (d.eq.1) then
  factr=factnm(n,1)
  return
end if
if (n.lt.0) then
  write(*,*)
  write(*,'("Error(factr): n < 0 : ",I8)') n
  write(*,*)
  stop
end if
if (d.lt.0) then
  write(*,*)
  write(*,'("Error(factr): d < 0 : ",I8)') d
  write(*,*)
  stop
end if
if (n.lt.d) then
  factr=dble(n+1)
  do i=n+2,d
    factr=factr*dble(i)
  end do
  factr=1.d0/factr
else if (n.eq.d) then
  factr=1.d0
else
  factr=dble(d+1)
  do i=d+2,n
    factr=factr*dble(i)
  end do
end if
return
end function
!EOC
!======================================================================
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: factnm
! !INTERFACE:
real(8) function factnm(n,m)
! !INPUT/OUTPUT PARAMETERS:
!   n : input (in,integer)
!   m : order of multifactorial (in,integer)
! !DESCRIPTION:
!   Returns the multifactorial
!   $$ n\underbrace{!!\,\cdots\,!}_{m\,{\rm times}}=
!    \prod_{\substack{i\ge 0\\ n-im>0}}(n-im) $$
!   for $n,\,m \ge 0$. $n$ should be less than 150.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n,m
! local variables
integer i,j
real(8) f1(24),f2(38)
data f1 / &
                       1.d0,                        2.d0,  &
                       6.d0,                       24.d0,  &
                     120.d0,                      720.d0,  &
                    5040.d0,                    40320.d0,  &
                  362880.d0,                  3628800.d0,  &
                39916800.d0,                479001600.d0,  &
              6227020800.d0,              87178291200.d0,  &
           1307674368000.d0,           20922789888000.d0,  &
         355687428096000.d0,         6402373705728000.d0,  &
      121645100408832000.d0,      2432902008176640000.d0,  &
    51090942171709440000.d0,   1124000727777607680000.d0,  &
 25852016738884976640000.d0, 620448401733239439360000.d0 /
data f2 / &
                       1.d0,                        2.d0,  &
                       3.d0,                        8.d0,  &
                      15.d0,                       48.d0,  &
                     105.d0,                      384.d0,  &
                     945.d0,                     3840.d0,  &
                   10395.d0,                    46080.d0,  &
                  135135.d0,                   645120.d0,  &
                 2027025.d0,                 10321920.d0,  &
                34459425.d0,                185794560.d0,  &
               654729075.d0,               3715891200.d0,  &
             13749310575.d0,              81749606400.d0,  &
            316234143225.d0,            1961990553600.d0,  &
           7905853580625.d0,           51011754393600.d0,  &
         213458046676875.d0,         1428329123020800.d0,  &
        6190283353629375.d0,        42849873690624000.d0,  &
      191898783962510625.d0,      1371195958099968000.d0,  &
     6332659870762850625.d0,     46620662575398912000.d0,  &
   221643095476699771875.d0,   1678343852714360832000.d0,  &
  8200794532637891559375.d0,  63777066403145711616000.d0 /
! fast return if possible
if (n.eq.0) then
  factnm=1.d0
  return
end if
if (m.eq.1) then
  if ((n.ge.1).and.(n.le.24)) then
    factnm=f1(n)
    return
  end if
end if
if (m.eq.2) then
  if ((n.ge.1).and.(n.le.38)) then
    factnm=f2(n)
    return
  end if
end if
if (n.lt.0) then
  write(*,*)
  write(*,'("Error(factnm): n < 0 : ",I8)') n
  write(*,*)
  stop
end if
if (m.le.0) then
  write(*,*)
  write(*,'("Error(factnm): m <= 0 : ",I8)') m
  write(*,*)
  stop
end if
if (n.gt.150) then
  write(*,*)
  write(*,'("Error(factnm): n out of range : ",I8)') n
  write(*,*)
  stop
end if
if (m.eq.1) then
  factnm=f1(24)
  do i=25,n
    factnm=factnm*dble(i)
  end do
else
  j=n/m
  if (mod(n,m).eq.0) j=j-1
  factnm=dble(n)
  do i=1,j
    factnm=factnm*dble(n-i*m)
  end do
end if
return
end function
!EOC
!======================================================================


end 	module modgaunt
