	module mscg
	implicit none
	contains
	

      subroutine scg(lmax,cg,indxc,l3cg)
C- Gaunt coefficients for real harmonics.  See scg0
C ----------------------------------------------------------------------
Ci Inputs
Ci   lmax  :maximum l
Co Outputs
Co   indxc : For a given L1, L2 (note L = compound lm index)
Co         : indxc holds limits ic1,ic2 for cg and L3, below.
Co         : Let indx = (ii*(ii-1))/2+min(L1,L2), ii=max(L1,L2)
Co         : Then icg1=indxc(indx) and icg2=indxc(indx+1)-1
Co   cg    : Gaunt coefficients strung together as one vector.
Co         : cg(ic1:ic2) = <L1 L2 | L3>, ic1,ic2 span all L3
Co           for which there is a nonvanishing <L1 L2 | L3> (see indxc)
Co   l3cg  :l3cg(ic1:ic2) is L3 corresponding to cg(ic1:ic2)
Cr Remarks
Cr   Calls scg0 in mode=1, which see
Cu Updates
Cu   23 Dec 17
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lmax
      integer indxc(*),l3cg(*)
      double precision cg(*)
C ... Local parameters
      integer lnjcg,lnxcg

      call scg0(1,lmax,cg,indxc,l3cg,lnjcg,lnxcg)
      end

      subroutine scg0(mode,lmax,cg,indxc,l3cg,lnjcg,lnxcg)
C- Gaunt coefficients for real harmonics
C ----------------------------------------------------------------
Ci Inputs
Ci   mode  :0 return lnjcg,lnxcg only
Ci         :1 Also return cg,indxc,l3cg
Ci   lmax  :maximum l
Co Outputs
Co   indxc : For a given L1, L2 (note L = compound lm index)
Co         : indxc holds limits ic1,ic2 for cg and L3, below.
Co         : Let indx = (ii*(ii-1))/2+min(L1,L2), ii=max(L1,L2)
Co         : Then icg1=indxc(indx) and icg2=indxc(indx+1)-1
Co   cg    : Gaunt coefficients strung together as one vector.
Co         : cg(ic1:ic2) = <L1 L2 | L3>, ic1,ic2 span all L3
Co           for which there is a nonvanishing <L1 L2 | L3> (see indxc)
Co   l3cg  :l3cg(ic1:ic2) is L3 corresponding to cg(ic1:ic2)
Co   lnjcg :number of coefficients cg
Co   lnxcg :number of coefficients in indxc
Cr Remarks
Cr   The following segment expands cg into array form cga(L3,L1,L2)
Cr      do  L1 = 1, nlm1
Cr        do  L2 = 1, nlm2
Cr          ii = max0(L1,L2)
Cr          indx = (ii*(ii-1))/2 + min0(L1,L2)
Cr          ic1 = indxc(indx)
Cr          ic2 = indxc(indx+1)-1
Cr          do  icg = ic1, ic2
Cr            L3 = jcg(icg)
Cr            cga(L3,L1,L2) = cg(icg)
Cr          enddo
Cr        enddo
Cr      enddo
Cr
Cr   (FORMERLY S104 IN ASW)
Cu Updates
Cu   23 Dec 17  (MvS) adapted so that it can return lnjcg,lnxcg for dimensioning
Cu   15 Mar 15  (DP)  Revised documentation
Cu   26 Sep 14  (MvS) Increased precision of sqrt(pi)
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,lmax,lnjcg,lnxcg
      integer indxc(*),l3cg(*)
      double precision cg(*)
C ... Local parameters
      integer i,i1,i2,i3,i31,i32,ic,j1,j1s,j2,j2s,k2,l1,l2,l3,lmindx,
     .        m1,m2,m3,mb,n1,n2,n3,nl,nm3,s1,s2,s3,t1,t2,t3
      double precision q1,sr2,t,srpi,fs !,f100,f102
      double precision fac(161)
      !external f100,f102
C     data srpi /1.772453851d0/  ! Prior to Sep 2014
      data srpi /1.7724538509055159d0/
      fs(i) = 1 + 4*(i/2) - 2*i

      lnjcg = 0
      lnxcg = 0

      mb = 999999
      nl = lmax+1
      sr2 = dsqrt(2d0)
      fac(1) = 1d0
      do  i = 1, 160
        fac(i+1) = i*fac(i)
      enddo
      ic = 0
      lmindx = 0
      do  i1 = 1, nl
        l1 = i1-1
        j1s = 2*l1+1
        do  j1 = 1, j1s
          m1 = j1-i1
          n1 = iabs(m1)
          s1 = 0
          if (m1 < 0) s1 = 1
          t1 = 0
          if (m1 == 0) t1 = 1
          do  i2 = 1, i1
            l2 = i2-1
            i31 = l1 - l2 + 1
            i32 = l1 + l2 + 1
            j2s = 2*l2 + 1
            k2 = j1s*j2s
            if (i2 == i1) j2s = j1
            do  j2 = 1, j2s
              lmindx = lmindx + 1
              if (mode == 1) indxc(lmindx) = ic+1
              m2 = j2-i2
              n2 = iabs(m2)
              s2 = 0
              if (m2 < 0) s2 = 1
              t2 = 0
              if (m2 == 0) t2 = 1
              if (m1*m2 < 0) then
                m3 = -n1-n2
                mb = -iabs(n1-n2)
                if (mb == 0) then
                  nm3 = 1
                else
                  nm3 = 2
                endif
              elseif (m1*m2 == 0) then
                m3 = m1+m2
                nm3 = 1
              else
                m3 = n1+n2
                mb = iabs(n1-n2)
                nm3 = 2
              endif
    1         continue
              n3 = iabs(m3)
              s3 = 0
              if (m3 < 0) s3 = 1
              t3 = 0
              if (m3 == 0) t3 = 1
              q1 = dsqrt(dble(k2))
     .           *fs(n3+(s1+s2+s3)/2)/(2*sr2**(1+t1+t2+t3))
              do  i3 = i31, i32, 2
                l3 = i3-1
                if (n3 > l3) cycle
                if (mode == 1) then
                t = 0d0
                if (n1+n2 == -n3) t = t + f102(fac,l1,l2,l3)
                if (n1+n2 == n3)
     .            t = t + f100(fac,l1,l2,l3,n1,n2,n3)*fs(n3+s3)
                if (n1-n2 == -n3)
     .            t = t + f100(fac,l1,l2,l3,n1,-n2,-n3)*fs(n2+s2)
                if (n1-n2 == n3)
     .            t = t + f100(fac,l1,l2,l3,-n1,n2,-n3)*fs(n1+s1)
                endif
                ic = ic+1
                lnjcg = max(ic,lnjcg)
                if (mode == 1) then
	       cg(ic) = q1*t*f102(fac,l1,l2,l3)/(srpi*dsqrt(dble(2*l3+1)))
                l3cg(ic) = l3*(l3+1) + m3 + 1
                endif
              enddo
              nm3 = nm3-1
              m3 = mb
              if (nm3 > 0) goto 1
            enddo
          enddo
        enddo
      enddo
      if (mode == 1) indxc(lmindx+1) = ic+1
      lnxcg = max(lnxcg,lmindx+1)
      end

      double precision function f100(fac,j1,j2,j3,m1,m2,m3)
C- Clebsch-Gordan coefficients
C ----------------------------------------------------------------
Ci Inputs
Ci   FAC,J1,J2,J3,M1,M2,M3
Co Outputs
Co   F100
Cr Remarks
Cr   selection rules: m3 == m1+m2
Cr                    |j2-j1| <= j3 <= j1+j2 (not checked here, correct input expected!)
Cr   <j1 m1, j2 m2| j3 m3> = (2*j3+1)
Cr       *sqrt((j1+j2-j3)!(j1-j2+j3)!(-j1+j2+j3)!/(j1+j2+j3+1)!)
Cr       *sqrt((j1-m1)!(j1+m1)!(j2-m2)!(j2+m2)!(j2-m2)!(j3+m3)!)
Cr       *sum_{k=(max(j2-j3-m1,j1-j3+m2,0):min(j1+j2-j3,j1-m1,j2+m2))}
Cr            ((-1)**k/(k!(j1+j2-j3-k)!(j1-m1-k)!(j2+m2-k)!(j3-j2+m1+k)!(j3-j1-m2+k)!))
Cu Updates
Cu   15 Mar 15  (DP)  Revised documentation
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer j1,j2,j3,m1,m2,m3
      double precision fac(50)
C Local parameters
      integer m,n,n1,n2
      double precision t,t1

      f100 = 0d0
      if (m3 /= m1+m2) return
      t = (2*j3+1)*fac(j1+j2-j3+1)*fac(j3+j1-j2+1)*fac(j3+j2-j1+1)/
     .    fac(j1+j2+j3+2)
      t = dsqrt(t*fac(j1+m1+1)*fac(j1-m1+1)*fac(j2+m2+1)*fac(j2-m2+1)
     .    *fac(j3+m3+1)*fac(j3-m3+1))
      n1 = max0(j2-j3-m1,j1-j3+m2,0) + 1
      n2 = min0(j1+j2-j3,j1-m1,j2+m2) + 1
      if (n1 > n2) return
      t1 = 0d0
      do   m = n1, n2
        n = m-1
        t1 = t1 +
     .    dble(1+4*(n/2)-2*n)/(fac(m)*fac(j1+j2-j3-n+1)*fac(j1-m1-n+1)
     .    *fac(j2+m2-n+1)*fac(j3-j2+m1+n+1)*fac(j3-j1-m2+n+1))
      enddo
      f100 = t*t1
      end

      double precision function f102(fac,l1,l2,l3)
C- Clebsch-Gordan coefficients for m1==m2==m3==0 (simplified expression).
C ----------------------------------------------------------------
Ci Inputs
Ci  fac,l1,l2,l3
Co Outputs
Co   f102
Cr Remarks
Cr   selection rules: l1+l2+l3 == 2n
Cr                    |l2-l1| <= l3 <= l1+l2 (not checked here, correct input expected!)
Cr   <l1 0, l2 0| l3 0>
Cu Updates
Cu   15 Mar 15  (DP)  Revised documentation
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer l1,l2,l3
      double precision fac(50)
C Local parameters
      integer lt,p,x

      lt = l1 + l2 + l3
      p = lt/2
      f102 = 0d0
      if (2*p /= lt) return
      f102 = dsqrt(dble(2*l3+1)/dble(lt+1))
      f102 = f102*fac(p+1)/dsqrt(fac(2*p+1))
      x = p-l1
      f102 = f102*dsqrt(fac(2*x+1))/fac(x+1)
      x = p-l2
      f102 = f102*dsqrt(fac(2*x+1))/fac(x+1)
      x = p-l3
      f102 = f102*dsqrt(fac(2*x+1))/fac(x+1)
      if (x > 2*(x/2)) f102 = -f102
      end
	end 	module mscg
