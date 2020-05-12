      logical function isanrg(i,i1,i2,t1,t2,lreqd)
C- Sanity check for integer value
C ----------------------------------------------------------------------
Ci Inputs
Ci   i     :quantity to be checked
Ci   i1    :test passes if i>=i1
Ci   i2    :test passes if i<=i2
Ci   t1    :first part of error message string, if check fails
Ci         :Note: if t1 begins with string 'file:' that part is
Ci         :stripped, and the printout message reads
Ci         :'unexpected file value' rather than
Ci   t2    :second part of error message string, if check fails
Ci   lreqd :F, error message printed as warning; isanrg returns
Ci         :T, error message printed; isanrg aborts
Co Outputs
Co   message output to lgunit(1) if check fails
Co   isanrg:(returns only if lreqd=F)
Co         :.false. if test passes
Co         :.true. if test fails
Cr Remarks
Cr   Error or or warning message reads (single value)
Cr    t1 'unexpected value (i) for' t2 'expected (i)'
Cr   or
Cr    t1 'unexpected file value (i) for' t2 'expected (i)'
Cr   In case of range of values the warning message reads one of
Cr    t1 'unexpected value (i) for' t2 '(valid range (i1) to (i2))'
Cr    t1 'unexpected file value (i) for' t2 '(valid range (i1) to (i2))'
Cr   Special case i1>i2 and i2 is NULLI:  message prints one of the following:
Cr    t1 'unexpected value (i) for' t2 'expected at least (i)'
Cr    t1 'unexpected file value (i) for' t2 'expected at least (i)'
Cu Updates
Cu   05 May 12 Extra functionality in message printout: i2 = NULLI
Cu   04 Aug 06 Extra functionality in message printout
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical lreqd
      integer i,i1,i2
      character*(*) t1,t2
C ... Local parameters
      integer lgunit,iprint,k1,k2,it1,NULLI
      parameter (NULLI=-99999)
      character strn*120,strn2*120,t3*30
      isanrg = .false.
      if (i >= i1 .and. i <= i2) return
      if (i >= i1 .and. i2 == NULLI) return
      isanrg = .true.
      it1 = 1
      t3 = ' unexpected value %i for %-1f'
      if (t1(1:5) == 'file:') then
        it1 = 6
        t3 = ' unexpected file value %i for '
      endif
      if (i1 == i2) then
        strn = t2 // ' ... expected %i'
      elseif (i1 > i2 .and. i2 == NULLI) then
        strn = t2 // ' ... expected at least %i'
      else
        strn = t2 // ' (valid range %i to %i)'
      endif
      strn2 = t1(it1:) // t3 // strn
      call awrit3(strn2,strn,len(strn),0,i,i1,i2)
      if (lreqd) then
        call strip(strn,k1,k2)
        call setpr(30)
        call rx(strn(k1:k2))
      endif
      if (iprint() > 0) call awrit0(strn,' ',-len(strn),lgunit(1))
      end

      subroutine fsanrg(f,f1,f2,tol,t1,t2,lreqd)
C- Sanity check for double precision value
C ----------------------------------------------------------------------
Ci Inputs
Ci   f     :quantity to be checked
Ci   f1    :(case f2-f1>0) test passes if f>=f1
Ci         :(case f2-f1=0) test passes if f>=f1-tol/2
Ci   f2    :(case f2-f1>0) test passes if f<=f2
Ci         :(case f2-f1=0) test passes if f<=f2+tol/2
Ci   tol   :allowed tolerance in f, if bounds f2-f1=0
Ci   t1    :first part of error message string, if check fails
Ci   t2    :second part of error message string, if check fails
Ci   lreqd :F, error message printed as warning; fsanrg returns
Ci         :T, error message printed; fsanrg aborts
Co Outputs
Co   message output to lgunit(1) if check fails
Cu Updates
Cu   06 Aug 01 Added argument tol (altered argument list)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical lreqd
      double precision f,f1,f2,tol
      character*(*) t1,t2
C ... Local parameters
      character strn*120,strn2*120,t3*31
      integer lgunit,k1,k2,iprint,it1

      if (f >= f1 .and. f <= f2) return
      if (f1 == f2 .and. f >= f1-tol/2 .and. f <= f2+tol/2) return

      it1 = 1
      t3 = ' unexpected value %g for%f'
      if (len(t2) == 0) t3 = ' unexpected value %g'
      if (t1(1:5) == 'file:') then
        it1 = 6
        t3 = ' unexpected file value %g for%f'
        if (len(t2) == 0) t3 = ' unexpected file value %g'
      endif
      if (f1 == f2) then
        strn = t2 // ' ... expected %g'
      else
        strn = t2 // ' (valid range %g to %g)'
      endif
      strn2 = t1(it1:) // trim(t3) // strn
      call awrit3(strn2,strn,len(strn),0,f,f1,f2)
      if (lreqd) then
        call strip(strn,k1,k2)
        call setpr(30)
        call rx(strn(k1:k2))
      endif
      if (iprint() > 0) call awrit0(strn,' ',-len(strn),lgunit(1))
      end

      subroutine sanrg(cond,i,i1,i2,t1,t2)
C- Conditional sanity check for integer value
C ----------------------------------------------------------------------
Cr Remarks
Cr   If cond is false, this routine does nothing.
Cr   If cond is true, this routine calls isanrg with lreqd = .true.
Cu Updates
Cu   13 Sep 17 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical cond
      integer i,i1,i2
      character*(*) t1,t2
C ... Local parameters
      logical, parameter:: lreqd = .true.
      logical ltmp
      procedure(logical) isanrg

      if (.not. cond) return
      ltmp = isanrg(i,i1,i2,t1,t2,lreqd)

      end
