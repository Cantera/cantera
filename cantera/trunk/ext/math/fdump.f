*DECK FDUMP
      SUBROUTINE FDUMP
C***BEGIN PROLOGUE  FDUMP
C***PURPOSE  Symbolic dump (should be locally written).
C***LIBRARY   SLATEC (XERMSG)
C***CATEGORY  R3
C***TYPE      ALL (FDUMP-A)
C***KEYWORDS  ERROR, XERMSG
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C        ***Note*** Machine Dependent Routine
C        FDUMP is intended to be replaced by a locally written
C        version which produces a symbolic dump.  Failing this,
C        it should be replaced by a version which prints the
C        subprogram nesting list.  Note that this dump must be
C        printed on each of up to five files, as indicated by the
C        XGETUA routine.  See XSETUA and XGETUA for details.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  FDUMP
C***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END

c$$$
c$$$      integer function isamax(n,sx,incx)
c$$$c
c$$$c     finds the index of element having max. absolute value.
c$$$c     jack dongarra, linpack, 3/11/78.
c$$$c     modified 3/93 to return if incx .le. 0.
c$$$c
c$$$      real sx(1),smax
c$$$      integer i,incx,ix,n
c$$$c
c$$$      isamax = 0
c$$$      if( n.lt.1 .or. incx.le.0 ) return
c$$$      isamax = 1
c$$$      if(n.eq.1)return
c$$$      if(incx.eq.1)go to 20
c$$$c
c$$$c        code for increment not equal to 1
c$$$c
c$$$      ix = 1
c$$$      smax = abs(sx(1))
c$$$      ix = ix + incx
c$$$      do 10 i = 2,n
c$$$         if(abs(sx(ix)).le.smax) go to 5
c$$$         isamax = i
c$$$         smax = abs(sx(ix))
c$$$    5    ix = ix + incx
c$$$   10 continue
c$$$      return
c$$$c
c$$$c        code for increment equal to 1
c$$$c
c$$$   20 smax = abs(sx(1))
c$$$      do 30 i = 2,n
c$$$         if(abs(sx(i)).le.smax) go to 30
c$$$         isamax = i
c$$$         smax = abs(sx(i))
c$$$   30 continue
c$$$      return
c$$$      end
c$$$
