! ------------------------------------------------------------------------
! 
! Copyright 2002 California Institute of Technology
! 
! Version 1.2.0.1
! Mon Jan  7 07:13:18 2002
! This file generated automatically.
! Fortran 90 module implementing Cantera class Utils
! ------------------------------------------------------------------------
      module ctutils
      use cttypes
      interface

      subroutine cutil_report(mix_hndl, n, txt)
!DEC$ attributes c, reference ::  cutil_report
!DEC$ attributes alias: '_cutil_report_' ::  cutil_report
!DEC$ attributes dllimport ::  cutil_report
      integer, intent(in) :: mix_hndl
      integer, intent(in) :: n
      character*(*), intent(out) :: txt
      end subroutine
      subroutine cutil_addDirectory(dirname)
!DEC$ attributes c, reference ::  cutil_addDirectory
!DEC$ attributes alias: '_cutil_adddirectory_' ::  cutil_addDirectory
!DEC$ attributes dllimport ::  cutil_addDirectory
      character*(*), intent(in) :: dirname
      end subroutine

      end interface
      contains

      subroutine util_printSummary(mix, lu)
      type(Mixture),  intent(in)  :: mix
      integer, intent(in) :: lu
      character*4000 txt
      n = 4000 
      call cutil_report(mix%hndl, n, txt)
      istrt = 1
      ifin = 100
      do i = 1,40
         write(lu,'(a$)') txt(istrt:ifin)
         istrt = i*100 + 1
         ifin = istrt + 99
         if (ifin .gt. n) ifin = n
         if (istrt .ge. n) return
      end do
      return
      end subroutine
!
!     addDirectory
!
      subroutine util_addDirectory(dirname)
      character*(*), intent(in) :: dirname
      call cutil_addDirectory(dirname)
      return
      end subroutine
      end module
