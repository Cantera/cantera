! ------------------------------------------------------------------------
! 
! Copyright 2002 California Institute of Technology
! 
! Version 1.2.0.1
! Mon Jan  7 07:13:18 2002
! This file generated automatically.
! Fortran 90 module implementing Cantera class Transport
! ------------------------------------------------------------------------
      module cttransport
      use cttypes
      interface

      integer function ctrans_Transport(mix_hndl, 
     &                           type, trans_db, logLevel)
!DEC$ attributes c, reference ::  ctrans_Transport
!DEC$ attributes alias: '_ctrans_transport_' ::  ctrans_Transport
!DEC$ attributes dllimport ::  ctrans_Transport
      integer, intent(inout) :: mix_hndl
      integer, intent(in) :: type
      character*(*), intent(in) :: trans_db
      integer, intent(in) :: logLevel
      end function
      subroutine ctrans_delete(transportMgr)
!DEC$ attributes c, reference ::  ctrans_delete
!DEC$ attributes alias: '_ctrans_delete_' ::  ctrans_delete
!DEC$ attributes dllimport ::  ctrans_delete
      type(Transport), intent(inout) :: transportMgr
      end subroutine


      subroutine ctrans_setopt(itr, ijob,
     &         iopt, dopt)
!DEC$ attributes c, reference ::  ctrans_setopt
!DEC$ attributes alias: '_ctrans_setopt_' ::  ctrans_setopt
!DEC$ attributes dllimport ::  ctrans_setopt
      integer, intent(inout) :: itr
      integer, intent(in) :: ijob
      integer, intent(in) :: iopt
      double precision, intent(in) :: dopt
      end subroutine
      end interface
      contains

      type(Transport) function MultiTransport(mix, 
     &                                     trans_db, logLevel)
      type(Mixture), intent(inout) :: mix
      character*(*), intent(in) :: trans_db
      integer, intent(in) :: logLevel
      type(Transport) tr
      tr%hndl=ctrans_transport(mix%hndl, 
     &                  Multicomponent, trans_db, logLevel)
      tr%ierr = 0
      MultiTransport = tr
      return
      end function


      type(Transport) function MixTransport(mix, 
     &                                     trans_db, logLevel)
      type(Mixture), intent(inout) :: mix
      character*(*), intent(in) :: trans_db
      integer, intent(in) :: logLevel
      type(Transport) tr
      tr%hndl=ctrans_transport(mix%hndl, 
     &                  MixtureAveraged, trans_db, logLevel)
      tr%ierr = 0
      MixTransport = tr
      return
      end function


      subroutine trans_delete(transportMgr)
      type(Transport), intent(inout) :: transportMgr
      call ctrans_delete(transportMgr%hndl)
      transportMgr%hndl = 0
      transportMgr%ierr = 0
      return
      end subroutine

      subroutine trans_setopt(tr, linearSolver, GMRES_m, GMRES_eps)
      implicit double precision (a-h,o-z)
      type(Transport), intent(inout) :: tr
      character*(*), intent(in), optional :: linearSolver
      integer, intent(in), optional :: GMRES_m
      double precision, intent(in), optional :: GMRES_eps
      itr = tr%hndl
      dummy = -1.d0
      idummy = -1
      if (present(linearSolver)) then
          meth = 2
          if (linearSolver .eq. 'GMRES') meth = 1
          ijob = 0
          write(*,*) 'setting lin solver'
          call ctrans_setopt(itr,ijob, meth, dummy)
          write(*,*) 'ret from setting lin solver'
      end if
      if (present(GMRES_m)) then
          ijob = 1
          call ctrans_setopt(itr,ijob, GMRES_m, dummy)
      end if
      if (present(GMRES_eps)) then
          ijob = 2
          call ctrans_setopt(itr,ijob,  idummy, GMRES_eps)
      end if
      return
      end subroutine
      end module
