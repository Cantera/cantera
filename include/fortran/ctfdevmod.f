! ------------------------------------------------------------------------
! 
! Copyright 2002 California Institute of Technology
! 
! Version 1.2.0.1
! Mon Jan  7 07:13:18 2002
! This file generated automatically.
! Fortran 90 module implementing Cantera class FlowDev
! ------------------------------------------------------------------------
      module ctfdev
      use cttypes

      use ctreactor
      interface

      integer function cfdev_newmfc()
!DEC$ attributes c, reference :: cfdev_newmfc
!DEC$ attributes alias:'_cfdev_newmfc_' :: cfdev_newmfc
!DEC$ attributes dllimport :: cfdev_newmfc
      end function

      integer function cfdev_newpc()
!DEC$ attributes c, reference :: cfdev_newpc
!DEC$ attributes alias:'_cfdev_newpc_' :: cfdev_newpc
!DEC$ attributes dllimport :: cfdev_newpc
      end function

      subroutine cfdev_delete(i) 
!DEC$ attributes c, reference :: cfdev_delete
!DEC$ attributes alias:'_cfdev_delete_' :: cfdev_delete
!DEC$ attributes dllimport :: cfdev_delete
      integer, intent(in) :: i
      end subroutine

      subroutine cfdev_setspnt(flowDev, setpnt)
!DEC$ attributes c, reference ::  cfdev_setspnt
!DEC$ attributes alias: '_cfdev_setspnt_' ::  cfdev_setspnt
!DEC$ attributes dllimport ::  cfdev_setspnt
      type(Flowdev), intent(inout) :: flowDev
      double precision, intent(in) :: setpnt
      end subroutine


      subroutine cfdev_install(flowDev_hndl, up_hndl, down_hndl)
!DEC$ attributes c, reference ::  cfdev_install
!DEC$ attributes alias: '_cfdev_install_' ::  cfdev_install
!DEC$ attributes dllimport ::  cfdev_install
      integer, intent(inout) :: flowDev_hndl
      integer, intent(inout) :: up_hndl
      integer, intent(inout) :: down_hndl
      end subroutine
      type(Reactor) function cfdev_upstream(flowDev)
!DEC$ attributes c, reference ::  cfdev_upstream
!DEC$ attributes alias: '_cfdev_upstream_' ::  cfdev_upstream
!DEC$ attributes dllimport ::  cfdev_upstream
      type(Flowdev), intent(inout) :: flowDev
      end function

      type(Reactor) function cfdev_downstream(flowDev)
!DEC$ attributes c, reference ::  cfdev_downstream
!DEC$ attributes alias: '_cfdev_downstream_' ::  cfdev_downstream
!DEC$ attributes dllimport ::  cfdev_downstream
      type(Flowdev), intent(inout) :: flowDev
      end function

      double precision function cfdev_mfrate(flowDev)
!DEC$ attributes c, reference ::  cfdev_mfrate
!DEC$ attributes alias: '_cfdev_mfrate_' ::  cfdev_mfrate
!DEC$ attributes dllimport ::  cfdev_mfrate
      type(Flowdev), intent(inout) :: flowDev
      end function

      subroutine cfdev_update(flowDev)
!DEC$ attributes c, reference ::  cfdev_update
!DEC$ attributes alias: '_cfdev_update_' ::  cfdev_update
!DEC$ attributes dllimport ::  cfdev_update
      type(Flowdev), intent(inout) :: flowDev
      end subroutine

      subroutine cfdev_reset(flowDev)
!DEC$ attributes c, reference ::  cfdev_reset
!DEC$ attributes alias: '_cfdev_reset_' ::  cfdev_reset
!DEC$ attributes dllimport ::  cfdev_reset
      type(Flowdev), intent(inout) :: flowDev
      end subroutine

      integer function cfdev_ready(flowDev)
!DEC$ attributes c, reference ::  cfdev_ready
!DEC$ attributes alias: '_cfdev_ready_' ::  cfdev_ready
!DEC$ attributes dllimport ::  cfdev_ready
      type(Flowdev), intent(inout) :: flowDev
      end function


      integer function cfdev_setGains(flowDev_hndl, n, gains)
!DEC$ attributes c, reference ::  cfdev_setGains
!DEC$ attributes alias: '_cfdev_setgains_' ::  cfdev_setGains
!DEC$ attributes dllimport ::  cfdev_setGains
      integer, intent(inout) :: flowDev_hndl
      integer, intent(in) :: n
      double precision, intent(in) :: gains(n)
      end function

      integer function cfdev_getGains(flowDev_hndl, n, gains)
!DEC$ attributes c, reference ::  cfdev_getGains
!DEC$ attributes alias: '_cfdev_getgains_' ::  cfdev_getGains
!DEC$ attributes dllimport ::  cfdev_getGains
      integer, intent(inout) :: flowDev_hndl
      integer, intent(in) :: n
      double precision, intent(out) :: gains(n)
      end function
      double precision function cfdev_maxError(flowDev)
!DEC$ attributes c, reference ::  cfdev_maxError
!DEC$ attributes alias: '_cfdev_maxerror_' ::  cfdev_maxError
!DEC$ attributes dllimport ::  cfdev_maxError
      type(Flowdev), intent(inout) :: flowDev
      end function

      end interface
      integer, parameter :: MassFlowController_Type = constant: requested text line not found.D
      integer, parameter :: PressureController_Type = constant: requested text line not found.D
      contains

      type(Flowdev) function MassFlowController()
      type(Flowdev) mfc
      mfc%hndl = cfdev_newmfc()
      mfc%type   = MassFlowController_Type;
      MassFlowController = mfc
      return
      end function

      type(Flowdev) function PressureController()
      type(Flowdev) pc
      pc%hndl = cfdev_newpc()
      pc%type   = PressureController_Type;
      PressureController = pc
      return
      end function
 
      subroutine fdev_copy(dest, src)
      type (Flowdev), intent(out) :: dest
      type (Flowdev), intent(in) :: src
      call reac_copy(dest%upstream, src%upstream)
      call reac_copy(dest%downstream, src%downstream)
      dest%hndl = src%hndl
      dest%type = src%type
      end subroutine
!
!     setSetpoint
!
      subroutine fdev_setspnt(self, flowDev, setpnt)
      type(FlowDev), intent(inout) :: self
      type(Flowdev), intent(inout) :: flowDev
      double precision, intent(in) :: setpnt
      call cfdev_setspnt(self%hndl, flowDev, setpnt)
      return
      end subroutine

      subroutine fdev_install(object, upstream, downstream)
      type(Flowdev), intent(inout) :: object
      type(Reactor), intent(inout) :: upstream
      type(Reactor), intent(inout) :: downstream
      call cfdev_install(object%hndl, upstream%hndl,
     &     downstream%hndl)
      call reac_copy(object%upstream, upstream)
      call reac_copy(object%downstream, downstream)
      return
      end subroutine

      type(Reactor) function fdev_upstream(object)
      type(Flowdev), intent(in) :: object
      call reac_copy(fdev_upstream, object%upstream)
      end function

      type(Reactor) function fdev_downstream(object)
      type(Flowdev), intent(in) :: object
      call reac_copy(fdev_downstream, object%downstream)
      end function
!
!     massFlowRate
!
      double precision function fdev_mfrate(self, flowDev)
      type(FlowDev), intent(inout) :: self
      type(Flowdev), intent(inout) :: flowDev
      fdev_mfrate=cfdev_mfrate(self%hndl, flowDev)
      return
      end function
!
!     update
!
      subroutine fdev_update(self, flowDev)
      type(FlowDev), intent(inout) :: self
      type(Flowdev), intent(inout) :: flowDev
      call cfdev_update(self%hndl, flowDev)
      return
      end subroutine
!
!     reset
!
      subroutine fdev_reset(self, flowDev)
      type(FlowDev), intent(inout) :: self
      type(Flowdev), intent(inout) :: flowDev
      call cfdev_reset(self%hndl, flowDev)
      return
      end subroutine

      logical function fdev_ready(object)
      type(Flowdev), intent(inout) :: object
      if (cfdev_ready(object%hndl) .gt. 0) then
         fdev_ready = .true.
      else
         fdev_ready = .false.
      end if
      return
      end function

      subroutine fdev_setGains(flowDev, gains)
      type(Flowdev), intent(inout) :: flowDev
      double precision, intent(in) :: gains(4)
      iok = cfdev_setGains(flowDev%hndl, 4, gains)
      if (iok .eq. 0) write(*,*) 'error setting gains!'
      return
      end subroutine

      subroutine fdev_getGains(flowDev, gains)
      type(Flowdev), intent(inout) :: flowDev
      double precision, intent(out) :: gains(4)
      iok = cfdev_getGains(flowDev%hndl, 4, gains)
      if (iok .eq. 0) write(*,*) 'error getting gains!'
      return
      end subroutine
!
!     maxError
!
      double precision function fdev_maxError(self, flowDev)
      type(FlowDev), intent(inout) :: self
      type(Flowdev), intent(inout) :: flowDev
      fdev_maxError=cfdev_maxError(self%hndl, flowDev)
      return
      end function
      end module
