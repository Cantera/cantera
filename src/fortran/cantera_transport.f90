! This file is part of Cantera. See License.txt in the top-level directory or
! at https://cantera.org/license.txt for license and copyright information.

module cantera_transport

  use cantera_thermo
  use fct

contains

!!$  integer function newTransport(model, phase, loglevel)
!!$    implicit none
!!$    type(XML_Node), intent(in) :: xml_phase
!!$    type(phase_t), intent(inout) :: phase
!!$    integer, optional :: loglevel
!!$    integer :: llevel
!!$    llevel = 0
!!$    if (present(loglevel)) then
!!$       llevel = loglevel
!!$    end if
!!$    phase%tran_id = newtransport(xml_phase%xml_id, phase%phase_id, loglevel)
!!$  end function newTransport

  double precision function ctrans_viscosity(self)
    implicit none
    type(phase_t), intent(inout) :: self
    ctrans_viscosity = trans_viscosity(self%tran_id)
    self%err = 0
  end function ctrans_viscosity

  double precision function ctrans_electricalConductivity(self)
    implicit none
    type(phase_t), intent(inout) :: self
    ctrans_electricalConductivity = trans_electricalConductivity(self%tran_id)
    self%err = 0
  end function ctrans_electricalConductivity

  double precision function ctrans_thermalConductivity(self)
    implicit none
    type(phase_t), intent(inout) :: self
    ctrans_thermalConductivity = trans_thermalConductivity(self%tran_id)
    self%err = 0
  end function ctrans_thermalConductivity

  subroutine ctrans_getThermalDiffCoeffs(self, dt)
    implicit none
    type(phase_t), intent(inout) :: self
    double precision, intent(out) :: dt(*)
    self%err = trans_getThermalDiffCoeffs(self%tran_id, dt)
  end subroutine ctrans_getThermalDiffCoeffs

    subroutine ctrans_getMixDiffCoeffs(self, d)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(out) :: d(*)
      self%err = trans_getMixDiffCoeffs(self%tran_id, d)
    end subroutine ctrans_getMixDiffCoeffs

    subroutine ctrans_getMixDiffCoeffsMass(self, d)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(out) :: d(*)
      self%err = trans_getMixDiffCoeffsMass(self%tran_id, d)
    end subroutine ctrans_getMixDiffCoeffsMass

    subroutine ctrans_getMixDiffCoeffsMole(self, d)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(out) :: d(*)
      self%err = trans_getMixDiffCoeffsMole(self%tran_id, d)
    end subroutine ctrans_getMixDiffCoeffsMole

    subroutine ctrans_getBinDiffCoeffs(self, ld, d)
      implicit none
      type(phase_t), intent(inout) :: self
      integer, intent(in) :: ld
      double precision, intent(out) :: d(*)
      self%err = trans_getBinDiffCoeffs(self%tran_id, ld, d)
    end subroutine ctrans_getBinDiffCoeffs

    subroutine ctrans_getMultiDiffCoeffs(self, ld, d)
      implicit none
      type(phase_t), intent(inout) :: self
      integer, intent(in) :: ld
      double precision, intent(out) :: d(*)
      self%err = trans_getMultiDiffCoeffs(self%tran_id, ld, d)
    end subroutine ctrans_getMultiDiffCoeffs

    subroutine ctrans_setParameters(self, type, k, d)
      implicit none
      type(phase_t), intent(inout) :: self
      integer, intent(in) :: type
      integer, intent(in) :: k
      double precision, intent(in) :: d(*)
      self%err = trans_setParameters(self%tran_id, type, k, d)
    end subroutine ctrans_setParameters

end module cantera_transport
