! This file is part of Cantera. See License.txt in the top-level directory or
! at https://cantera.org/license.txt for license and copyright information.

module cantera_kinetics

  use cantera_thermo
  use cantera_xml
  use fct

  contains

    subroutine ctkin_newFromFile(phase, filename, id, neighbor1, neighbor2, &
                                 neighbor3, neighbor4)
      implicit none
      type(phase_t), intent(inout) :: phase
      character*(*), intent(in) :: filename
      character*(*), intent(in), optional :: id
      type(phase_t), intent(in), optional :: neighbor1
      type(phase_t), intent(in), optional :: neighbor2
      type(phase_t), intent(in), optional :: neighbor3
      type(phase_t), intent(in), optional :: neighbor4
      integer :: n1, n2, n3, n4

      if (present(neighbor1)) then
         n1 = neighbor1%thermo_id
      else
         n1 = -1
      end if
      if (present(neighbor2)) then
         n2 = neighbor2%thermo_id
      else
         n2 = -1
      end if
      if (present(neighbor3)) then
         n3 = neighbor3%thermo_id
      else
         n3 = -1
      end if
      if (present(neighbor4)) then
         n4 = neighbor4%thermo_id
      else
         n4 = -1
      end if

      if (present(id)) then
          phase%kin_id = kin_newfromfile(filename, id, phase%thermo_id, &
                                         n1, n2, n3, n4)
      else
          phase%kin_id = kin_newfromfile(filename, '', phase%thermo_id, &
                                         n1, n2, n3, n4)
      end if
      phase%nrxn = kin_nreactions(phase%kin_id)
    end subroutine ctkin_newFromFile

    subroutine newKinetics(xml_phase, phase, &
         neighbor1, neighbor2, neighbor3, neighbor4)
      implicit none
      type(XML_Node), intent(in) :: xml_phase
      type(phase_t), intent(inout) :: phase
      type(phase_t), intent(in), optional :: neighbor1
      type(phase_t), intent(in), optional :: neighbor2
      type(phase_t), intent(in), optional :: neighbor3
      type(phase_t), intent(in), optional :: neighbor4
      integer :: missing
      missing = -1

      if (present(neighbor1)) then
         if (present(neighbor2)) then
            if (present(neighbor3)) then
               if (present(neighbor4)) then
                  phase%kin_id = newkineticsfromxml(xml_phase%xml_id, phase%thermo_id, &
                       neighbor1%thermo_id, neighbor2%thermo_id, neighbor3%thermo_id, &
                       neighbor4%thermo_id)
               else
                  phase%kin_id = newkineticsfromxml(xml_phase%xml_id, phase%thermo_id, &
                       neighbor1%thermo_id, neighbor2%thermo_id, neighbor3%thermo_id, &
                       missing)
               end if
            else
               phase%kin_id = newkineticsfromxml(xml_phase%xml_id, phase%thermo_id, &
                    neighbor1%thermo_id, neighbor2%thermo_id, missing, missing)
            end if
         else
            phase%kin_id = newkineticsfromxml(xml_phase%xml_id, phase%thermo_id, &
                 neighbor1%thermo_id, missing, missing, missing)
         end if
      else
         phase%kin_id = newkineticsfromxml(xml_phase%xml_id, phase%thermo_id, &
              missing, missing, missing, missing)
      end if
      phase%nrxn = kin_nreactions(phase%kin_id)
    end subroutine newKinetics

    subroutine ctkin_getKineticsType(self, nm)
      implicit none
      type(phase_t), intent(inout) :: self
      character*(*), intent(out) :: nm
      self%err = kin_gettype(self%kin_id, nm)
    end subroutine ctkin_getKineticsType

    integer function ctkin_kineticsStart(self, p)
      implicit none
      type(phase_t), intent(in) :: self
      integer, intent(in) :: p
      ctkin_kineticsStart = kin_start(self%kin_id, p)
    end function ctkin_kineticsstart

    integer function ctkin_kineticsSpeciesIndex(self, name, phase)
      implicit none
      type(phase_t), intent(inout) :: self
      character*(*), intent(in) :: name
      character*(*), intent(in) :: phase
      ctkin_kineticsSpeciesIndex = kin_speciesindex(self%kin_id, name, phase)
    end function ctkin_kineticsSpeciesIndex

    integer function ctkin_nTotalSpecies(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctkin_ntotalspecies = kin_ntotalspecies(self%kin_id)
    end function ctkin_ntotalspecies

    integer function ctkin_nReactions(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctkin_nreactions = kin_nreactions(self%kin_id)
    end function ctkin_nreactions

    integer function ctkin_nPhases(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctkin_nphases = kin_nphases(self%kin_id)
    end function ctkin_nphases

    integer function ctkin_phaseIndex(self, name)
      implicit none
      type(phase_t), intent(inout) :: self
      character*(*), intent(in) :: name
      ctkin_phaseindex = kin_phaseindex(self%kin_id, name)
    end function ctkin_phaseindex

    double precision function ctkin_reactantStoichCoeff(self, k, i)
      implicit none
      type(phase_t), intent(in) :: self
      integer, intent(in) :: k
      integer, intent(in) :: i
      ctkin_reactantstoichcoeff = kin_reactantstoichcoeff(self%kin_id, k, i)
    end function ctkin_reactantstoichcoeff

    double precision function ctkin_productStoichCoeff(self, k, i)
      implicit none
      type(phase_t), intent(in) :: self
      integer, intent(in) :: k
      integer, intent(in) :: i
      ctkin_productstoichcoeff = kin_productstoichcoeff(self%kin_id, k, i)
    end function ctkin_productstoichcoeff

    integer function ctkin_reactionType(self, i)
      implicit none
      type(phase_t), intent(in) :: self
      integer, intent(in) :: i
      ctkin_reactiontype = kin_reactiontype(self%kin_id, i)
    end function ctkin_reactiontype

    subroutine ctkin_getFwdRatesOfProgress(self, fwdROP)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(out) :: fwdROP(*)
      self%err = kin_getfwdratesofprogress(self%kin_id, fwdROP)
    end subroutine ctkin_getfwdratesofprogress

    subroutine ctkin_getRevRatesOfProgress(self, revROP)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(out) :: revROP(self%nrxn)
      self%err = kin_getrevratesofprogress(self%kin_id, revROP)
    end subroutine ctkin_getrevratesofprogress

    integer function ctkin_isReversible(self, i)
      implicit none
      type(phase_t), intent(inout) :: self
      integer, intent(in) :: i
      ctkin_isreversible = kin_isreversible(self%kin_id, i)
    end function ctkin_isreversible

    subroutine ctkin_getNetRatesOfProgress(self, netROP)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(out) :: netROP(*)
      self%err = kin_getnetratesofprogress(self%kin_id, netROP)
    end subroutine ctkin_getnetratesofprogress

    subroutine ctkin_getCreationRates(self, cdot)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(out) :: cdot(*)
      self%err = kin_getcreationrates(self%kin_id, cdot)
    end subroutine ctkin_getcreationrates

    subroutine ctkin_getDestructionRates(self, ddot)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(out) :: ddot(*)
      self%err = kin_getdestructionrates(self%kin_id, ddot)
    end subroutine ctkin_getdestructionrates

    subroutine ctkin_getNetProductionRates(self, wdot)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(out) :: wdot(*)
      self%err = kin_getnetproductionrates(self%kin_id, wdot)
    end subroutine ctkin_getnetproductionrates

    double precision function ctkin_multiplier(self, i)
      implicit none
      type(phase_t), intent(inout) :: self
      integer, intent(in) :: i
      ctkin_multiplier = kin_multiplier(self%kin_id, i)
    end function ctkin_multiplier

    subroutine ctkin_getEquilibriumConstants(self, kc)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(out) :: kc(*)
      self%err = kin_getequilibriumconstants(self%kin_id, kc)
    end subroutine ctkin_getequilibriumconstants

    subroutine ctkin_getReactionString(self, i, buf)
      implicit none
      type(phase_t), intent(inout) :: self
      integer, intent(in) :: i
      character*(*), intent(out) :: buf
      self%err = kin_getreactionstring(self%kin_id, i, buf)
    end subroutine ctkin_getreactionstring

    subroutine ctkin_setMultiplier(self, i, v)
      implicit none
      type(phase_t), intent(inout) :: self
      integer, intent(in) :: i
      double precision, intent(inout) :: v
      self%err = kin_setmultiplier(self%kin_id, i, v)
    end subroutine ctkin_setmultiplier

    subroutine ctkin_advanceCoverages(self, tstep)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(in) :: tstep
      self%err = kin_advancecoverages(self%kin_id, tstep)
    end subroutine ctkin_advancecoverages

end module cantera_kinetics
