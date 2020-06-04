! This file is part of Cantera. See License.txt in the top-level directory or
! at https://cantera.org/license.txt for license and copyright information.

module cantera_thermo

  use fct
  use cantera_xml

  type phase_t
    integer :: thermo_id
    integer :: kin_id
    integer :: tran_id
    integer :: err
    integer :: nel
    integer :: nsp
    integer :: nrxn
  end type phase_t

! these definitions are for use with the equilibrate function.
!  integer, parameter :: TV = 100
!  integer, parameter :: HP = 101
!  integer, parameter :: SP = 102
!  integer, parameter :: PV = 103
!  integer, parameter :: TP = 104
!  integer, parameter :: UV = 105
!  integer, parameter :: SV = 107

!  integer, parameter :: VT = -100
!  integer, parameter :: PH = -101
!  integer, parameter :: PS = -102
!  integer, parameter :: VP = -103
!  integer, parameter :: PT = -104
!  integer, parameter :: VU = -105
!  integer, parameter :: VS = -107

contains

    type(phase_t) function ctthermo_newFromFile(filename, id)
      implicit none
      character*(*), intent(in) :: filename
      character*(*), intent(in), optional :: id
      type(phase_t) :: self
      if (present(id)) then
         self%thermo_id = th_newfromfile(filename, id)
      else
         self%thermo_id = th_newfromfile(filename, '')
      end if
      self%nel = phase_nelements(self%thermo_id)
      self%nsp = phase_nspecies(self%thermo_id)
      self%nrxn = 0
      self%err = 0
      self%kin_id = -1
      self%tran_id = -1
      ctthermo_newFromFile = self
    end function ctthermo_newFromFile

    type(phase_t) function newThermoPhase(xml_phase, index)
      implicit none
      type(XML_Node), intent(inout), optional :: xml_phase
      integer, intent(in), optional :: index
      type(phase_t) :: self
      if (present(index)) then
         self%thermo_id = index
         self%err = 0
      else if (present(xml_phase)) then
         self%thermo_id = newthermofromxml(xml_phase%xml_id)
         self%nel = phase_nelements(self%thermo_id)
         self%nsp = phase_nspecies(self%thermo_id)
         self%nrxn = 0
         self%err = 0
         self%kin_id = -1
         self%tran_id = -1
      else
         call cantera_error('newThermoPhase','xml_phase or id must be specified.')
      end if
      newThermoPhase = self
    end function newThermoPhase

    subroutine ctthermo_getName(self, name)
      implicit none
      type(phase_t), intent(inout) :: self
      character*(*), intent(out) :: name
      call phase_getname(self%thermo_id, name)
    end subroutine ctthermo_getName

    integer function ctthermo_nElements(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctthermo_nelements = phase_nelements(self%thermo_id)
    end function ctthermo_nElements

    integer function ctthermo_nSpecies(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctthermo_nspecies = phase_nspecies(self%thermo_id)
    end function ctthermo_nSpecies

    double precision function ctthermo_temperature(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctthermo_temperature = phase_temperature(self%thermo_id)
    end function ctthermo_temperature

    subroutine ctthermo_setTemperature(self, t)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(in) :: t
      self%err = phase_settemperature(self%thermo_id, t)
    end subroutine ctthermo_setTemperature

    double precision function ctthermo_density(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctthermo_density = phase_density(self%thermo_id)
    end function ctthermo_density

    subroutine ctthermo_setDensity(self, rho)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(in) :: rho
      self%err = phase_setdensity(self%thermo_id, rho)
    end subroutine ctthermo_setDensity

    double precision function ctthermo_molarDensity(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctthermo_molardensity = phase_molardensity(self%thermo_id)
    end function ctthermo_molarDensity

    double precision function ctthermo_meanMolecularWeight(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctthermo_meanmolecularweight = phase_meanmolecularweight(self%thermo_id)
    end function ctthermo_meanMolecularWeight

    integer function ctthermo_elementIndex(self, nm)
      implicit none
      type(phase_t), intent(inout) :: self
      character*(*), intent(in) :: nm
      ctthermo_elementindex = phase_elementindex(self%thermo_id, nm)
    end function ctthermo_elementIndex

    integer function ctthermo_speciesIndex(self, nm)
      implicit none
      type(phase_t), intent(inout) :: self
      character*(*), intent(in) :: nm
      ctthermo_speciesindex = phase_speciesindex(self%thermo_id, nm)
    end function ctthermo_speciesIndex

    subroutine ctthermo_getMoleFractions(self, x)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(out) :: x(self%nsp)
      self%err = phase_getmolefractions(self%thermo_id, x)
    end subroutine ctthermo_getMoleFractions

    double precision function ctthermo_moleFraction(self, k)
      implicit none
      type(phase_t), intent(inout) :: self
      integer, intent(in) :: k
      ctthermo_molefraction = phase_molefraction(self%thermo_id, k)
    end function ctthermo_moleFraction

    subroutine ctthermo_getMassFractions(self, y)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(out) :: y(self%nsp)
      self%err = phase_getmassfractions(self%thermo_id, y)
    end subroutine ctthermo_getMassFractions

    double precision function ctthermo_massFraction(self, k)
      implicit none
      type(phase_t), intent(inout) :: self
      integer, intent(in) :: k
      ctthermo_massfraction = phase_massfraction(self%thermo_id, k)
    end function ctthermo_massFraction

    subroutine ctthermo_setMoleFractions(self, x, norm)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(in) :: x(self%nsp)
      integer, intent(in), optional :: norm
      integer :: n
      if (present(norm)) then
         n = norm
      else
         n = 1
      end if
      self%err = phase_setmolefractions(self%thermo_id, x, n)
    end subroutine ctthermo_setmolefractions

    subroutine ctthermo_setMoleFractionsByName(self, x)
      implicit none
      type(phase_t), intent(inout) :: self
      character*(*), intent(in) :: x
      self%err = phase_setmolefractionsbyname(self%thermo_id, x)
    end subroutine ctthermo_setmolefractionsbyname

    subroutine ctthermo_setMassFractions(self, y, norm)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(in) :: y(self%nsp)
      integer, intent(in), optional :: norm
      integer :: n
      if (present(norm)) then
         n = norm
      else
         n = 1
      end if
      self%err = phase_setmassfractions(self%thermo_id, y, n)
    end subroutine ctthermo_setmassfractions

    subroutine ctthermo_setMassFractionsByName(self, y)
      implicit none
      type(phase_t), intent(inout) :: self
      character*(*), intent(in) :: y
      self%err = phase_setmassfractionsbyname(self%thermo_id, y)
    end subroutine ctthermo_setmassfractionsbyname

    subroutine ctthermo_getAtomicWeights(self, atw)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(out) :: atw(self%nel)
      self%err = phase_getatomicweights(self%thermo_id, atw)
    end subroutine ctthermo_getatomicweights

    subroutine ctthermo_getMolecularWeights(self, mw)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(out) :: mw(self%nsp)
      self%err = phase_getmolecularweights(self%thermo_id, mw)
    end subroutine ctthermo_getmolecularweights

    subroutine ctthermo_getSpeciesName(self, k, nm)
      implicit none
      type(phase_t), intent(inout) :: self
      integer, intent(in) :: k
      character*(*), intent(out) :: nm
      self%err = phase_getspeciesname(self%thermo_id, k, nm)
    end subroutine ctthermo_getspeciesname

    subroutine ctthermo_getElementName(self, m, nm)
      implicit none
      type(phase_t), intent(inout) :: self
      integer, intent(in) :: m
      character*(*), intent(out) :: nm
      self%err = phase_getelementname(self%thermo_id, m, nm)
    end subroutine ctthermo_getelementname

    double precision function ctthermo_nAtoms(self, k, m)
      implicit none
      type(phase_t), intent(inout) :: self
      integer, intent(in) :: k
      integer, intent(in) :: m
      ctthermo_natoms = phase_natoms(self%thermo_id, k, m)
    end function ctthermo_natoms

    subroutine ctthermo_getEosType(self, nm)
      implicit none
      type(phase_t), intent(inout) :: self
      character*(*), intent(out) :: nm
      self%err = th_geteostype(self%thermo_id, nm)
    end subroutine ctthermo_getEosType

    double precision function ctthermo_enthalpy_mole(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctthermo_enthalpy_mole = th_enthalpy_mole(self%thermo_id)
    end function ctthermo_enthalpy_mole

    double precision function ctthermo_intEnergy_mole(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctthermo_intenergy_mole = th_intenergy_mole(self%thermo_id)
    end function ctthermo_intenergy_mole

    double precision function ctthermo_entropy_mole(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctthermo_entropy_mole = th_entropy_mole(self%thermo_id)
    end function ctthermo_entropy_mole

    double precision function ctthermo_gibbs_mole(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctthermo_gibbs_mole = th_gibbs_mole(self%thermo_id)
    end function ctthermo_gibbs_mole

    double precision function ctthermo_cp_mole(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctthermo_cp_mole = th_cp_mole(self%thermo_id)
    end function ctthermo_cp_mole

    double precision function ctthermo_cv_mole(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctthermo_cv_mole = th_cv_mole(self%thermo_id)
    end function ctthermo_cv_mole

    double precision function ctthermo_pressure(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctthermo_pressure = th_pressure(self%thermo_id)
    end function ctthermo_pressure

    double precision function ctthermo_enthalpy_mass(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctthermo_enthalpy_mass = th_enthalpy_mass(self%thermo_id)
    end function ctthermo_enthalpy_mass

    double precision function ctthermo_intEnergy_mass(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctthermo_intEnergy_mass = th_intenergy_mass(self%thermo_id)
    end function ctthermo_intEnergy_mass

    double precision function ctthermo_entropy_mass(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctthermo_entropy_mass = th_entropy_mass(self%thermo_id)
    end function ctthermo_entropy_mass

    double precision function ctthermo_gibbs_mass(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctthermo_gibbs_mass = th_gibbs_mass(self%thermo_id)
    end function ctthermo_gibbs_mass

    double precision function ctthermo_cp_mass(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctthermo_cp_mass = th_cp_mass(self%thermo_id)
    end function ctthermo_cp_mass

    double precision function ctthermo_cv_mass(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctthermo_cv_mass = th_cv_mass(self%thermo_id)
    end function ctthermo_cv_mass

    subroutine ctthermo_chemPotentials(self, mu)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(out) :: mu(self%nsp)
      self%err = th_chempotentials(self%thermo_id, mu)
    end subroutine ctthermo_chempotentials

    subroutine ctthermo_setPressure(self, p)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(in) :: p
      self%err = th_setpressure(self%thermo_id, p)
    end subroutine ctthermo_setpressure

    subroutine ctthermo_setState_TPX(self, t, p, x)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      double precision, intent(in) :: x(*)
      call ctthermo_setTemperature(self, t)
      call ctthermo_setMoleFractions(self, x)
      call ctthermo_setPressure(self, p)
    end subroutine ctthermo_setState_TPX

    subroutine ctstring_setState_TPX(self, t, p, x)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      character*(*), intent(in) :: x
      call ctthermo_setTemperature(self, t)
      call ctthermo_setMoleFractionsByName(self, x)
      call ctthermo_setPressure(self, p)
    end subroutine ctstring_setState_TPX

    subroutine ctthermo_setState_TRX(self, t, rho, x)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: rho
      double precision, intent(in) :: x(*)
      call ctthermo_setTemperature(self, t)
      call ctthermo_setMoleFractions(self, x)
      call ctthermo_setDensity(self, rho)
    end subroutine ctthermo_setState_TRX

    subroutine ctstring_setState_TRX(self, t, rho, x)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: rho
      character*(*), intent(in) :: x
      call ctthermo_setTemperature(self, t)
      call ctthermo_setMoleFractionsByName(self, x)
      call ctthermo_setDensity(self, rho)
    end subroutine ctstring_setState_TRX

    subroutine ctthermo_setState_TRY(self, t, rho, y)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: rho
      double precision, intent(in) :: y(*)
      call ctthermo_setTemperature(self, t)
      call ctthermo_setMassFractions(self, y)
      call ctthermo_setDensity(self, rho)
    end subroutine ctthermo_setState_TRY

    subroutine ctstring_setState_TRY(self, t, rho, y)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: rho
      character*(*), intent(in) :: y
      call ctthermo_setTemperature(self, t)
      call ctthermo_setMassFractionsByName(self, y)
      call ctthermo_setDensity(self, rho)
    end subroutine ctstring_setState_TRY

    subroutine ctthermo_setState_TPY(self, t, p, y)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      double precision, intent(in) :: y(*)
      call ctthermo_setTemperature(self, t)
      call ctthermo_setMassFractions(self, y)
      call ctthermo_setPressure(self, p)
    end subroutine ctthermo_setState_TPY

    subroutine ctstring_setState_TPY(self, t, p, y)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      character*(*), intent(in) :: y
      call ctthermo_setTemperature(self, t)
      call ctthermo_setMassFractionsByName(self, y)
      call ctthermo_setPressure(self, p)
    end subroutine ctstring_setState_TPY


    subroutine ctthermo_setState_HP(self, h, p)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(in) :: h
      double precision, intent(in) :: p
      self%err = th_set_hp(self%thermo_id, h, p)
    end subroutine ctthermo_setstate_hp

    subroutine ctthermo_setState_UV(self, u, v)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(in) :: u
      double precision, intent(in) :: v
      self%err = th_set_uv(self%thermo_id, u, v)
    end subroutine ctthermo_setstate_uv

    subroutine ctthermo_setState_SV(self, s, v)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(in) :: s
      double precision, intent(in) :: v
      self%err = th_set_sv(self%thermo_id, s, v)
    end subroutine ctthermo_setstate_sv

    subroutine ctthermo_setState_SP(self, s, p)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(in) :: s
      double precision, intent(in) :: p
      self%err = th_set_sp(self%thermo_id, s, p)
    end subroutine ctthermo_setstate_sp

    subroutine ctthermo_equilibrate(self, XY, solver, rtol, max_steps, max_iter, estimate_equil, log_level)
      implicit none
      type(phase_t), intent(inout) :: self
      character*(*), intent(in) :: XY
      character*(*), intent(in), optional :: solver
      double precision, intent(in), optional :: rtol
      integer, intent(in), optional :: max_steps
      integer, intent(in), optional :: max_iter
      integer, intent(in), optional :: estimate_equil
      integer, intent(in), optional :: log_level
      character*(50) :: solver_
      double precision :: rtol_
      integer :: max_steps_
      integer :: max_iter_
      integer :: estimate_equil_
      integer :: log_level_
      solver_ = 'auto'
      rtol_ = 1e-9
      max_steps_ = 50000
      max_iter_ = 100
      estimate_equil_ = 0
      log_level_ = 0
      if(present(solver)) then
          solver_ = solver
      endif
      if(present(rtol)) then
          rtol_ = rtol
      endif
      if(present(max_steps)) then
          max_steps_ = max_steps
      endif
      if(present(max_iter)) then
          max_iter_ = max_iter
      endif
      if(present(estimate_equil)) then
          estimate_equil_ = estimate_equil
      endif
      if(present(log_level)) then
          log_level_ = log_level
      endif
      self%err = th_equil(self%thermo_id, XY, trim(solver_), rtol_, max_steps_, max_iter_, estimate_equil_, log_level_)
    end subroutine ctthermo_equilibrate

    double precision function ctthermo_refPressure(self)
      implicit none
      type(phase_t), intent(inout) :: self
      ctthermo_refpressure = th_refpressure(self%thermo_id)
    end function ctthermo_refpressure

    double precision function ctthermo_minTemp(self, k)
      implicit none
      type(phase_t), intent(inout) :: self
      integer, intent(in) :: k
      ctthermo_mintemp = th_mintemp(self%thermo_id, k)
    end function ctthermo_mintemp

    double precision function ctthermo_maxTemp(self, k)
      implicit none
      type(phase_t), intent(inout) :: self
      integer, intent(in) :: k
      ctthermo_maxtemp = th_maxtemp(self%thermo_id, k)
    end function ctthermo_maxtemp

    subroutine ctthermo_getEnthalpies_RT(self, h_rt)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(out) :: h_rt(self%nsp)
      self%err = th_getenthalpies_rt(self%thermo_id, h_rt)
    end subroutine ctthermo_getenthalpies_rt

    subroutine ctthermo_getEntropies_R(self, s_r)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(out) :: s_r(self%nsp)
      self%err = th_getentropies_r(self%thermo_id, s_r)
    end subroutine ctthermo_getentropies_r

    subroutine ctthermo_getCp_R(self, lenm, cp_r)
      implicit none
      type(phase_t), intent(inout) :: self
      integer, intent(out) :: lenm
      double precision, intent(out) :: cp_r(self%nsp)
      self%err = th_getcp_r(self%thermo_id, lenm, cp_r)
    end subroutine ctthermo_getcp_r

    subroutine ctthermo_getPartialMolarIntEnerg_R(self, ie)
      implicit none
      type(phase_t), intent(inout) :: self
      double precision, intent(out) :: ie(self%nsp)
      self%err = th_getpartialmolarintenergies_r(self%thermo_id, ie)
    end subroutine ctthermo_getPartialMolarIntEnerg_R

end module cantera_thermo
