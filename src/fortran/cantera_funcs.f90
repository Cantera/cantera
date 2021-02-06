! This file is part of Cantera. See License.txt in the top-level directory or
! at https://cantera.org/license.txt for license and copyright information.

module cantera_funcs

  use cantera_xml
  use cantera_thermo
  use cantera_kinetics
  use cantera_transport
  use cantera_iface

  contains

    type(phase_t) function ctfunc_importPhase(src, id, loglevel)
      implicit none
      character*(*), intent(in) :: src
      character*(*), intent(in), optional :: id
      integer, intent(in), optional :: loglevel

      character(20) :: model
      type(phase_t) self

      if (present(id)) then
         self = ctthermo_newFromFile(src, id)
      else
         self = ctthermo_newFromFile(src)
      end if

      call ctkin_newFromFile(self, src, id)

      if (present(loglevel)) then
         self%tran_id = trans_newdefault(self%thermo_id, loglevel)
      else
         self%tran_id = trans_newdefault(self%thermo_id, 0)
      end if

      ctfunc_importPhase = self
      return
    end function ctfunc_importphase


    type(interface_t) function ctfunc_importInterface(src, id, gas, bulk, loglevel)
      implicit none
      character*(*), intent(in) :: src
      character*(*), intent(in), optional :: id
      type(phase_t) gas, bulk, surf
      integer, intent(in), optional :: loglevel

      character(20) :: model
      type(XML_Node) root, s, str
      type(interface_t) self

      root = new_XML_Node(src = src)
      if (present(id)) then
         s = ctxml_child(root, id = id)
      else
         s = ctxml_child(root, 'phase')
      end if

      surf = newThermoPhase(s)
      !write(*,*) 'from importInterface: nSpecies = ',ctthermo_nSpecies(surf)
      self = newInterface(s,surf,gas,bulk)
      ! call newKinetics(s, self)

      ctfunc_importInterface = self
      return
    end function ctfunc_importInterface


    subroutine ctfunc_phase_report(self, buf, show_thermo)
      implicit none
      type(phase_t), intent(inout) :: self
      character*(*), intent(out) :: buf
      integer, intent(in), optional :: show_thermo
      if (present(show_thermo)) then
         self%err = ctphase_report(self%thermo_id, buf, show_thermo)
      else
         self%err = ctphase_report(self%thermo_id, buf, 0)
      end if
    end subroutine ctfunc_phase_report

    subroutine ctfunc_getCanteraError(buf)
      implicit none
      integer :: ierr
      character*(*), intent(out) :: buf
      ierr = ctgetCanteraError(buf)
    end subroutine ctfunc_getCanteraError

    subroutine ctfunc_addCanteraDirectory(self, buf)
      implicit none
      type(phase_t), intent(inout) :: self
      character*(*), intent(in) :: buf
      self%err = ctaddCanteraDirectory(self%thermo_id, buf)
    end subroutine ctfunc_addCanteraDirectory

end module cantera_funcs
