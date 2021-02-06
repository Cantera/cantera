! This file is part of Cantera. See License.txt in the top-level directory or
! at https://cantera.org/license.txt for license and copyright information.

module cantera_iface

  use fct
  use cantera_xml
  use cantera_thermo
  use cantera_kinetics

  type interface_t
    type(phase_t) :: gas
    type(phase_t) :: bulk
    type(phase_t) :: surf
  end type interface_t

contains

    type(interface_t) function newInterface(xml_phase, surf, gas, bulk)
      implicit none
      type(XML_Node), intent(inout):: xml_phase
      type(phase_t) :: gas
      type(phase_t) :: bulk
      type(phase_t) :: surf
      type(interface_t) :: iface
      iface%gas = gas
      iface%bulk = bulk
      iface%surf = surf
      call newKinetics(xml_phase, surf, gas, bulk)
      newInterface = iface
    end function newInterface

end module cantera_iface

