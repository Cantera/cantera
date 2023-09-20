! Fortran 90 demo
! ===============
!
! This program illustrates using Cantera in Fortran 90 to compute
! thermodynamic, kinetic, and transport properties of a gas mixture.
!
! .. tags:: Fortran 90, tutorial, equilibrium, thermodynamics, kinetics, transport

! This file is part of Cantera. See License.txt in the top-level directory or
! at https://cantera.org/license.txt for license and copyright information.

program main

  ! use the Cantera module
  use cantera

  implicit none

  ! objects representing phases of matter have type 'phase_t'
  type(phase_t) :: gas, surf
  type(interface_t) iface
  integer nsp, nrxns
  double precision :: t, p
  character*1000 err_buf

  write(*,*)
  write(*,*) '********   Fortran 90 Test Program   ********'


  ! Read in a definition of the 'gas' phase.
  ! This will take the definition with name 'ohmech' from file
  ! 'h2o2.yaml', located in the Cantera data directory
  gas = importPhase('h2o2.yaml','ohmech')

  t = 1200.0     ! K
  p = 101325.0   ! Pa

  ! set the temperature, pressure, and mole fractions.
  call setState_TPX(gas, t, p, 'H2:1, O2:1, AR:2')


  nsp = nSpecies(gas)      ! number of species
  nrxns = nReactions(gas)  ! number of reactions

  call demo(gas, nsp, nrxns)

  ! Demo of importing a phase with surface kinetics
  gas = importPhase('ptcombust.yaml', 'gas')
  iface = importInterface('ptcombust.yaml', 'Pt_surf', gas)
  surf = iface%surf
  call setState_TPX(gas, t, p, 'H2:0.1, O2:0.65, H2O:0.2, CO:0.05')
  call setState_TPX(surf, t, p, 'O(S):0.01, PT(S): 0.8, CO(S): 0.19')
  call demo_surf(surf, nSpecies(surf), nReactions(surf))

end program main


!--------------------------------------------------------

subroutine demo(gas, MAXSP, MAXRXNS)

  ! use the Cantera module
  use cantera

  implicit none

  ! declare the arguments
  type(phase_t), intent(inout) :: gas
  integer, intent(in) :: MAXSP
  integer, intent(in) :: MAXRXNS

  double precision q(MAXRXNS), qf(MAXRXNS), qr(MAXRXNS)
  double precision diff(MAXSP)
  double precision molar_cp(MAXSP), molar_h(MAXSP), g_rt(MAXSP)

  character*80 eq
  character*20 name
  double precision :: dnu, dlam
  integer :: i, irxns, nsp, k

  write(*,*) 'Initial state properties:'
  write(*,10) temperature(gas), pressure(gas), density(gas), &
       enthalpy_mole(gas), entropy_mole(gas), cp_mole(gas)

  ! compute the equilibrium state holding the specific
  ! enthalpy and pressure constant
  call equilibrate(gas, 'HP')

  write(*,*) 'Equilibrium state properties:'
  write(*,10) temperature(gas), pressure(gas), density(gas), &
       enthalpy_mole(gas), entropy_mole(gas), cp_mole(gas)

10 format(//'Temperature:   ',g14.5,' K'/ &
       'Pressure:      ',g14.5,' Pa'/ &
       'Density:       ',g14.5,' kg/m3'/ &
       'Molar Enthalpy:',g14.5,' J/kmol'/ &
       'Molar Entropy: ',g14.5,' J/kmol-K'/ &
       'Molar cp:      ',g14.5,' J/kmol-K'//)


  !     Reaction information
  irxns = nReactions(gas)

  ! forward and reverse rates of progress should be equal
  ! in equilibrium states
  call getFwdRatesOfProgress(gas, qf)
  call getRevRatesOfProgress(gas, qr)

  ! net rates of progress should be zero in equilibrium states
  call getNetRatesOfProgress(gas, q)

  ! for each reaction, print the equation and the rates of progress
  do i = 1,irxns
     call getReactionString(gas, i,eq)
     write(*,20) eq,qf(i),qr(i),q(i)
20   format(a27,3e14.5,' kmol/m3/s')
  end do

  ! transport properties
  dnu = viscosity(gas)
  dlam = thermalConductivity(gas)
  call getMixDiffCoeffs(gas, diff)

  write(*,30) dnu, dlam
30 format(//'Viscosity:             ',g14.5,'  Pa-s'/ &
       'Thermal conductivity:  ',g14.5,'  W/m/K'/)

  write(*,*) 'Species Diffusion Coefficient'
  nsp = nSpecies(gas)
  do k = 1, nsp
     call getSpeciesName(gas, k, name)
     write(*,40) name, diff(k)
40   format(' ',a20,e14.5,' m2/s')
  end do

  ! Thermodynamic properties
  CALL getPartialMolarCp(gas, molar_cp)
  CALL getPartialMolarEnthalpies(gas, molar_h)
  CALL getGibbs_RT(gas, g_rt)

  write(*,*) 'Species molar cp, molar enthalpy, normalized Gibbs free energy'
  do k = 1, nsp
     call getSpeciesName(gas, k, name)
     write(*,50) name, molar_cp(k), molar_h(k), g_rt(k)
50   format(' ',a20,g14.5,' J/kmol-K',g14.5,' J/kmol',g14.5,'-')
  end do

  return

end subroutine demo

subroutine demo_surf(surf, nsp, nrxn)
  use cantera
  implicit none

  type(phase_t) surf
  integer :: nsp, nrxn, i
  character*40 equation
  double precision :: ropf(nrxn), ropnet(nrxn)

  write(*,*)
  write(*,*) '********   Interface Kinetics Test   ********'
  write(*,*)

  call getFwdRatesOfProgress(surf, ropf)
  call getNetRatesOfProgress(surf, ropnet)

  write(*,*) 'Equation                                    Fwd rate      Net rate'
  do i = 1, nrxn
    call getReactionString(surf, i, equation)
    write(*,60) equation, ropf(i), ropnet(i)
60  format(' ',a40, e14.5, e14.5)
  end do
end subroutine demo_surf
