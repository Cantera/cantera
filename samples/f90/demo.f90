! This program illustrates using Cantera in Fortran 90 to compute
! thermodynamic, kinetic, and transport properties of a gas mixture.
!
program main

  ! use the Cantera module
  use cantera

  implicit none

  ! objects representing phases of matter have type 'phase_t'
  type(phase_t) gas
  integer nsp, nrxns
  double precision :: t, p

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

  return

end subroutine demo
