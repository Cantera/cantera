!
!     Replace this sample main program with your program
!
!     This program uses functions defined in demo_ftnlib.cpp to create
!     an ideal gas mixture and print some its properties. 
!
!     For a C++ version of this program, see ../cxx/demo.cpp.
!
program main

  call demo(100, 500)

stop
end

subroutine demo(maxsp, maxrxns)
  use cantera
  implicit none 
  integer, intent(in) :: maxsp
  integer, intent(in) :: maxrxns

  type(phase_t) gas
  double precision q(MAXRXNS), qf(MAXRXNS), qr(MAXRXNS)
  double precision diff(MAXSP)
  character*80 eq
  character*20 name
  double precision :: t, p, dnu, dlam
  integer :: i, irxns, nsp, k

write(*,*)
write(*,*) '********   Fortran 90 Test Program   ********'

gas = importPhase('h2o2.cti','ohmech')

t = 1200.0
p = 101325.0

call setState_TPX(gas, t, p, 'H2:1, O2:1, AR:2')

write(*,*) 'Initial state properties:'
write(*,10) temperature(gas), pressure(gas), density(gas), &
     enthalpy_mole(gas), entropy_mole(gas), cp_mole(gas)
     
! compute the equilibrium state holding the specifi! 
! enthalpy and pressure constant
call equilibrate(gas, HP)

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

!     forward and reverse rates of progress should be equal
!     in equilibrium states
call getFwdRatesOfProgress(gas, qf)
call getRevRatesOfProgress(gas, qr)

!     net rates of progress should be zero in equilibrium states
call getNetRatesOfProgress(gas, q)

!     for each reaction, print the equation and the rates of progress
do i = 1,irxns
   call getReactionString(gas, i,eq)
   write(*,20) eq,qf(i),qr(i),q(i)
20 format(a27,3e14.5,' kmol/m3/s')
end do

!     Transport properties
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
40 format(' ',a20,e14.5,' m2/s')
end do

stop
end subroutine demo
