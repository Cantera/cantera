c
c     Replace this sample main program with your program
c
c     This program uses functions defined in demo_ftnlib.cpp to create
c     an ideal gas mixture and print some of its properties.
c
c     For a C++ version of this program, see ../cxx/demo.cpp.
c
      program demo
      implicit double precision (a-h,o-z)
      parameter (MAXSP = 20, MAXRXNS = 100)
      double precision q(MAXRXNS), qf(MAXRXNS), qr(MAXRXNS)
      double precision diff(MAXSP)
      character*80 eq
      character*20 name
c
      write(*,*)
      write(*,*) '********   Fortran 77 Test Program   ********'

      call newIdealGasMix('h2o2.yaml','ohmech','Mix')
      t = 1200.0
      p = 101325.0
      call setState_TPX_String(t, p,
     $     'H2:2, O2:1, OH:0.01, H:0.01, O:0.01')

c
      write(*,*) 'Initial state properties:'
      write(*,10) temperature(), pressure(), density(),
     $     enthalpy_mole(), entropy_mole(), cp_mole()

c     compute the equilibrium state holding the specific
c     enthalpy and pressure constant
      call equilibrate('HP')

      write(*,*) 'Equilibrium state properties:'
      write(*,10) temperature(), pressure(), density(),
     $     enthalpy_mole(), entropy_mole(), cp_mole()

 10   format(//'Temperature:   ',g14.5,' K'/
     $         'Pressure:      ',g14.5,' Pa'/
     $         'Density:       ',g14.5,' kg/m3'/
     $         'Molar Enthalpy:',g14.5,' J/kmol'/
     $         'Molar Entropy: ',g14.5,' J/kmol-K'/
     $         'Molar cp:      ',g14.5,' J/kmol-K'//)


c
c     Reaction information
c
      irxns = nReactions()

c     forward and reverse rates of progress should be equal
c     in equilibrium states
      call getFwdRatesOfProgress(qf)
      call getRevRatesOfProgress(qr)

c     net rates of progress should be zero in equilibrium states
      call getNetRatesOfProgress(q)

c     for each reaction, print the equation and the rates of progress
      do i = 1,irxns
         call getReactionEqn(i,eq)
         write(*,20) eq,qf(i),qr(i),q(i)
 20      format(a27,3e14.5,' kmol/m3/s')
      end do

c
c     Transport properties
c
      dnu = viscosity()
      dlam = thermalConductivity()
      call getMixDiffCoeffs(diff)

      write(*,30) dnu, dlam
 30   format(//'Viscosity:             ',g14.5,'  Pa-s'/
     $         'Thermal conductivity:  ',g14.5,'  W/m/K'/)
      write(*,*) 'Species            ',
     $     '    Diffusion Coefficient'
      nsp = nSpecies()
      do k = 1, nsp
         call getSpeciesName(k, name)
         write(*,40) name, diff(k)
 40      format(' ',a20,e14.5,' m2/s')
      end do

      stop
      end
