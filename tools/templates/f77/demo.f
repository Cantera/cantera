c
c     Replace this sample main program with your program
c
c     This program uses functions defined in demo_ftnlib.cpp.
c
      program demo
      implicit double precision (a-h,o-z)
      parameter (MAXSP = 20, MAXRXNS = 100)
      double precision q(MAXRXNS), qf(MAXRXNS), qr(MAXRXNS)
      double precision x(MAXSP), y(MAXSP), wdot(MAXSP)
      character*80 eq
c
      call newIdealGasMix('h2o2.xml','')
      t = 1200.0
      p = 101325.0
      call setState_TPX_String(t, p,
     $     'H2:2, O2:1, OH:0.01, H:0.01, O:0.01')
      call equilibrate('HP')

      write(*,*) '**** Fortran 77 Test Program ****'


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
      call getFwdRatesOfProgress(qf)
      call getRevRatesOfProgress(qr)
      call getNetRatesOfProgress(q)
      do i = 1,irxns
         call getReactionEqn(i,eq)
         write(*,20) eq,qf(i),qr(i),q(i)
 20      format(a20,3g14.5,' kmol/m3/s')
      end do
      stop
      end

      
