c
c     Replace this sample main program with your program
c
c     This program uses functions defined in demo_ftnlib.cpp.
c
      program demo
      implicit double precision (a-h,o-z)
      parameter (MAXSP = 20, MAXRXNS = 100)
      double precision x(MAXSP), y(MAXSP), wdot(MAXRXNS)
c
      call readmechanism('h2o2.xml','')
      t = 1200.0
      p = 101325.0
      call setState_TPX_String(t, p, 'H2:2, O2:1')
      write(*,*) ' **** Test Program ****'
      write(*,10) temperature(), pressure(), density(),
     $     enthalpy_mole(), entropy_mole(), cp_mole()
 10   format(//'Temperature:',g14.5,' K'/'Pressure:',g14.5,' Pa'
     $     /'Density:',g14.5,' kg/m**3'//)
c
      stop
      end

      
