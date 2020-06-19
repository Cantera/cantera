c
c     Mach number vs. area for an isentropic flow.  See also the Python
c     version of this problem in the Python demos.
c
      program isentropic
      implicit double precision (a-h,o-z)
      parameter (oneatm = 1.01325d5, NPTS = 200)
      double precision a(NPTS), dmach(NPTS), t(NPTS),
     $     ratio(NPTS)

      call newIdealGasMix('gri30.yaml','gri30','')
      temp = 1200.d0
      pres = 10.d0*oneatm
      call setState_TPX_String(temp, pres,'H2:1, N2:0.1')

c     stagnation state properties
      s0 = entropy_mass()
      h0 = enthalpy_mass()
      p0 = pressure()

      dmdot = 1.0d0
      amin = 1.0d14

      do n = 1, NPTS
         p = p0*n/(NPTS+1)
         call setState_SP(s0,p)
         h = enthalpy_mass()
         rho = density()

         v2 = 2.0*(h0 - h)
         v = sqrt(v2)
         area = dmdot/(rho*v)
         if (area .lt. amin) then
            amin = area
         end if

         a(n) = area
         dmach(n) = v/soundspeed()
         t(n) = temperature()
         ratio(n) = p/p0
      end do

      do n = 1, NPTS
         a(n) = a(n)/amin
         write(*,30) a(n), dmach(n), t(n), ratio(n)
 30      format(4e16.5)
      end do
      end


      double precision function soundspeed()
      implicit double precision (a-h,o-z)
      double precision meanMolarMass
      parameter (R = 8314.4621d0)
      gamma = cp_mass()/cv_mass()
      soundspeed = sqrt(gamma * R * temperature()
     $     / meanMolarMass())
      return
      end
