c
c     This example shows how to implement subroutines that emulate those
c     of the Chemkin CKLIB library. This may be useful to port an
c     existing Chemkin-based application to Cantera. As shown here, the
c     subroutine names begin with 'ct' instead of 'ck', so that Cantera
c     and CKLIB subroutines can be both used in an application, if
c     desired. It is also possible to rename these subroutines with the
c     'ck' prefix if the application is not linked to the Chemkin CKLIB
c     library. In this case, application programs do not need to be
c     modified or recompiled - they only need to be relinked.
c
c     Only a few subroutines are implemented here, but the same idea can
c     be applied to create Cantera-based versions of any other
c     subroutines in the CKLIB library.
c

c-----------------------------------------------------------------------
c     example driver program

      program ctck
      implicit double precision (a-h,o-z)
      double precision rckwrk(1)
      integer ickwrk(1)
      parameter (MAXSP = 100)
      double precision y(MAXSP), wdot(MAXSP)

c     Read in the reaction mechanism. Since this is done differently
c     than in Chemkin, this function does not correspond to any CKLIB
c     subroutine.
      call newIdealGasMix('gri30.yaml','gri30','')

c     get the number of elements, species, and reactions
      call ctindx(ickwrk, rckwrk, mm, kk, ii)

      do k = 1, kk
         y(k) = 1.0/kk
      end do

c     compute the net production rates in cgs units
      p = 1.0d6
      t = 2500.0d0

      call ctwyp(p, t, y, ickwrk, rckwrk, wdot)
      do k = 1, kk
         write(*,*) k, y(k), wdot(k)
      end do

      end

c----------------------------------------------------------------------
c
c     The subroutines below emulate ones in the Chemkin CKLIB
c     library. They are implemented in terms of the procedures in
c     demo_ftnlib. It would also be possible to rewrite demo_ftnlib to
c     implement a Chemkin-like interface directly. Note that the arrays
c     ickwrk and rckwrk are passed in for consistency with the Chemkin
c     interface specification, but the are not used. These may simply be
c     dummy arrays, as in the main program above.
c

c     CTINDX: get the number of elements, species, and reactions

      subroutine ctindx(ickwrk, rckwrk, mm, kk, ii)
      implicit double precision (a-h,o-z)
      integer ickwrk(*)
      double precision rckwrk(*)
      mm = nElements()
      kk = nSpecies()
      ii = nReactions()
      return
      end


c     CTWYP: get the net molar production rates, given the pressure,
c     temperature, and array of mass fractions.

      subroutine ctwyp(p,t,y,ickwrk,rckwrk,wdot)
      implicit double precision (a-h,o-z)
      double precision y(*), rckwrk(*), wdot(*)
      integer ickwrk(*)

c     set the state
      psi = 0.1*p
      call setState_TPY(t, psi, y)

c     get the net production rates
      call getNetProductionRates(wdot)

c     convert SI -> cgs
      nsp = nSpecies()
      do k = 1, nsp
         wdot(k) = 1.0d3*wdot(k)
      end do
      return
      end
