! ------------------------------------------------------------------------
! 
! Copyright 2002 California Institute of Technology
! 
! Version 1.2.0.1
! Mon Jan  7 07:13:18 2002
! This file generated automatically.
! Fortran 90 module implementing Cantera class Reactor
! ------------------------------------------------------------------------
      module ctreactor
      use cttypes
      interface

      integer function creac_newreactor()
!DEC$ attributes c, reference :: creac_newreactor
!DEC$ attributes alias:'_creac_newreactor_' :: creac_newreactor
!DEC$ attributes dllimport :: creac_newreactor
      end function creac_newreactor

      integer function creac_newreservoir()
!DEC$ attributes c, reference :: creac_newreservoir
!DEC$ attributes alias:'_creac_newreservoir_' :: creac_newreservoir
!DEC$ attributes dllimport :: creac_newreservoir
      end function creac_newreservoir
      type(Reactor) function creac_StirredReactor()
!DEC$ attributes c, reference ::  creac_StirredReactor
!DEC$ attributes alias: '_creac_stirredreactor_' ::  creac_StirredReactor
!DEC$ attributes dllimport ::  creac_StirredReactor
      end function

      type(Reactor) function creac_Reservoir()
!DEC$ attributes c, reference ::  creac_Reservoir
!DEC$ attributes alias: '_creac_reservoir_' ::  creac_Reservoir
!DEC$ attributes dllimport ::  creac_Reservoir
      end function

      subroutine creac_setivol(reac, vol)
!DEC$ attributes c, reference ::  creac_setivol
!DEC$ attributes alias: '_creac_setivol_' ::  creac_setivol
!DEC$ attributes dllimport ::  creac_setivol
      type(Reactor), intent(inout) :: reac
      double precision, intent(in) :: vol
      end subroutine

      subroutine creac_setitm(reac, time)
!DEC$ attributes c, reference ::  creac_setitm
!DEC$ attributes alias: '_creac_setitm_' ::  creac_setitm
!DEC$ attributes dllimport ::  creac_setitm
      type(Reactor), intent(inout) :: reac
      double precision, intent(in) :: time
      end subroutine

      subroutine creac_setMaxStep(reac, maxstep)
!DEC$ attributes c, reference ::  creac_setMaxStep
!DEC$ attributes alias: '_creac_setmaxstep_' ::  creac_setMaxStep
!DEC$ attributes dllimport ::  creac_setMaxStep
      type(Reactor), intent(inout) :: reac
      double precision, intent(in) :: maxstep
      end subroutine

      subroutine creac_setArea(reac, area)
!DEC$ attributes c, reference ::  creac_setArea
!DEC$ attributes alias: '_creac_setarea_' ::  creac_setArea
!DEC$ attributes dllimport ::  creac_setArea
      type(Reactor), intent(inout) :: reac
      double precision, intent(in) :: area
      end subroutine

      subroutine creac_setExtTemp(reac, ts)
!DEC$ attributes c, reference ::  creac_setExtTemp
!DEC$ attributes alias: '_creac_setexttemp_' ::  creac_setExtTemp
!DEC$ attributes dllimport ::  creac_setExtTemp
      type(Reactor), intent(inout) :: reac
      double precision, intent(in) :: ts
      end subroutine

      subroutine creac_setExtRadTemp(reac, trad)
!DEC$ attributes c, reference ::  creac_setExtRadTemp
!DEC$ attributes alias: '_creac_setextradtemp_' ::  creac_setExtRadTemp
!DEC$ attributes dllimport ::  creac_setExtRadTemp
      type(Reactor), intent(inout) :: reac
      double precision, intent(in) :: trad
      end subroutine

      subroutine creac_seth(reac, h)
!DEC$ attributes c, reference ::  creac_seth
!DEC$ attributes alias: '_creac_seth_' ::  creac_seth
!DEC$ attributes dllimport ::  creac_seth
      type(Reactor), intent(inout) :: reac
      double precision, intent(in) :: h
      end subroutine

      subroutine creac_setVDotCoeff(reac, k)
!DEC$ attributes c, reference ::  creac_setVDotCoeff
!DEC$ attributes alias: '_creac_setvdotcoeff_' ::  creac_setVDotCoeff
!DEC$ attributes dllimport ::  creac_setVDotCoeff
      type(Reactor), intent(inout) :: reac
      double precision, intent(in) :: k
      end subroutine

      subroutine creac_setEmissivity(reac, emis)
!DEC$ attributes c, reference ::  creac_setEmissivity
!DEC$ attributes alias: '_creac_setemissivity_' ::  creac_setEmissivity
!DEC$ attributes dllimport ::  creac_setEmissivity
      type(Reactor), intent(inout) :: reac
      double precision, intent(in) :: emis
      end subroutine

      subroutine creac_setepr(reac, p0)
!DEC$ attributes c, reference ::  creac_setepr
!DEC$ attributes alias: '_creac_setepr_' ::  creac_setepr
!DEC$ attributes dllimport ::  creac_setepr
      type(Reactor), intent(inout) :: reac
      double precision, intent(in) :: p0
      end subroutine

      subroutine creac_insmix(reac, mix)
!DEC$ attributes c, reference ::  creac_insmix
!DEC$ attributes alias: '_creac_insmix_' ::  creac_insmix
!DEC$ attributes dllimport ::  creac_insmix
      type(Reactor), intent(inout) :: reac
      type(Mixture), intent(inout) :: mix
      end subroutine

      type(Mixture) function creac_contents(reac)
!DEC$ attributes c, reference ::  creac_contents
!DEC$ attributes alias: '_creac_contents_' ::  creac_contents
!DEC$ attributes dllimport ::  creac_contents
      type(Reactor), intent(inout) :: reac
      end function

      subroutine creac_advance(reac, time)
!DEC$ attributes c, reference ::  creac_advance
!DEC$ attributes alias: '_creac_advance_' ::  creac_advance
!DEC$ attributes dllimport ::  creac_advance
      type(Reactor), intent(inout) :: reac
      double precision, intent(in) :: time
      end subroutine

      double precision function creac_residenceTime(reac)
!DEC$ attributes c, reference ::  creac_residenceTime
!DEC$ attributes alias: '_creac_residencetime_' ::  creac_residenceTime
!DEC$ attributes dllimport ::  creac_residenceTime
      type(Reactor), intent(inout) :: reac
      end function

      double precision function creac_time(reac)
!DEC$ attributes c, reference ::  creac_time
!DEC$ attributes alias: '_creac_time_' ::  creac_time
!DEC$ attributes dllimport ::  creac_time
      type(Reactor), intent(inout) :: reac
      end function

      double precision function creac_volume(reac)
!DEC$ attributes c, reference ::  creac_volume
!DEC$ attributes alias: '_creac_volume_' ::  creac_volume
!DEC$ attributes dllimport ::  creac_volume
      type(Reactor), intent(inout) :: reac
      end function

      double precision function creac_density(reac)
!DEC$ attributes c, reference ::  creac_density
!DEC$ attributes alias: '_creac_density_' ::  creac_density
!DEC$ attributes dllimport ::  creac_density
      type(Reactor), intent(inout) :: reac
      end function

      double precision function creac_temperature(reac)
!DEC$ attributes c, reference ::  creac_temperature
!DEC$ attributes alias: '_creac_temperature_' ::  creac_temperature
!DEC$ attributes dllimport ::  creac_temperature
      type(Reactor), intent(inout) :: reac
      end function

      double precision function creac_enthalpy_mass(reac)
!DEC$ attributes c, reference ::  creac_enthalpy_mass
!DEC$ attributes alias: '_creac_enthalpy_mass_' ::  creac_enthalpy_mass
!DEC$ attributes dllimport ::  creac_enthalpy_mass
      type(Reactor), intent(inout) :: reac
      end function

      double precision function creac_umass(reac)
!DEC$ attributes c, reference ::  creac_umass
!DEC$ attributes alias: '_creac_umass_' ::  creac_umass
!DEC$ attributes dllimport ::  creac_umass
      type(Reactor), intent(inout) :: reac
      end function

      double precision function creac_pressure(reac)
!DEC$ attributes c, reference ::  creac_pressure
!DEC$ attributes alias: '_creac_pressure_' ::  creac_pressure
!DEC$ attributes dllimport ::  creac_pressure
      type(Reactor), intent(inout) :: reac
      end function

      double precision function creac_mass(reac)
!DEC$ attributes c, reference ::  creac_mass
!DEC$ attributes alias: '_creac_mass_' ::  creac_mass
!DEC$ attributes dllimport ::  creac_mass
      type(Reactor), intent(inout) :: reac
      end function

      subroutine creac_chemon(reac)
!DEC$ attributes c, reference ::  creac_chemon
!DEC$ attributes alias: '_creac_chemon_' ::  creac_chemon
!DEC$ attributes dllimport ::  creac_chemon
      type(Reactor), intent(inout) :: reac
      end subroutine

      subroutine creac_chemoff(reac)
!DEC$ attributes c, reference ::  creac_chemoff
!DEC$ attributes alias: '_creac_chemoff_' ::  creac_chemoff
!DEC$ attributes dllimport ::  creac_chemoff
      type(Reactor), intent(inout) :: reac
      end subroutine

      end interface
      contains

      type(Reactor) function StirredReactor() 
      type(Reactor) :: reac
      reac%hndl = creac_newreactor()
      StirredReactor = reac
      return
      end function

      type(Reactor) function Reservoir()
      type(Reactor) :: r
      r%hndl = creac_newreservoir()
      Reservoir = r
      return
      end function
 
      subroutine reac_copy(dest, src)
      type (Reactor), intent(out) :: dest
      type (Reactor), intent(in) :: src
      call mix_copy(dest%mix, src%mix)
      dest%hndl = src%hndl
      end subroutine
!
!     setInitialVolume
!
      subroutine reac_setivol(self, reac, vol)
      type(Reactor), intent(inout) :: self
      type(Reactor), intent(inout) :: reac
      double precision, intent(in) :: vol
      call creac_setivol(self%hndl, reac%hndl, vol)
      return
      end subroutine
!
!     setInitialTime
!
      subroutine reac_setitm(self, reac, time)
      type(Reactor), intent(inout) :: self
      type(Reactor), intent(inout) :: reac
      double precision, intent(in) :: time
      call creac_setitm(self%hndl, reac%hndl, time)
      return
      end subroutine
!
!     setMaxStep
!
      subroutine reac_setMaxStep(self, reac, maxstep)
      type(Reactor), intent(inout) :: self
      type(Reactor), intent(inout) :: reac
      double precision, intent(in) :: maxstep
      call creac_setMaxStep(self%hndl, reac%hndl, maxstep)
      return
      end subroutine
!
!     setArea
!
      subroutine reac_setArea(self, reac, area)
      type(Reactor), intent(inout) :: self
      type(Reactor), intent(inout) :: reac
      double precision, intent(in) :: area
      call creac_setArea(self%hndl, reac%hndl, area)
      return
      end subroutine
!
!     setExtTemp
!
      subroutine reac_setExtTemp(self, reac, ts)
      type(Reactor), intent(inout) :: self
      type(Reactor), intent(inout) :: reac
      double precision, intent(in) :: ts
      call creac_setExtTemp(self%hndl, reac%hndl, ts)
      return
      end subroutine
!
!     setExtRadTemp
!
      subroutine reac_setExtRadTemp(self, reac, trad)
      type(Reactor), intent(inout) :: self
      type(Reactor), intent(inout) :: reac
      double precision, intent(in) :: trad
      call creac_setExtRadTemp(self%hndl, reac%hndl, trad)
      return
      end subroutine
!
!     setHeatTransferCoeff
!
      subroutine reac_seth(self, reac, h)
      type(Reactor), intent(inout) :: self
      type(Reactor), intent(inout) :: reac
      double precision, intent(in) :: h
      call creac_seth(self%hndl, reac%hndl, h)
      return
      end subroutine
!
!     setVDotCoeff
!
      subroutine reac_setVDotCoeff(self, reac, k)
      type(Reactor), intent(inout) :: self
      type(Reactor), intent(inout) :: reac
      double precision, intent(in) :: k
      call creac_setVDotCoeff(self%hndl, reac%hndl, k)
      return
      end subroutine
!
!     setEmissivity
!
      subroutine reac_setEmissivity(self, reac, emis)
      type(Reactor), intent(inout) :: self
      type(Reactor), intent(inout) :: reac
      double precision, intent(in) :: emis
      call creac_setEmissivity(self%hndl, reac%hndl, emis)
      return
      end subroutine
!
!     setExtPressure
!
      subroutine reac_setepr(self, reac, p0)
      type(Reactor), intent(inout) :: self
      type(Reactor), intent(inout) :: reac
      double precision, intent(in) :: p0
      call creac_setepr(self%hndl, reac%hndl, p0)
      return
      end subroutine

      subroutine reac_insmix(reac, mix)
      type(Reactor), intent(inout) :: reac
      type(Mixture), intent(inout) :: mix
      call creac_insmix(reac%hndl, mix%hndl)
      reac%mix = mix
      return
      end subroutine

!
!     contents
!
      type(Mixture) function reac_contents(reac)
      type(Reactor), intent(inout) :: reac
      type(Mixture) mix
      mix%hndl = creac_contents(reac%hndl)
      reac_contents = mix
      return
      end function
!
!     advance
!
      subroutine reac_advance(self, reac, time)
      type(Reactor), intent(inout) :: self
      type(Reactor), intent(inout) :: reac
      double precision, intent(in) :: time
      call creac_advance(self%hndl, reac%hndl, time)
      return
      end subroutine
!
!     residenceTime
!
      double precision function reac_residenceTime(self, reac)
      type(Reactor), intent(inout) :: self
      type(Reactor), intent(inout) :: reac
      reac_residenceTime=creac_residenceTime(self%hndl, reac%hndl)
      return
      end function
!
!     time
!
      double precision function reac_time(self, reac)
      type(Reactor), intent(inout) :: self
      type(Reactor), intent(inout) :: reac
      reac_time=creac_time(self%hndl, reac%hndl)
      return
      end function
!
!     volume
!
      double precision function reac_volume(self, reac)
      type(Reactor), intent(inout) :: self
      type(Reactor), intent(inout) :: reac
      reac_volume=creac_volume(self%hndl, reac%hndl)
      return
      end function
!
!     density
!
      double precision function reac_density(self, reac)
      type(Reactor), intent(inout) :: self
      type(Reactor), intent(inout) :: reac
      reac_density=creac_density(self%hndl, reac%hndl)
      return
      end function
!
!     temperature
!
      double precision function reac_temperature(self, reac)
      type(Reactor), intent(inout) :: self
      type(Reactor), intent(inout) :: reac
      reac_temperature=creac_temperature(self%hndl, reac%hndl)
      return
      end function
!
!     enthalpy_mass
!
      double precision function reac_enthalpy_mass(self, reac)
      type(Reactor), intent(inout) :: self
      type(Reactor), intent(inout) :: reac
      reac_enthalpy_mass=creac_enthalpy_mass(self%hndl, reac%hndl)
      return
      end function
!
!     intEnergy_mass
!
      double precision function reac_umass(self, reac)
      type(Reactor), intent(inout) :: self
      type(Reactor), intent(inout) :: reac
      reac_umass=creac_umass(self%hndl, reac%hndl)
      return
      end function
!
!     pressure
!
      double precision function reac_pressure(self, reac)
      type(Reactor), intent(inout) :: self
      type(Reactor), intent(inout) :: reac
      reac_pressure=creac_pressure(self%hndl, reac%hndl)
      return
      end function
!
!     mass
!
      double precision function reac_mass(self, reac)
      type(Reactor), intent(inout) :: self
      type(Reactor), intent(inout) :: reac
      reac_mass=creac_mass(self%hndl, reac%hndl)
      return
      end function
!
!     enableChemistry
!
      subroutine reac_chemon(self, reac)
      type(Reactor), intent(inout) :: self
      type(Reactor), intent(inout) :: reac
      call creac_chemon(self%hndl, reac%hndl)
      return
      end subroutine
!
!     disableChemistry
!
      subroutine reac_chemoff(self, reac)
      type(Reactor), intent(inout) :: self
      type(Reactor), intent(inout) :: reac
      call creac_chemoff(self%hndl, reac%hndl)
      return
      end subroutine
      end module
