! ------------------------------------------------------------------------
! 
! Copyright 2002 California Institute of Technology
! 
! Version 1.2.0.1
! Mon Jan  7 07:13:18 2002
! This file generated automatically.
! Fortran 90 module implementing Cantera class Mixture
! ------------------------------------------------------------------------
      module ctmixture
      use cttypes
      interface

      integer function cmix_newmixture()
!DEC$ attributes alias:'_cmix_newmixture_' :: cmix_newmixture
!DEC$ attributes c, reference :: cmix_newmixture
      end function

      integer function cmix_newidealgas()
!DEC$ attributes alias:'_cmix_newidealgas_' :: cmix_newidealgas
!DEC$ attributes c, reference :: cmix_newidealgas
      end function cmix_newidealgas
      
      integer function cmix_import(i, inputfile, thermofile, ivalidate)
!DEC$ attributes c, reference, alias:'_cmix_import_' :: cmix_import
      character*(*), intent(in) :: inputfile, thermofile
      integer, intent(in) :: ivalidate
      end function cmix_import

      integer function cmix_newsub(i)
!DEC$ attributes alias:'_cmix_newsub_' :: cmix_newsub
!DEC$ attributes c, reference :: cmix_newsub
      integer, intent(in) :: i
      end function 
      type(Mixture) function cmix_CKGas(inputFile, thermoFile, iValidate)
!DEC$ attributes c, reference ::  cmix_CKGas
!DEC$ attributes alias: '_cmix_ckgas_' ::  cmix_CKGas
!DEC$ attributes dllimport ::  cmix_CKGas
      character*(*), intent(in) :: inputFile
      character*(*), intent(in) :: thermoFile
      integer, intent(in) :: iValidate
      end function

      type(Mixture) function cmix_GRIMech30()
!DEC$ attributes c, reference ::  cmix_GRIMech30
!DEC$ attributes alias: '_cmix_grimech30_' ::  cmix_GRIMech30
!DEC$ attributes dllimport ::  cmix_GRIMech30
      end function

      type(Mixture) function cmix_Water()
!DEC$ attributes c, reference ::  cmix_Water
!DEC$ attributes alias: '_cmix_water_' ::  cmix_Water
!DEC$ attributes dllimport ::  cmix_Water
      end function

      type(Mixture) function cmix_Nitrogen()
!DEC$ attributes c, reference ::  cmix_Nitrogen
!DEC$ attributes alias: '_cmix_nitrogen_' ::  cmix_Nitrogen
!DEC$ attributes dllimport ::  cmix_Nitrogen
      end function

      type(Mixture) function cmix_Methane()
!DEC$ attributes c, reference ::  cmix_Methane
!DEC$ attributes alias: '_cmix_methane_' ::  cmix_Methane
!DEC$ attributes dllimport ::  cmix_Methane
      end function

      type(Mixture) function cmix_Hydrogen()
!DEC$ attributes c, reference ::  cmix_Hydrogen
!DEC$ attributes alias: '_cmix_hydrogen_' ::  cmix_Hydrogen
!DEC$ attributes dllimport ::  cmix_Hydrogen
      end function

      type(Mixture) function cmix_Oxygen()
!DEC$ attributes c, reference ::  cmix_Oxygen
!DEC$ attributes alias: '_cmix_oxygen_' ::  cmix_Oxygen
!DEC$ attributes dllimport ::  cmix_Oxygen
      end function

      subroutine cmix_delete()
!DEC$ attributes c, reference ::  cmix_delete
!DEC$ attributes alias: '_cmix_delete_' ::  cmix_delete
!DEC$ attributes dllimport ::  cmix_delete
      end subroutine

      logical function cmix_ready()
!DEC$ attributes c, reference ::  cmix_ready
!DEC$ attributes alias: '_cmix_ready_' ::  cmix_ready
!DEC$ attributes dllimport ::  cmix_ready
      end function

      subroutine cmix_saveState(state)
!DEC$ attributes c, reference ::  cmix_saveState
!DEC$ attributes alias: '_cmix_savestate_' ::  cmix_saveState
!DEC$ attributes dllimport ::  cmix_saveState
      double precision, intent(out) :: state(*)
      end subroutine

      subroutine cmix_restoreState(state)
!DEC$ attributes c, reference ::  cmix_restoreState
!DEC$ attributes alias: '_cmix_restorestate_' ::  cmix_restoreState
!DEC$ attributes dllimport ::  cmix_restoreState
      double precision, intent(out) :: state(*)
      end subroutine

      integer function cmix_nElements()
!DEC$ attributes c, reference ::  cmix_nElements
!DEC$ attributes alias: '_cmix_nelements_' ::  cmix_nElements
!DEC$ attributes dllimport ::  cmix_nElements
      end function

      integer function cmix_nSpecies()
!DEC$ attributes c, reference ::  cmix_nSpecies
!DEC$ attributes alias: '_cmix_nspecies_' ::  cmix_nSpecies
!DEC$ attributes dllimport ::  cmix_nSpecies
      end function

      integer function cmix_nReactions()
!DEC$ attributes c, reference ::  cmix_nReactions
!DEC$ attributes alias: '_cmix_nreactions_' ::  cmix_nReactions
!DEC$ attributes dllimport ::  cmix_nReactions
      end function

      subroutine cmix_addElement(name, atomicWt)
!DEC$ attributes c, reference ::  cmix_addElement
!DEC$ attributes alias: '_cmix_addelement_' ::  cmix_addElement
!DEC$ attributes dllimport ::  cmix_addElement
      character*(*), intent(in) :: name
      double precision, intent(in) :: atomicWt
      end subroutine


      subroutine cmix_getename(i, m, ename)
!DEC$ attributes c, reference, alias:'_cmix_getename_' :: cmix_getename
      integer, intent(in) :: i, m
      character*(*), intent(out) :: ename
      end subroutine
      integer function cmix_eindex(name)
!DEC$ attributes c, reference ::  cmix_eindex
!DEC$ attributes alias: '_cmix_eindex_' ::  cmix_eindex
!DEC$ attributes dllimport ::  cmix_eindex
      character*(*), intent(in) :: name
      end function

      double precision function cmix_atwt(m)
!DEC$ attributes c, reference ::  cmix_atwt
!DEC$ attributes alias: '_cmix_atwt_' ::  cmix_atwt
!DEC$ attributes dllimport ::  cmix_atwt
      integer, intent(in) :: m
      end function

      double precision function cmix_nAtoms(k, m)
!DEC$ attributes c, reference ::  cmix_nAtoms
!DEC$ attributes alias: '_cmix_natoms_' ::  cmix_nAtoms
!DEC$ attributes dllimport ::  cmix_nAtoms
      integer, intent(in) :: k
      integer, intent(in) :: m
      end function

      double precision function cmix_charge(k)
!DEC$ attributes c, reference ::  cmix_charge
!DEC$ attributes alias: '_cmix_charge_' ::  cmix_charge
!DEC$ attributes dllimport ::  cmix_charge
      integer, intent(in) :: k
      end function

      integer function cmix_speciesIndex(name)
!DEC$ attributes c, reference ::  cmix_speciesIndex
!DEC$ attributes alias: '_cmix_speciesindex_' ::  cmix_speciesIndex
!DEC$ attributes dllimport ::  cmix_speciesIndex
      character*(*), intent(in) :: name
      end function


      subroutine cmix_getspnm(i, k, sname)
!DEC$ attributes c, reference, alias:'_cmix_getspnm_' :: cmix_getspnm
      integer, intent(in) :: i, k
      character*(*), intent(out) :: sname
      end subroutine
      double precision function cmix_molwt(k)
!DEC$ attributes c, reference ::  cmix_molwt
!DEC$ attributes alias: '_cmix_molwt_' ::  cmix_molwt
!DEC$ attributes dllimport ::  cmix_molwt
      integer, intent(in) :: k
      end function

      subroutine cmix_gmolwts(molWts)
!DEC$ attributes c, reference ::  cmix_gmolwts
!DEC$ attributes alias: '_cmix_gmolwts_' ::  cmix_gmolwts
!DEC$ attributes dllimport ::  cmix_gmolwts
      double precision, intent(out) :: molWts(*)
      end subroutine

      subroutine cmix_setx(moleFracs)
!DEC$ attributes c, reference ::  cmix_setx
!DEC$ attributes alias: '_cmix_setx_' ::  cmix_setx
!DEC$ attributes dllimport ::  cmix_setx
      double precision, intent(in) :: moleFracs(*)
      end subroutine

      subroutine cmix_sxnonorm(moleFracs)
!DEC$ attributes c, reference ::  cmix_sxnonorm
!DEC$ attributes alias: '_cmix_sxnonorm_' ::  cmix_sxnonorm
!DEC$ attributes dllimport ::  cmix_sxnonorm
      double precision, intent(in) :: moleFracs(*)
      end subroutine

      subroutine cmix_sety(massFracs)
!DEC$ attributes c, reference ::  cmix_sety
!DEC$ attributes alias: '_cmix_sety_' ::  cmix_sety
!DEC$ attributes dllimport ::  cmix_sety
      double precision, intent(in) :: massFracs(*)
      end subroutine

      subroutine cmix_synonorm(y)
!DEC$ attributes c, reference ::  cmix_synonorm
!DEC$ attributes alias: '_cmix_synonorm_' ::  cmix_synonorm
!DEC$ attributes dllimport ::  cmix_synonorm
      double precision, intent(in) :: y(*)
      end subroutine

      subroutine cmix_sconc(conc)
!DEC$ attributes c, reference ::  cmix_sconc
!DEC$ attributes alias: '_cmix_sconc_' ::  cmix_sconc
!DEC$ attributes dllimport ::  cmix_sconc
      double precision, intent(in) :: conc(*)
      end subroutine

      double precision function cmix_xbyname(name)
!DEC$ attributes c, reference ::  cmix_xbyname
!DEC$ attributes alias: '_cmix_xbyname_' ::  cmix_xbyname
!DEC$ attributes dllimport ::  cmix_xbyname
      character*(*), intent(in) :: name
      end function

      double precision function cmix_ybyname(name)
!DEC$ attributes c, reference ::  cmix_ybyname
!DEC$ attributes alias: '_cmix_ybyname_' ::  cmix_ybyname
!DEC$ attributes dllimport ::  cmix_ybyname
      character*(*), intent(in) :: name
      end function

      subroutine cmix_getx(moleFracs)
!DEC$ attributes c, reference ::  cmix_getx
!DEC$ attributes alias: '_cmix_getx_' ::  cmix_getx
!DEC$ attributes dllimport ::  cmix_getx
      double precision, intent(inout) :: moleFracs(*)
      end subroutine

      subroutine cmix_gety(massFracs)
!DEC$ attributes c, reference ::  cmix_gety
!DEC$ attributes alias: '_cmix_gety_' ::  cmix_gety
!DEC$ attributes dllimport ::  cmix_gety
      double precision, intent(inout) :: massFracs(*)
      end subroutine

      subroutine cmix_getconc(concentrations)
!DEC$ attributes c, reference ::  cmix_getconc
!DEC$ attributes alias: '_cmix_getconc_' ::  cmix_getconc
!DEC$ attributes dllimport ::  cmix_getconc
      double precision, intent(inout) :: concentrations(*)
      end subroutine

      double precision function cmix_mean_X(Q)
!DEC$ attributes c, reference ::  cmix_mean_X
!DEC$ attributes alias: '_cmix_mean_x_' ::  cmix_mean_X
!DEC$ attributes dllimport ::  cmix_mean_X
      double precision, intent(in) :: Q(*)
      end function

      double precision function cmix_mean_Y(Q)
!DEC$ attributes c, reference ::  cmix_mean_Y
!DEC$ attributes alias: '_cmix_mean_y_' ::  cmix_mean_Y
!DEC$ attributes dllimport ::  cmix_mean_Y
      double precision, intent(in) :: Q(*)
      end function

      double precision function cmix_meanmw()
!DEC$ attributes c, reference ::  cmix_meanmw
!DEC$ attributes alias: '_cmix_meanmw_' ::  cmix_meanmw
!DEC$ attributes dllimport ::  cmix_meanmw
      end function

      double precision function cmix_sum_xlogQ(Q)
!DEC$ attributes c, reference ::  cmix_sum_xlogQ
!DEC$ attributes alias: '_cmix_sum_xlogq_' ::  cmix_sum_xlogQ
!DEC$ attributes dllimport ::  cmix_sum_xlogQ
      double precision, intent(out) :: Q(*)
      end function

      subroutine cmix_setTPXarray(t, p, moleFracs)
!DEC$ attributes c, reference ::  cmix_setTPXarray
!DEC$ attributes alias: '_cmix_settpxarray_' ::  cmix_setTPXarray
!DEC$ attributes dllimport ::  cmix_setTPXarray
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      double precision, intent(in) :: moleFracs(*)
      end subroutine

      subroutine cmix_setTPXstr(t, p, moleFracs)
!DEC$ attributes c, reference ::  cmix_setTPXstr
!DEC$ attributes alias: '_cmix_settpxstr_' ::  cmix_setTPXstr
!DEC$ attributes dllimport ::  cmix_setTPXstr
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      character*(*), intent(in) :: moleFracs
      end subroutine

      subroutine cmix_setTPYarray(t, p, massFracs)
!DEC$ attributes c, reference ::  cmix_setTPYarray
!DEC$ attributes alias: '_cmix_settpyarray_' ::  cmix_setTPYarray
!DEC$ attributes dllimport ::  cmix_setTPYarray
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      double precision, intent(in) :: massFracs(*)
      end subroutine

      subroutine cmix_setTPYstr(t, p, massFracs)
!DEC$ attributes c, reference ::  cmix_setTPYstr
!DEC$ attributes alias: '_cmix_settpystr_' ::  cmix_setTPYstr
!DEC$ attributes dllimport ::  cmix_setTPYstr
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      character*(*), intent(in) :: massFracs
      end subroutine

      subroutine cmix_setState_TRX(t, density, moleFracs)
!DEC$ attributes c, reference ::  cmix_setState_TRX
!DEC$ attributes alias: '_cmix_setstate_trx_' ::  cmix_setState_TRX
!DEC$ attributes dllimport ::  cmix_setState_TRX
      double precision, intent(in) :: t
      double precision, intent(in) :: density
      double precision, intent(in) :: moleFracs(*)
      end subroutine

      subroutine cmix_setState_TRY(t, density, massFracs)
!DEC$ attributes c, reference ::  cmix_setState_TRY
!DEC$ attributes alias: '_cmix_setstate_try_' ::  cmix_setState_TRY
!DEC$ attributes dllimport ::  cmix_setState_TRY
      double precision, intent(in) :: t
      double precision, intent(in) :: density
      double precision, intent(in) :: massFracs(*)
      end subroutine

      subroutine cmix_settemp(temp)
!DEC$ attributes c, reference ::  cmix_settemp
!DEC$ attributes alias: '_cmix_settemp_' ::  cmix_settemp
!DEC$ attributes dllimport ::  cmix_settemp
      double precision, intent(in) :: temp
      end subroutine

      subroutine cmix_setDensity(density)
!DEC$ attributes c, reference ::  cmix_setDensity
!DEC$ attributes alias: '_cmix_setdensity_' ::  cmix_setDensity
!DEC$ attributes dllimport ::  cmix_setDensity
      double precision, intent(in) :: density
      end subroutine

      subroutine cmix_setState_TP(t, p)
!DEC$ attributes c, reference ::  cmix_setState_TP
!DEC$ attributes alias: '_cmix_setstate_tp_' ::  cmix_setState_TP
!DEC$ attributes dllimport ::  cmix_setState_TP
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      end subroutine

      subroutine cmix_setState_PX(p, x)
!DEC$ attributes c, reference ::  cmix_setState_PX
!DEC$ attributes alias: '_cmix_setstate_px_' ::  cmix_setState_PX
!DEC$ attributes dllimport ::  cmix_setState_PX
      double precision, intent(in) :: p
      double precision, intent(in) :: x(*)
      end subroutine

      subroutine cmix_setState_PY(p, massFracs)
!DEC$ attributes c, reference ::  cmix_setState_PY
!DEC$ attributes alias: '_cmix_setstate_py_' ::  cmix_setState_PY
!DEC$ attributes dllimport ::  cmix_setState_PY
      double precision, intent(in) :: p
      double precision, intent(in) :: massFracs(*)
      end subroutine

      subroutine cmix_setPressure(p)
!DEC$ attributes c, reference ::  cmix_setPressure
!DEC$ attributes alias: '_cmix_setpressure_' ::  cmix_setPressure
!DEC$ attributes dllimport ::  cmix_setPressure
      double precision, intent(in) :: p
      end subroutine

      subroutine cmix_setState_TR(t, rho)
!DEC$ attributes c, reference ::  cmix_setState_TR
!DEC$ attributes alias: '_cmix_setstate_tr_' ::  cmix_setState_TR
!DEC$ attributes dllimport ::  cmix_setState_TR
      double precision, intent(in) :: t
      double precision, intent(in) :: rho
      end subroutine

      subroutine cmix_setState_TX(t, moleFracs)
!DEC$ attributes c, reference ::  cmix_setState_TX
!DEC$ attributes alias: '_cmix_setstate_tx_' ::  cmix_setState_TX
!DEC$ attributes dllimport ::  cmix_setState_TX
      double precision, intent(in) :: t
      double precision, intent(inout) :: moleFracs(*)
      end subroutine

      subroutine cmix_setState_TY(t, massFracs)
!DEC$ attributes c, reference ::  cmix_setState_TY
!DEC$ attributes alias: '_cmix_setstate_ty_' ::  cmix_setState_TY
!DEC$ attributes dllimport ::  cmix_setState_TY
      double precision, intent(in) :: t
      double precision, intent(in) :: massFracs(*)
      end subroutine

      subroutine cmix_setState_RX(density, moleFracs)
!DEC$ attributes c, reference ::  cmix_setState_RX
!DEC$ attributes alias: '_cmix_setstate_rx_' ::  cmix_setState_RX
!DEC$ attributes dllimport ::  cmix_setState_RX
      double precision, intent(in) :: density
      double precision, intent(in) :: moleFracs(*)
      end subroutine

      subroutine cmix_setState_RY(density, y)
!DEC$ attributes c, reference ::  cmix_setState_RY
!DEC$ attributes alias: '_cmix_setstate_ry_' ::  cmix_setState_RY
!DEC$ attributes dllimport ::  cmix_setState_RY
      double precision, intent(in) :: density
      double precision, intent(in) :: y(*)
      end subroutine

      subroutine cmix_setState_HP(t, p)
!DEC$ attributes c, reference ::  cmix_setState_HP
!DEC$ attributes alias: '_cmix_setstate_hp_' ::  cmix_setState_HP
!DEC$ attributes dllimport ::  cmix_setState_HP
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      end subroutine

      subroutine cmix_setState_UV(t, p)
!DEC$ attributes c, reference ::  cmix_setState_UV
!DEC$ attributes alias: '_cmix_setstate_uv_' ::  cmix_setState_UV
!DEC$ attributes dllimport ::  cmix_setState_UV
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      end subroutine

      subroutine cmix_setState_SP(t, p)
!DEC$ attributes c, reference ::  cmix_setState_SP
!DEC$ attributes alias: '_cmix_setstate_sp_' ::  cmix_setState_SP
!DEC$ attributes dllimport ::  cmix_setState_SP
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      end subroutine

      subroutine cmix_setState_SV(t, p)
!DEC$ attributes c, reference ::  cmix_setState_SV
!DEC$ attributes alias: '_cmix_setstate_sv_' ::  cmix_setState_SV
!DEC$ attributes dllimport ::  cmix_setState_SV
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      end subroutine

      subroutine cmix_gh_rt(h_RT)
!DEC$ attributes c, reference ::  cmix_gh_rt
!DEC$ attributes alias: '_cmix_gh_rt_' ::  cmix_gh_rt
!DEC$ attributes dllimport ::  cmix_gh_rt
      double precision, intent(out) :: h_RT(*)
      end subroutine

      subroutine cmix_ggibbs_rt(g_RT)
!DEC$ attributes c, reference ::  cmix_ggibbs_rt
!DEC$ attributes alias: '_cmix_ggibbs_rt_' ::  cmix_ggibbs_rt
!DEC$ attributes dllimport ::  cmix_ggibbs_rt
      double precision, intent(out) :: g_RT(*)
      end subroutine

      subroutine cmix_gcp_r(cp_R)
!DEC$ attributes c, reference ::  cmix_gcp_r
!DEC$ attributes alias: '_cmix_gcp_r_' ::  cmix_gcp_r
!DEC$ attributes dllimport ::  cmix_gcp_r
      double precision, intent(out) :: cp_R(*)
      end subroutine

      subroutine cmix_gs_r(s_R)
!DEC$ attributes c, reference ::  cmix_gs_r
!DEC$ attributes alias: '_cmix_gs_r_' ::  cmix_gs_r
!DEC$ attributes dllimport ::  cmix_gs_r
      double precision, intent(out) :: s_R(*)
      end subroutine

      double precision function cmix_minTemp()
!DEC$ attributes c, reference ::  cmix_minTemp
!DEC$ attributes alias: '_cmix_mintemp_' ::  cmix_minTemp
!DEC$ attributes dllimport ::  cmix_minTemp
      end function

      double precision function cmix_maxTemp()
!DEC$ attributes c, reference ::  cmix_maxTemp
!DEC$ attributes alias: '_cmix_maxtemp_' ::  cmix_maxTemp
!DEC$ attributes dllimport ::  cmix_maxTemp
      end function

      double precision function cmix_refp()
!DEC$ attributes c, reference ::  cmix_refp
!DEC$ attributes alias: '_cmix_refp_' ::  cmix_refp
!DEC$ attributes dllimport ::  cmix_refp
      end function

      double precision function cmix_temperature()
!DEC$ attributes c, reference ::  cmix_temperature
!DEC$ attributes alias: '_cmix_temperature_' ::  cmix_temperature
!DEC$ attributes dllimport ::  cmix_temperature
      end function

      double precision function cmix_density()
!DEC$ attributes c, reference ::  cmix_density
!DEC$ attributes alias: '_cmix_density_' ::  cmix_density
!DEC$ attributes dllimport ::  cmix_density
      end function

      double precision function cmix_moldens()
!DEC$ attributes c, reference ::  cmix_moldens
!DEC$ attributes alias: '_cmix_moldens_' ::  cmix_moldens
!DEC$ attributes dllimport ::  cmix_moldens
      end function

      double precision function cmix_pressure()
!DEC$ attributes c, reference ::  cmix_pressure
!DEC$ attributes alias: '_cmix_pressure_' ::  cmix_pressure
!DEC$ attributes dllimport ::  cmix_pressure
      end function

      double precision function cmix_enthalpy_mole()
!DEC$ attributes c, reference ::  cmix_enthalpy_mole
!DEC$ attributes alias: '_cmix_enthalpy_mole_' ::  cmix_enthalpy_mole
!DEC$ attributes dllimport ::  cmix_enthalpy_mole
      end function

      double precision function cmix_umole()
!DEC$ attributes c, reference ::  cmix_umole
!DEC$ attributes alias: '_cmix_umole_' ::  cmix_umole
!DEC$ attributes dllimport ::  cmix_umole
      end function

      double precision function cmix_entropy_mole()
!DEC$ attributes c, reference ::  cmix_entropy_mole
!DEC$ attributes alias: '_cmix_entropy_mole_' ::  cmix_entropy_mole
!DEC$ attributes dllimport ::  cmix_entropy_mole
      end function

      double precision function cmix_gibbs_mole()
!DEC$ attributes c, reference ::  cmix_gibbs_mole
!DEC$ attributes alias: '_cmix_gibbs_mole_' ::  cmix_gibbs_mole
!DEC$ attributes dllimport ::  cmix_gibbs_mole
      end function

      double precision function cmix_cp_mole()
!DEC$ attributes c, reference ::  cmix_cp_mole
!DEC$ attributes alias: '_cmix_cp_mole_' ::  cmix_cp_mole
!DEC$ attributes dllimport ::  cmix_cp_mole
      end function

      double precision function cmix_cv_mole()
!DEC$ attributes c, reference ::  cmix_cv_mole
!DEC$ attributes alias: '_cmix_cv_mole_' ::  cmix_cv_mole
!DEC$ attributes dllimport ::  cmix_cv_mole
      end function

      subroutine cmix_gchempot(mu)
!DEC$ attributes c, reference ::  cmix_gchempot
!DEC$ attributes alias: '_cmix_gchempot_' ::  cmix_gchempot
!DEC$ attributes dllimport ::  cmix_gchempot
      double precision, intent(out) :: mu(*)
      end subroutine

      double precision function cmix_enthalpy_mass()
!DEC$ attributes c, reference ::  cmix_enthalpy_mass
!DEC$ attributes alias: '_cmix_enthalpy_mass_' ::  cmix_enthalpy_mass
!DEC$ attributes dllimport ::  cmix_enthalpy_mass
      end function

      double precision function cmix_umass()
!DEC$ attributes c, reference ::  cmix_umass
!DEC$ attributes alias: '_cmix_umass_' ::  cmix_umass
!DEC$ attributes dllimport ::  cmix_umass
      end function

      double precision function cmix_entropy_mass()
!DEC$ attributes c, reference ::  cmix_entropy_mass
!DEC$ attributes alias: '_cmix_entropy_mass_' ::  cmix_entropy_mass
!DEC$ attributes dllimport ::  cmix_entropy_mass
      end function

      double precision function cmix_gibbs_mass()
!DEC$ attributes c, reference ::  cmix_gibbs_mass
!DEC$ attributes alias: '_cmix_gibbs_mass_' ::  cmix_gibbs_mass
!DEC$ attributes dllimport ::  cmix_gibbs_mass
      end function

      double precision function cmix_cp_mass()
!DEC$ attributes c, reference ::  cmix_cp_mass
!DEC$ attributes alias: '_cmix_cp_mass_' ::  cmix_cp_mass
!DEC$ attributes dllimport ::  cmix_cp_mass
      end function

      double precision function cmix_cv_mass()
!DEC$ attributes c, reference ::  cmix_cv_mass
!DEC$ attributes alias: '_cmix_cv_mass_' ::  cmix_cv_mass
!DEC$ attributes dllimport ::  cmix_cv_mass
      end function

      subroutine cmix_setpe_k(k, pe)
!DEC$ attributes c, reference ::  cmix_setpe_k
!DEC$ attributes alias: '_cmix_setpe_k_' ::  cmix_setpe_k
!DEC$ attributes dllimport ::  cmix_setpe_k
      integer, intent(in) :: k
      double precision, intent(in) :: pe
      end subroutine

      double precision function cmix_pe(k)
!DEC$ attributes c, reference ::  cmix_pe
!DEC$ attributes alias: '_cmix_pe_' ::  cmix_pe
!DEC$ attributes dllimport ::  cmix_pe
      integer, intent(in) :: k
      end function

      double precision function cmix_crittemp()
!DEC$ attributes c, reference ::  cmix_crittemp
!DEC$ attributes alias: '_cmix_crittemp_' ::  cmix_crittemp
!DEC$ attributes dllimport ::  cmix_crittemp
      end function

      double precision function cmix_critpres()
!DEC$ attributes c, reference ::  cmix_critpres
!DEC$ attributes alias: '_cmix_critpres_' ::  cmix_critpres
!DEC$ attributes dllimport ::  cmix_critpres
      end function

      double precision function cmix_sattemp(p)
!DEC$ attributes c, reference ::  cmix_sattemp
!DEC$ attributes alias: '_cmix_sattemp_' ::  cmix_sattemp
!DEC$ attributes dllimport ::  cmix_sattemp
      double precision, intent(in) :: p
      end function

      double precision function cmix_satpres(t)
!DEC$ attributes c, reference ::  cmix_satpres
!DEC$ attributes alias: '_cmix_satpres_' ::  cmix_satpres
!DEC$ attributes dllimport ::  cmix_satpres
      double precision, intent(in) :: t
      end function

      type(EOS) function cmix_get_eos()
!DEC$ attributes c, reference ::  cmix_get_eos
!DEC$ attributes alias: '_cmix_get_eos_' ::  cmix_get_eos
!DEC$ attributes dllimport ::  cmix_get_eos
      end function

      subroutine cmix_set_eos(eos)
!DEC$ attributes c, reference ::  cmix_set_eos
!DEC$ attributes alias: '_cmix_set_eos_' ::  cmix_set_eos
!DEC$ attributes dllimport ::  cmix_set_eos
      type(EOS), intent(inout) :: eos
      end subroutine

      subroutine cmix_equilib(propPair)
!DEC$ attributes c, reference ::  cmix_equilib
!DEC$ attributes alias: '_cmix_equilib_' ::  cmix_equilib
!DEC$ attributes dllimport ::  cmix_equilib
      character*(*), intent(in) :: propPair
      end subroutine

      subroutine cmix_set_trans(transportMgr)
!DEC$ attributes c, reference ::  cmix_set_trans
!DEC$ attributes alias: '_cmix_set_trans_' ::  cmix_set_trans
!DEC$ attributes dllimport ::  cmix_set_trans
      type(Transport), intent(inout) :: transportMgr
      end subroutine

      type(Transport) function cmix_transportMgr()
!DEC$ attributes c, reference ::  cmix_transportMgr
!DEC$ attributes alias: '_cmix_transportmgr_' ::  cmix_transportMgr
!DEC$ attributes dllimport ::  cmix_transportMgr
      end function

      double precision function cmix_visc()
!DEC$ attributes c, reference ::  cmix_visc
!DEC$ attributes alias: '_cmix_visc_' ::  cmix_visc
!DEC$ attributes dllimport ::  cmix_visc
      end function

      subroutine cmix_gspvisc(visc)
!DEC$ attributes c, reference ::  cmix_gspvisc
!DEC$ attributes alias: '_cmix_gspvisc_' ::  cmix_gspvisc
!DEC$ attributes dllimport ::  cmix_gspvisc
      double precision, intent(inout) :: visc(*)
      end subroutine


      subroutine ctrans_gflux(mix_hndl, ndim, ldx, grad_X, 
     &                        grad_T, ldf, fluxes)
!DEC$ attributes c, reference ::  ctrans_gflux
!DEC$ attributes alias: '_ctrans_gflux_' ::  ctrans_gflux
!DEC$ attributes dllimport ::  ctrans_gflux
      integer, intent(in) :: mix_hndl
      integer, intent(in) :: ndim
      integer, intent(in) :: ldx
      double precision, intent(inout) :: grad_X(ldx,*)
      double precision, intent(in) :: grad_T(*)
      integer, intent(in) :: ldf
      double precision, intent(out) :: fluxes(ldf,*)
      end subroutine
      subroutine cmix_gmultidiff(ldim, multiDiff)
!DEC$ attributes c, reference ::  cmix_gmultidiff
!DEC$ attributes alias: '_cmix_gmultidiff_' ::  cmix_gmultidiff
!DEC$ attributes dllimport ::  cmix_gmultidiff
      integer, intent(in) :: ldim
      double precision, intent(in) :: multiDiff(ldim, *)
      end subroutine

      double precision function cmix_tcon()
!DEC$ attributes c, reference ::  cmix_tcon
!DEC$ attributes alias: '_cmix_tcon_' ::  cmix_tcon
!DEC$ attributes dllimport ::  cmix_tcon
      end function

      subroutine cmix_gtdiff(dt)
!DEC$ attributes c, reference ::  cmix_gtdiff
!DEC$ attributes alias: '_cmix_gtdiff_' ::  cmix_gtdiff
!DEC$ attributes dllimport ::  cmix_gtdiff
      double precision, intent(out) :: dt(*)
      end subroutine

      subroutine cmix_getrxnstr(i, rxnString)
!DEC$ attributes c, reference ::  cmix_getrxnstr
!DEC$ attributes alias: '_cmix_getrxnstr_' ::  cmix_getrxnstr
!DEC$ attributes dllimport ::  cmix_getrxnstr
      integer, intent(in) :: i
      character*(*), intent(out) :: rxnString
      end subroutine

      double precision function cmix_rstoich(k, i)
!DEC$ attributes c, reference ::  cmix_rstoich
!DEC$ attributes alias: '_cmix_rstoich_' ::  cmix_rstoich
!DEC$ attributes dllimport ::  cmix_rstoich
      integer, intent(in) :: k
      integer, intent(in) :: i
      end function

      double precision function cmix_pstoich(k, i)
!DEC$ attributes c, reference ::  cmix_pstoich
!DEC$ attributes alias: '_cmix_pstoich_' ::  cmix_pstoich
!DEC$ attributes dllimport ::  cmix_pstoich
      integer, intent(in) :: k
      integer, intent(in) :: i
      end function

      double precision function cmix_nstoich(k, i)
!DEC$ attributes c, reference ::  cmix_nstoich
!DEC$ attributes alias: '_cmix_nstoich_' ::  cmix_nstoich
!DEC$ attributes dllimport ::  cmix_nstoich
      integer, intent(in) :: k
      integer, intent(in) :: i
      end function

      subroutine cmix_get_fwdrop(fwdrop)
!DEC$ attributes c, reference ::  cmix_get_fwdrop
!DEC$ attributes alias: '_cmix_get_fwdrop_' ::  cmix_get_fwdrop
!DEC$ attributes dllimport ::  cmix_get_fwdrop
      double precision, intent(out) :: fwdrop(*)
      end subroutine

      subroutine cmix_get_revrop(revrop)
!DEC$ attributes c, reference ::  cmix_get_revrop
!DEC$ attributes alias: '_cmix_get_revrop_' ::  cmix_get_revrop
!DEC$ attributes dllimport ::  cmix_get_revrop
      double precision, intent(out) :: revrop(*)
      end subroutine

      subroutine cmix_get_netrop(netrop)
!DEC$ attributes c, reference ::  cmix_get_netrop
!DEC$ attributes alias: '_cmix_get_netrop_' ::  cmix_get_netrop
!DEC$ attributes dllimport ::  cmix_get_netrop
      double precision, intent(out) :: netrop(*)
      end subroutine

      subroutine cmix_get_kc(kc)
!DEC$ attributes c, reference ::  cmix_get_kc
!DEC$ attributes alias: '_cmix_get_kc_' ::  cmix_get_kc
!DEC$ attributes dllimport ::  cmix_get_kc
      double precision, intent(out) :: kc(*)
      end subroutine

      subroutine cmix_get_cdot(cdot)
!DEC$ attributes c, reference ::  cmix_get_cdot
!DEC$ attributes alias: '_cmix_get_cdot_' ::  cmix_get_cdot
!DEC$ attributes dllimport ::  cmix_get_cdot
      double precision, intent(out) :: cdot(*)
      end subroutine

      subroutine cmix_get_ddot(ddot)
!DEC$ attributes c, reference ::  cmix_get_ddot
!DEC$ attributes alias: '_cmix_get_ddot_' ::  cmix_get_ddot
!DEC$ attributes dllimport ::  cmix_get_ddot
      double precision, intent(out) :: ddot(*)
      end subroutine

      subroutine cmix_get_wdot(wdot)
!DEC$ attributes c, reference ::  cmix_get_wdot
!DEC$ attributes alias: '_cmix_get_wdot_' ::  cmix_get_wdot
!DEC$ attributes dllimport ::  cmix_get_wdot
      double precision, intent(out) :: wdot(*)
      end subroutine

      end interface
      contains

      type (Mixture) function CKGas(inputfile, thermofile,            &
     &     ivalidate)
      character*(*), intent(in), optional :: inputfile
      character*(*), intent(in), optional :: thermofile
      integer, intent(in), optional :: ivalidate
      type (Mixture) :: mix
      mix%hndl = cmix_newidealgas()
      if (present(inputfile)) then
         if (present(thermofile)) then
            if (present(ivalidate)) then
               call importCK(mix, inputfile, thermofile, ivalidate)
            else
               call importCK(mix, inputfile, thermofile)
            end if
         else
            call importCK(mix, inputfile)
         end if
      end if
      CKGas = mix
      end function CKGas

      type(Mixture) function GRIMech30()
      type(Mixture) :: mix
      mix = CKGas('gri30.inp')
      GRIMech30 = mix
      return
      end function

      subroutine importCK(mix, inputfile, thermofile, ivalidate)
      type (Mixture), intent(inout) :: mix
      character*(*), intent(in) :: inputfile
      character*(*), intent(in), optional :: thermofile
      integer, intent(in), optional :: ivalidate
      character*1 th
      integer iv, iok
      th = ''
      iv = 0
      if (present(thermofile) .and. present(ivalidate)) then
         iok = cmix_import(mix%hndl, inputfile, thermofile, ivalidate)
      else if (present(thermofile)) then
         iok = cmix_import(mix%hndl, inputfile, thermofile, iv)
      else
         iok = cmix_import(mix%hndl, inputfile, th, iv)
      end if
      if (iok .lt. 0) then
         mix%nel = 0
         mix%nsp = 0
         mix%nrxn = 0
      else
         mix%nel = mix_nElements(mix)
         mix%nsp = mix_nSpecies(mix)
         mix%nrxn = mix_nReactions(mix)
      end if
      end subroutine importCK

      type(Mixture) function PureSub(i)
      type(Mixture) :: mix
      mix%hndl = cmix_newsub(i)
      mix%nel = cmix_nElements(mix%hndl)
      mix%nsp = 1
      mix%nrxn = 0
      ! mix%model = -1;
      PureSub = mix
      return
      end function

      type(Mixture) function Water()
      Water = PureSub(0)
      return
      end function

      type(Mixture) function Nitrogen()
      Nitrogen = PureSub(1)
      return
      end function

      type(Mixture) function Methane()
      Methane = PureSub(2)
      return
      end function

      type(Mixture) function Hydrogen()
      Hydrogen = PureSub(3)
      return
      end function

      type(Mixture) function Oxygen()
      Oxygen = PureSub(4)
      return
      end function
      subroutine mix_delete(mix)
       type (Mixture), intent(inout) :: mix
       if (mix%hndl .ne. 0) then
       call cmix_delete(mix%hndl)
       mix%nel = 0
       mix%nsp = 0
       mix%nrxn = 0
       mix%hndl = 0
       end if
       end subroutine
 
      logical function mix_ready(mix)
      type (Mixture), intent(in) :: mix
      if (mix%hndl .eq. 0 .or. mix%nel .eq. 0) then
         mix_ready = .false.
      else
         mix_ready = .true.
      end if
      return
      end function
 
      subroutine mix_copy(dest, src)
      type (Mixture), intent(out) :: dest
      type (Mixture), intent(in) :: src
      dest%nel = src%nel
      dest%nsp = src%nsp
      dest%nrxn = src%nrxn
      dest%hndl = src%hndl
      !call incref(src%hndl)
      end subroutine
!
!     saveState
!
      subroutine mix_saveState(self, state)
      type(Mixture), intent(in) :: self
      double precision, intent(out) :: state(state%nsp + 2)
      call cmix_saveState(self%hndl, state)
      return
      end subroutine
!
!     restoreState
!
      subroutine mix_restoreState(self, state)
      type(Mixture), intent(inout) :: self
      double precision, intent(out) :: state(state%nsp + 2)
      call cmix_restoreState(self%hndl, state)
      return
      end subroutine
!
!     nElements
!
      integer function mix_nElements(self)
      type(Mixture), intent(in) :: self
      mix_nElements=cmix_nElements(self%hndl)
      return
      end function
!
!     nSpecies
!
      integer function mix_nSpecies(self)
      type(Mixture), intent(in) :: self
      mix_nSpecies=cmix_nSpecies(self%hndl)
      return
      end function
!
!     nReactions
!
      integer function mix_nReactions(self)
      type(Mixture), intent(inout) :: self
      mix_nReactions=cmix_nReactions(self%hndl)
      return
      end function

      type (Mixture) function Mixture()
      type(Mixture) mix
      mix%hndl = cmix_newmixture()
      Mixture = mix
      end function
!
!     addElement
!
      subroutine mix_addElement(self, name, atomicWt)
      type(Mixture), intent(inout) :: self
      character*(*), intent(in) :: name
      double precision, intent(in) :: atomicWt
      call cmix_addElement(self%hndl, name, atomicWt)
      return
      end subroutine

      subroutine mix_getElementNames(mix, enames)
      type(Mixture), intent(in) :: mix
      character*(*), intent(out) :: enames(mix%nel)
      do i = 1, mix%nel
         call cmix_getename(mix%hndl, i, enames(i))
      end do
      end subroutine

      integer function mix_eindex(mix, name)
      type(Mixture), intent(in) :: mix
      character*(*), intent(in) :: name   ! element name
      mix_eindex = cmix_eindex(mix%hndl, name)
      end function
!
!     atomicWeight
!
      double precision function mix_atwt(self, m)
      type(Mixture), intent(inout) :: self
      integer, intent(in) :: m
      mix_atwt=cmix_atwt(self%hndl, m)
      return
      end function
!
!     nAtoms
!
      double precision function mix_nAtoms(self, k, m)
      type(Mixture), intent(inout) :: self
      integer, intent(in) :: k
      integer, intent(in) :: m
      mix_nAtoms=cmix_nAtoms(self%hndl, k, m)
      return
      end function
!
!     charge
!
      double precision function mix_charge(self, k)
      type(Mixture), intent(inout) :: self
      integer, intent(in) :: k
      mix_charge=cmix_charge(self%hndl, k)
      return
      end function
!
!     speciesIndex
!
      integer function mix_speciesIndex(self, name)
      type(Mixture), intent(inout) :: self
      character*(*), intent(in) :: name
      mix_speciesIndex=cmix_speciesIndex(self%hndl, name)
      return
      end function

      subroutine mix_getspnm(mix, snames)
      type(Mixture), intent(in) :: mix
      character*(*), intent(out) :: snames(mix%nsp)
      do i = 1, mix%nsp
         call cmix_getspnm(mix%hndl, i, snames(i))
      end do
      end subroutine
!
!     molecularWeight
!
      double precision function mix_molwt(self, k)
      type(Mixture), intent(inout) :: self
      integer, intent(in) :: k
      mix_molwt=cmix_molwt(self%hndl, k)
      return
      end function
!
!     getMolecularWeights
!
      subroutine mix_gmolwts(self, molWts)
      type(Mixture), intent(inout) :: self
      double precision, intent(out) :: molWts(molWts%nsp)
      call cmix_gmolwts(self%hndl, molWts)
      return
      end subroutine
!
!     setMoleFractions
!
      subroutine mix_setx(self, moleFracs)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: moleFracs(moleFracs%nsp)
      call cmix_setx(self%hndl, moleFracs)
      return
      end subroutine
!
!     setMoleFractions_NoNorm
!
      subroutine mix_sxnonorm(self, moleFracs)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: moleFracs(moleFracs%nsp)
      call cmix_sxnonorm(self%hndl, moleFracs)
      return
      end subroutine
!
!     setMassFractions
!
      subroutine mix_sety(self, massFracs)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: massFracs(massFracs%nsp)
      call cmix_sety(self%hndl, massFracs)
      return
      end subroutine
!
!     setMassFractions_NoNorm
!
      subroutine mix_synonorm(self, y)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: y(y%nsp)
      call cmix_synonorm(self%hndl, y)
      return
      end subroutine
!
!     setConcentrations
!
      subroutine mix_sconc(self, conc)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: conc(conc%nsp)
      call cmix_sconc(self%hndl, conc)
      return
      end subroutine
!
!     moleFraction
!
      double precision function mix_xbyname(self, name)
      type(Mixture), intent(inout) :: self
      character*(*), intent(in) :: name
      mix_xbyname=cmix_xbyname(self%hndl, name)
      return
      end function
!
!     massFraction
!
      double precision function mix_ybyname(self, name)
      type(Mixture), intent(inout) :: self
      character*(*), intent(in) :: name
      mix_ybyname=cmix_ybyname(self%hndl, name)
      return
      end function
!
!     getMoleFractions
!
      subroutine mix_getx(self, moleFracs)
      type(Mixture), intent(inout) :: self
      double precision, intent(inout) :: moleFracs(moleFracs%nsp)
      call cmix_getx(self%hndl, moleFracs)
      return
      end subroutine
!
!     getMassFractions
!
      subroutine mix_gety(self, massFracs)
      type(Mixture), intent(inout) :: self
      double precision, intent(inout) :: massFracs(massFracs%nsp)
      call cmix_gety(self%hndl, massFracs)
      return
      end subroutine
!
!     getConcentrations
!
      subroutine mix_getconc(self, concentrations)
      type(Mixture), intent(inout) :: self
      double precision, intent(inout) :: concentrations(concentrations%nsp)
      call cmix_getconc(self%hndl, concentrations)
      return
      end subroutine
!
!     mean_X
!
      double precision function mix_mean_X(self, Q)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: Q(Q%nsp)
      mix_mean_X=cmix_mean_X(self%hndl, Q)
      return
      end function
!
!     mean_Y
!
      double precision function mix_mean_Y(self, Q)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: Q(Q%nsp)
      mix_mean_Y=cmix_mean_Y(self%hndl, Q)
      return
      end function
!
!     meanMolecularWeight
!
      double precision function mix_meanmw(self)
      type(Mixture), intent(inout) :: self
      mix_meanmw=cmix_meanmw(self%hndl)
      return
      end function
!
!     sum_xlogQ
!
      double precision function mix_sum_xlogQ(self, Q)
      type(Mixture), intent(inout) :: self
      double precision, intent(out) :: Q(Q%nsp)
      mix_sum_xlogQ=cmix_sum_xlogQ(self%hndl, Q)
      return
      end function
!
!     setState_TPX
!
      subroutine mix_setTPXarray(self, t, p, moleFracs)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      double precision, intent(in) :: moleFracs(t%nsp)
      call cmix_setTPXarray(self%hndl, t, p, moleFracs)
      return
      end subroutine
!
!     setState_TPX
!
      subroutine mix_setTPXstr(self, t, p, moleFracs)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      character*(*), intent(in) :: moleFracs
      call cmix_setTPXstr(self%hndl, t, p, moleFracs)
      return
      end subroutine
!
!     setState_TPY
!
      subroutine mix_setTPYarray(self, t, p, massFracs)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      double precision, intent(in) :: massFracs(t%nsp)
      call cmix_setTPYarray(self%hndl, t, p, massFracs)
      return
      end subroutine
!
!     setState_TPY
!
      subroutine mix_setTPYstr(self, t, p, massFracs)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      character*(*), intent(in) :: massFracs
      call cmix_setTPYstr(self%hndl, t, p, massFracs)
      return
      end subroutine
!
!     setState_TRX
!
      subroutine mix_setState_TRX(self, t, density, moleFracs)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: density
      double precision, intent(in) :: moleFracs(t%nsp)
      call cmix_setState_TRX(self%hndl, t, density, moleFracs)
      return
      end subroutine
!
!     setState_TRY
!
      subroutine mix_setState_TRY(self, t, density, massFracs)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: density
      double precision, intent(in) :: massFracs(t%nsp)
      call cmix_setState_TRY(self%hndl, t, density, massFracs)
      return
      end subroutine
!
!     setTemperature
!
      subroutine mix_settemp(self, temp)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: temp
      call cmix_settemp(self%hndl, temp)
      return
      end subroutine
!
!     setDensity
!
      subroutine mix_setDensity(self, density)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: density
      call cmix_setDensity(self%hndl, density)
      return
      end subroutine
!
!     setState_TP
!
      subroutine mix_setState_TP(self, t, p)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      call cmix_setState_TP(self%hndl, t, p)
      return
      end subroutine
!
!     setState_PX
!
      subroutine mix_setState_PX(self, p, x)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: p
      double precision, intent(in) :: x(p%nsp)
      call cmix_setState_PX(self%hndl, p, x)
      return
      end subroutine
!
!     setState_PY
!
      subroutine mix_setState_PY(self, p, massFracs)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: p
      double precision, intent(in) :: massFracs(p%nsp)
      call cmix_setState_PY(self%hndl, p, massFracs)
      return
      end subroutine
!
!     setPressure
!
      subroutine mix_setPressure(self, p)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: p
      call cmix_setPressure(self%hndl, p)
      return
      end subroutine
!
!     setState_TR
!
      subroutine mix_setState_TR(self, t, rho)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: rho
      call cmix_setState_TR(self%hndl, t, rho)
      return
      end subroutine
!
!     setState_TX
!
      subroutine mix_setState_TX(self, t, moleFracs)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(inout) :: moleFracs(t%nsp)
      call cmix_setState_TX(self%hndl, t, moleFracs)
      return
      end subroutine
!
!     setState_TY
!
      subroutine mix_setState_TY(self, t, massFracs)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: massFracs(t%nsp)
      call cmix_setState_TY(self%hndl, t, massFracs)
      return
      end subroutine
!
!     setState_RX
!
      subroutine mix_setState_RX(self, density, moleFracs)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: density
      double precision, intent(in) :: moleFracs(density%nsp)
      call cmix_setState_RX(self%hndl, density, moleFracs)
      return
      end subroutine
!
!     setState_RY
!
      subroutine mix_setState_RY(self, density, y)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: density
      double precision, intent(in) :: y(density%nsp)
      call cmix_setState_RY(self%hndl, density, y)
      return
      end subroutine
!
!     setState_HP
!
      subroutine mix_setState_HP(self, t, p)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      call cmix_setState_HP(self%hndl, t, p)
      return
      end subroutine
!
!     setState_UV
!
      subroutine mix_setState_UV(self, t, p)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      call cmix_setState_UV(self%hndl, t, p)
      return
      end subroutine
!
!     setState_SP
!
      subroutine mix_setState_SP(self, t, p)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      call cmix_setState_SP(self%hndl, t, p)
      return
      end subroutine
!
!     setState_SV
!
      subroutine mix_setState_SV(self, t, p)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      call cmix_setState_SV(self%hndl, t, p)
      return
      end subroutine
!
!     getEnthalpy_RT
!
      subroutine mix_gh_rt(self, h_RT)
      type(Mixture), intent(inout) :: self
      double precision, intent(out) :: h_RT(h_RT%nsp)
      call cmix_gh_rt(self%hndl, h_RT)
      return
      end subroutine
!
!     getGibbs_RT
!
      subroutine mix_ggibbs_rt(self, g_RT)
      type(Mixture), intent(inout) :: self
      double precision, intent(out) :: g_RT(g_RT%nsp)
      call cmix_ggibbs_rt(self%hndl, g_RT)
      return
      end subroutine
!
!     getCp_R
!
      subroutine mix_gcp_r(self, cp_R)
      type(Mixture), intent(inout) :: self
      double precision, intent(out) :: cp_R(cp_R%nsp)
      call cmix_gcp_r(self%hndl, cp_R)
      return
      end subroutine
!
!     getEntropy_R
!
      subroutine mix_gs_r(self, s_R)
      type(Mixture), intent(inout) :: self
      double precision, intent(out) :: s_R(s_R%nsp)
      call cmix_gs_r(self%hndl, s_R)
      return
      end subroutine
!
!     minTemp
!
      double precision function mix_minTemp(self)
      type(Mixture), intent(inout) :: self
      mix_minTemp=cmix_minTemp(self%hndl)
      return
      end function
!
!     maxTemp
!
      double precision function mix_maxTemp(self)
      type(Mixture), intent(inout) :: self
      mix_maxTemp=cmix_maxTemp(self%hndl)
      return
      end function
!
!     refPressure
!
      double precision function mix_refp(self)
      type(Mixture), intent(inout) :: self
      mix_refp=cmix_refp(self%hndl)
      return
      end function
!
!     temperature
!
      double precision function mix_temperature(self)
      type(Mixture), intent(inout) :: self
      mix_temperature=cmix_temperature(self%hndl)
      return
      end function
!
!     density
!
      double precision function mix_density(self)
      type(Mixture), intent(inout) :: self
      mix_density=cmix_density(self%hndl)
      return
      end function
!
!     molarDensity
!
      double precision function mix_moldens(self)
      type(Mixture), intent(inout) :: self
      mix_moldens=cmix_moldens(self%hndl)
      return
      end function
!
!     pressure
!
      double precision function mix_pressure(self)
      type(Mixture), intent(inout) :: self
      mix_pressure=cmix_pressure(self%hndl)
      return
      end function
!
!     enthalpy_mole
!
      double precision function mix_enthalpy_mole(self)
      type(Mixture), intent(inout) :: self
      mix_enthalpy_mole=cmix_enthalpy_mole(self%hndl)
      return
      end function
!
!     intEnergy_mole
!
      double precision function mix_umole(self)
      type(Mixture), intent(inout) :: self
      mix_umole=cmix_umole(self%hndl)
      return
      end function
!
!     entropy_mole
!
      double precision function mix_entropy_mole(self)
      type(Mixture), intent(inout) :: self
      mix_entropy_mole=cmix_entropy_mole(self%hndl)
      return
      end function
!
!     gibbs_mole
!
      double precision function mix_gibbs_mole(self)
      type(Mixture), intent(inout) :: self
      mix_gibbs_mole=cmix_gibbs_mole(self%hndl)
      return
      end function
!
!     cp_mole
!
      double precision function mix_cp_mole(self)
      type(Mixture), intent(inout) :: self
      mix_cp_mole=cmix_cp_mole(self%hndl)
      return
      end function
!
!     cv_mole
!
      double precision function mix_cv_mole(self)
      type(Mixture), intent(inout) :: self
      mix_cv_mole=cmix_cv_mole(self%hndl)
      return
      end function
!
!     getChemPotentials_RT
!
      subroutine mix_gchempot(self, mu)
      type(Mixture), intent(inout) :: self
      double precision, intent(out) :: mu(mu%nsp)
      call cmix_gchempot(self%hndl, mu)
      return
      end subroutine
!
!     enthalpy_mass
!
      double precision function mix_enthalpy_mass(self)
      type(Mixture), intent(inout) :: self
      mix_enthalpy_mass=cmix_enthalpy_mass(self%hndl)
      return
      end function
!
!     intEnergy_mass
!
      double precision function mix_umass(self)
      type(Mixture), intent(inout) :: self
      mix_umass=cmix_umass(self%hndl)
      return
      end function
!
!     entropy_mass
!
      double precision function mix_entropy_mass(self)
      type(Mixture), intent(inout) :: self
      mix_entropy_mass=cmix_entropy_mass(self%hndl)
      return
      end function
!
!     gibbs_mass
!
      double precision function mix_gibbs_mass(self)
      type(Mixture), intent(inout) :: self
      mix_gibbs_mass=cmix_gibbs_mass(self%hndl)
      return
      end function
!
!     cp_mass
!
      double precision function mix_cp_mass(self)
      type(Mixture), intent(inout) :: self
      mix_cp_mass=cmix_cp_mass(self%hndl)
      return
      end function
!
!     cv_mass
!
      double precision function mix_cv_mass(self)
      type(Mixture), intent(inout) :: self
      mix_cv_mass=cmix_cv_mass(self%hndl)
      return
      end function
!
!     setPotentialEnergy
!
      subroutine mix_setpe_k(self, k, pe)
      type(Mixture), intent(inout) :: self
      integer, intent(in) :: k
      double precision, intent(in) :: pe
      call cmix_setpe_k(self%hndl, k, pe)
      return
      end subroutine
!
!     potentialEnergy
!
      double precision function mix_pe(self, k)
      type(Mixture), intent(inout) :: self
      integer, intent(in) :: k
      mix_pe=cmix_pe(self%hndl, k)
      return
      end function
!
!     critTemperature
!
      double precision function mix_crittemp(self)
      type(Mixture), intent(inout) :: self
      mix_crittemp=cmix_crittemp(self%hndl)
      return
      end function
!
!     critPressure
!
      double precision function mix_critpres(self)
      type(Mixture), intent(inout) :: self
      mix_critpres=cmix_critpres(self%hndl)
      return
      end function
!
!     satTemperature
!
      double precision function mix_sattemp(self, p)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: p
      mix_sattemp=cmix_sattemp(self%hndl, p)
      return
      end function
!
!     satPressure
!
      double precision function mix_satpres(self, t)
      type(Mixture), intent(inout) :: self
      double precision, intent(in) :: t
      mix_satpres=cmix_satpres(self%hndl, t)
      return
      end function

      type(EOS) function thrm_get_eos(mix)
      type(Mixture), intent(inout) :: mix
      type(EOS) e
      e%hndl = cthrm_get_eos(mix%hndl)
      thrm_get_eos = e
      return
      end function
!
!     setEquationOfState
!
      subroutine mix_set_eos(self, eos)
      type(Mixture), intent(inout) :: self
      type(EOS), intent(inout) :: eos
      call cmix_set_eos(self%hndl, eos%hndl)
      return
      end subroutine
!
!     equilibrate
!
      subroutine mix_equilib(self, propPair)
      type(Mixture), intent(inout) :: self
      character*(*), intent(in) :: propPair
      call cmix_equilib(self%hndl, propPair)
      return
      end subroutine

      subroutine trans_set_trans(self, transportMgr)
      type(Mixture), intent(inout) :: self
      type(Transport), intent(inout) :: transportMgr
      call ctrans_set_trans(self%hndl, transportMgr%hndl)
      return
      end subroutine

      type(Transport) function trans_transportMgr(self)
      type(Mixture), intent(inout) :: self
      type(Transport) tr
      tr%hndl = ctrans_transportmgr(self%hndl)
      tr%ierr  = 0
      trans_transportMgr = tr
      return
      end function
!
!     viscosity
!
      double precision function mix_visc(self)
      type(Mixture), intent(inout) :: self
      mix_visc=cmix_visc(self%hndl)
      return
      end function
!
!     getSpeciesViscosities
!
      subroutine mix_gspvisc(self, visc)
      type(Mixture), intent(inout) :: self
      double precision, intent(inout) :: visc(visc%nsp)
      call cmix_gspvisc(self%hndl, visc)
      return
      end subroutine

!
!     getSpeciesFluxes
!
      subroutine trans_gflux(mix, ndim, grad_X, grad_T, fluxes)
      type(Mixture), intent(in) :: mix
      integer, intent(in) :: ndim
      double precision, intent(inout) :: grad_X(:,:)
      double precision, intent(in) :: grad_T(:)
      double precision, intent(out) :: fluxes(:,:)
      ldx = size(grad_X, dim=1)
      ldf = size(fluxes, dim=1)
      if (ldx .lt. mix%nsp .or. ldf .lt. mix%nsp) then
	write(*,*) '### Error in getSpeciesFluxes ###'
	write(*,*) 'An array is too small.'
        write(*,*) 'minimum array sizes:'
        write(*,*) '  grad_X(',mix%nsp,',',ndim,')'
        write(*,*) '  grad_T(',ndim,')'
        write(*,*) '  fluxes(',mix%nsp,',',ndim,')'
	stop
      end if
      call ctrans_gflux(mix%hndl, ndim, ldx, grad_X, 
     &                            grad_T, ldf, fluxes)
      return
      end subroutine
!
!     getMultiDiffCoeffs
!
      subroutine mix_gmultidiff(self, ldim, multiDiff)
      type(Mixture), intent(inout) :: self
      integer, intent(in) :: ldim
      double precision, intent(in) :: multiDiff(ldim, ldim%nsp)
      call cmix_gmultidiff(self%hndl, ldim, multiDiff)
      return
      end subroutine
!
!     thermalConductivity
!
      double precision function mix_tcon(self)
      type(Mixture), intent(inout) :: self
      mix_tcon=cmix_tcon(self%hndl)
      return
      end function
!
!     getThermalDiffCoeffs
!
      subroutine mix_gtdiff(self, dt)
      type(Mixture), intent(inout) :: self
      double precision, intent(out) :: dt(dt%nsp)
      call cmix_gtdiff(self%hndl, dt)
      return
      end subroutine
!
!     getReactionString
!
      subroutine mix_getrxnstr(self, i, rxnString)
      type(Mixture), intent(inout) :: self
      integer, intent(in) :: i
      character*(*), intent(out) :: rxnString
      call cmix_getrxnstr(self%hndl, i, rxnString)
      return
      end subroutine
!
!     reactantStoichCoeff
!
      double precision function mix_rstoich(self, k, i)
      type(Mixture), intent(inout) :: self
      integer, intent(in) :: k
      integer, intent(in) :: i
      mix_rstoich=cmix_rstoich(self%hndl, k, i)
      return
      end function
!
!     productStoichCoeff
!
      double precision function mix_pstoich(self, k, i)
      type(Mixture), intent(inout) :: self
      integer, intent(in) :: k
      integer, intent(in) :: i
      mix_pstoich=cmix_pstoich(self%hndl, k, i)
      return
      end function
!
!     netStoichCoeff
!
      double precision function mix_nstoich(self, k, i)
      type(Mixture), intent(inout) :: self
      integer, intent(in) :: k
      integer, intent(in) :: i
      mix_nstoich=cmix_nstoich(self%hndl, k, i)
      return
      end function
!
!     getFwdRatesOfProgress
!
      subroutine mix_get_fwdrop(self, fwdrop)
      type(Mixture), intent(inout) :: self
      double precision, intent(out) :: fwdrop(fwdrop%nrxn)
      call cmix_get_fwdrop(self%hndl, fwdrop)
      return
      end subroutine
!
!     getRevRatesOfProgress
!
      subroutine mix_get_revrop(self, revrop)
      type(Mixture), intent(inout) :: self
      double precision, intent(out) :: revrop(revrop%nrxn)
      call cmix_get_revrop(self%hndl, revrop)
      return
      end subroutine
!
!     getNetRatesOfProgress
!
      subroutine mix_get_netrop(self, netrop)
      type(Mixture), intent(inout) :: self
      double precision, intent(out) :: netrop(netrop%nrxn)
      call cmix_get_netrop(self%hndl, netrop)
      return
      end subroutine
!
!     getEquilibriumConstants
!
      subroutine mix_get_kc(self, kc)
      type(Mixture), intent(inout) :: self
      double precision, intent(out) :: kc(kc%nsp)
      call cmix_get_kc(self%hndl, kc)
      return
      end subroutine
!
!     getCreationRates
!
      subroutine mix_get_cdot(self, cdot)
      type(Mixture), intent(inout) :: self
      double precision, intent(out) :: cdot(cdot%nsp)
      call cmix_get_cdot(self%hndl, cdot)
      return
      end subroutine
!
!     getDestructionRates
!
      subroutine mix_get_ddot(self, ddot)
      type(Mixture), intent(inout) :: self
      double precision, intent(out) :: ddot(ddot%nsp)
      call cmix_get_ddot(self%hndl, ddot)
      return
      end subroutine
!
!     getNetProductionRates
!
      subroutine mix_get_wdot(self, wdot)
      type(Mixture), intent(inout) :: self
      double precision, intent(out) :: wdot(wdot%nsp)
      call cmix_get_wdot(self%hndl, wdot)
      return
      end subroutine
      end module
