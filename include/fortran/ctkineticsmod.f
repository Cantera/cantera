      module ctkinetics

      use ctmixture


      use ctmixture
      interface
      subroutine ckin_getrxnstr(mix_hndl, i, rxnString)
!DEC$ attributes c, reference ::  ckin_getrxnstr
!DEC$ attributes alias: '_ckin_getrxnstr_' ::  ckin_getrxnstr
!DEC$ attributes dllimport ::  ckin_getrxnstr
      integer, intent(in) :: mix_hndl
      integer, intent(in) :: i
      character*(*), intent(out) :: rxnString
      end subroutine

      double precision function ckin_rstoich(mix_hndl, k, i)
!DEC$ attributes c, reference ::  ckin_rstoich
!DEC$ attributes alias: '_ckin_rstoich_' ::  ckin_rstoich
!DEC$ attributes dllimport ::  ckin_rstoich
      integer, intent(in) :: mix_hndl
      integer, intent(in) :: k
      integer, intent(in) :: i
      end function

      double precision function ckin_pstoich(mix_hndl, k, i)
!DEC$ attributes c, reference ::  ckin_pstoich
!DEC$ attributes alias: '_ckin_pstoich_' ::  ckin_pstoich
!DEC$ attributes dllimport ::  ckin_pstoich
      integer, intent(in) :: mix_hndl
      integer, intent(in) :: k
      integer, intent(in) :: i
      end function

      double precision function ckin_nstoich(mix_hndl, k, i)
!DEC$ attributes c, reference ::  ckin_nstoich
!DEC$ attributes alias: '_ckin_nstoich_' ::  ckin_nstoich
!DEC$ attributes dllimport ::  ckin_nstoich
      integer, intent(in) :: mix_hndl
      integer, intent(in) :: k
      integer, intent(in) :: i
      end function

      subroutine ckin_get_fwdrop(mix_hndl, fwdrop)
!DEC$ attributes c, reference ::  ckin_get_fwdrop
!DEC$ attributes alias: '_ckin_get_fwdrop_' ::  ckin_get_fwdrop
!DEC$ attributes dllimport ::  ckin_get_fwdrop
      integer, intent(in) :: mix_hndl
      double precision, intent(out) :: fwdrop(*)
      end subroutine

      subroutine ckin_get_revrop(mix_hndl, revrop)
!DEC$ attributes c, reference ::  ckin_get_revrop
!DEC$ attributes alias: '_ckin_get_revrop_' ::  ckin_get_revrop
!DEC$ attributes dllimport ::  ckin_get_revrop
      integer, intent(in) :: mix_hndl
      double precision, intent(out) :: revrop(*)
      end subroutine

      subroutine ckin_get_netrop(mix_hndl, netrop)
!DEC$ attributes c, reference ::  ckin_get_netrop
!DEC$ attributes alias: '_ckin_get_netrop_' ::  ckin_get_netrop
!DEC$ attributes dllimport ::  ckin_get_netrop
      integer, intent(in) :: mix_hndl
      double precision, intent(out) :: netrop(*)
      end subroutine

      subroutine ckin_get_kc(mix_hndl, kc)
!DEC$ attributes c, reference ::  ckin_get_kc
!DEC$ attributes alias: '_ckin_get_kc_' ::  ckin_get_kc
!DEC$ attributes dllimport ::  ckin_get_kc
      integer, intent(in) :: mix_hndl
      double precision, intent(out) :: kc(*)
      end subroutine

      subroutine ckin_get_cdot(mix_hndl, cdot)
!DEC$ attributes c, reference ::  ckin_get_cdot
!DEC$ attributes alias: '_ckin_get_cdot_' ::  ckin_get_cdot
!DEC$ attributes dllimport ::  ckin_get_cdot
      integer, intent(in) :: mix_hndl
      double precision, intent(out) :: cdot(*)
      end subroutine

      subroutine ckin_get_ddot(mix_hndl, ddot)
!DEC$ attributes c, reference ::  ckin_get_ddot
!DEC$ attributes alias: '_ckin_get_ddot_' ::  ckin_get_ddot
!DEC$ attributes dllimport ::  ckin_get_ddot
      integer, intent(in) :: mix_hndl
      double precision, intent(out) :: ddot(*)
      end subroutine

      subroutine ckin_get_wdot(mix_hndl, wdot)
!DEC$ attributes c, reference ::  ckin_get_wdot
!DEC$ attributes alias: '_ckin_get_wdot_' ::  ckin_get_wdot
!DEC$ attributes dllimport ::  ckin_get_wdot
      integer, intent(in) :: mix_hndl
      double precision, intent(out) :: wdot(*)
      end subroutine

      end interface
      contains
!
!     getReactionString
!
      subroutine kin_getrxnstr(mix, i, rxnString)
      type(mixture_t), intent(in) :: mix
      integer, intent(in) :: i
      character*(*), intent(out) :: rxnString
      call ckin_getrxnstr(mix%hndl, i, rxnString)
      return
      end subroutine
!
!     reactantStoichCoeff
!
      double precision function kin_rstoich(mix, k, i)
      type(mixture_t), intent(in) :: mix
      integer, intent(in) :: k
      integer, intent(in) :: i
      kin_rstoich=ckin_rstoich(mix%hndl, k, i)
      return
      end function
!
!     productStoichCoeff
!
      double precision function kin_pstoich(mix, k, i)
      type(mixture_t), intent(in) :: mix
      integer, intent(in) :: k
      integer, intent(in) :: i
      kin_pstoich=ckin_pstoich(mix%hndl, k, i)
      return
      end function
!
!     netStoichCoeff
!
      double precision function kin_nstoich(mix, k, i)
      type(mixture_t), intent(in) :: mix
      integer, intent(in) :: k
      integer, intent(in) :: i
      kin_nstoich=ckin_nstoich(mix%hndl, k, i)
      return
      end function
!
!     getFwdRatesOfProgress
!
      subroutine kin_get_fwdrop(mix, fwdrop)
      type(mixture_t), intent(in) :: mix
      double precision, intent(out) :: fwdrop(mix%nrxn)
      call ckin_get_fwdrop(mix%hndl, fwdrop)
      return
      end subroutine
!
!     getRevRatesOfProgress
!
      subroutine kin_get_revrop(mix, revrop)
      type(mixture_t), intent(in) :: mix
      double precision, intent(out) :: revrop(mix%nrxn)
      call ckin_get_revrop(mix%hndl, revrop)
      return
      end subroutine
!
!     getNetRatesOfProgress
!
      subroutine kin_get_netrop(mix, netrop)
      type(mixture_t), intent(in) :: mix
      double precision, intent(out) :: netrop(mix%nrxn)
      call ckin_get_netrop(mix%hndl, netrop)
      return
      end subroutine
!
!     getEquilibriumConstants
!
      subroutine kin_get_kc(mix, kc)
      type(mixture_t), intent(in) :: mix
      double precision, intent(out) :: kc(mix%nsp)
      call ckin_get_kc(mix%hndl, kc)
      return
      end subroutine
!
!     getCreationRates
!
      subroutine kin_get_cdot(mix, cdot)
      type(mixture_t), intent(in) :: mix
      double precision, intent(out) :: cdot(mix%nsp)
      call ckin_get_cdot(mix%hndl, cdot)
      return
      end subroutine
!
!     getDestructionRates
!
      subroutine kin_get_ddot(mix, ddot)
      type(mixture_t), intent(in) :: mix
      double precision, intent(out) :: ddot(mix%nsp)
      call ckin_get_ddot(mix%hndl, ddot)
      return
      end subroutine
!
!     getNetProductionRates
!
      subroutine kin_get_wdot(mix, wdot)
      type(mixture_t), intent(in) :: mix
      double precision, intent(out) :: wdot(mix%nsp)
      call ckin_get_wdot(mix%hndl, wdot)
      return
      end subroutine
      end module
