
#SUBDIRS    = surfkin VPsilance_test \
	     diamondSurf_dupl surfSolverTest

SUBDIRS  = cathermo mixGasTransport ChemEquil_gri_matrix ChemEquil_gri_pairs
SUBDIRS += ChemEquil_ionizedGas ChemEquil_red1 CpJump cxx_ex diamondSurf 
SUBDIRS += diamondSurf_dupl fracCoeff multiGasTransport NASA9poly_test
SUBDIRS += negATest printUtilUnitTest pureFluidTest silane_equil 


# skipped:
# ck2cti_test
# nasa9_reader (uses ck2cti)
# python
# rankine_democxx

EXTRA_DIST =

dist-hook:
	rm -rf `find $(distdir)/ -name .svn`
	rm -rf `find $(distdir)/ -name .deps`