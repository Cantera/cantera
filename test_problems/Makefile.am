
SUBDIRS  = statmech cathermo mixGasTransport ChemEquil_gri_matrix ChemEquil_gri_pairs
SUBDIRS += ChemEquil_ionizedGas ChemEquil_red1 CpJump cxx_ex diamondSurf 
SUBDIRS += diamondSurf_dupl fracCoeff multiGasTransport NASA9poly_test
SUBDIRS += negATest printUtilUnitTest pureFluidTest silane_equil 
SUBDIRS += spectroscopy surfkin surfSolverTest VPsilane_test 
SUBDIRS += VCSnonideal PecosTransport

#f77test

# skipped (uses ck2cti or python):
# python
# rankine_democxx
# VCSnonideal (unsure of problem)

EXTRA_DIST =

dist-hook:
	rm -rf `find $(distdir)/ -name .svn`
	rm -rf `find $(distdir)/ -name .deps`