
#SUBDIRS    = cxx_ex surfkin fracCoeff silane_equil VPsilance_test \
	     negATest diamondSurf diamondSurf_dupl surfSolverTest \
             ChemEquil_gri_matrix ChemEquil_gri_pairs             \
             ChemEquil_ionizedGas ChemEquil_red1 CpJump           \
	     mixGasTransport multiGasTransport printUtilUnitTest 

SUBDIRS    = cathermo
#
# would be good to add:
#
# pureFluidTest
# rankine_democxx
# min_python
# ck2cti_test
# nasa9_reader
# cathermo
# VCSnonideal

EXTRA_DIST =

dist-hook:
	rm -rf `find $(distdir)/ -name .svn`
	rm -rf `find $(distdir)/ -name .deps`