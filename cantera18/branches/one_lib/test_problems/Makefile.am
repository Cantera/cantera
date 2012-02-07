
#SUBDIRS    = cxx_ex surfkin fracCoeff silane_equil VPsilance_test \
	     negATest diamondSurf diamondSurf_dupl surfSolverTest \
             ChemEquil_gri_matrix ChemEquil_gri_pairs             \
	     mixGasTransport multiGasTransport printUtilUnitTest 

SUBDIRS  = cathermo mixGasTransport ChemEquil_gri_matrix ChemEquil_gri_pairs
SUBDIRS += ChemEquil_ionizedGas ChemEquil_red1 CpJump


#
# would be good to add:
#
# pureFluidTest
# rankine_democxx
# min_python
# 
# nasa9_reader
# cathermo
# VCSnonideal


# skipped:
# ck2cti_test


EXTRA_DIST =

dist-hook:
	rm -rf `find $(distdir)/ -name .svn`
	rm -rf `find $(distdir)/ -name .deps`