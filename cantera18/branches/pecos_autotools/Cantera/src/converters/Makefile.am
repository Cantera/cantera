h_sources  = ck2ct.h CKParser.h CKReader.h Constituent.h Reaction.h \
	     Species.h writelog.h ck2ctml.h ckr_defs.h ckr_utils.h  \
	     Element.h RxnSpecies.h thermoFunctions.h

cc_sources = atomicWeightDB.cpp ck2ctml.cpp CKReader.cpp filter.cpp  \
	     Reaction.cpp thermoFunctions.cpp ck2ct.cpp CKParser.cpp \
	     ckr_utils.cpp NASA9Parser.cpp Species.cpp writelog.cpp

AM_CPPFLAGS             = -I../base -I../../../
AM_CXXFLAGS             = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libconverters.la
library_includedir      = $(includedir)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libconverters_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libconverters_la_SOURCES      = $(cc_sources) $(h_sources)

# header file accumulation
all:
	@echo copying headers $<
	cp -u $(h_sources) $(top_builddir)/build/include/
