h_sources  = ck2ct.h CKParser.h CKReader.h Constituent.h Reaction.h \
	     Species.h writelog.h ck2ctml.h ckr_defs.h ckr_utils.h  \
	     Element.h  RxnSpecies.h thermoFunctions.h

cc_sources = atomicWeightDB.cpp ck2ctml.cpp CKReader.cpp filter.cpp  \
	     Reaction.cpp thermoFunctions.cpp ck2ct.cpp CKParser.cpp \
	     ckr_utils.cpp NASA9Parser.cpp Species.cpp writelog.cpp

INCLUDES                = -I../base -I../../../
AM_CXXFLAGS             = $(INCLUDES)

lib_LTLIBRARIES         = libconverters.la
library_includedir      = $(includedir)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

libconverters_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
libconverters_la_SOURCES      = $(cc_sources) $(h_sources)