cc_sources = fct.cpp fctxml.cpp

fc_sources = fct_interface.f90 fctxml_interface.f90 cantera_xml.f90 cantera_thermo.f90 \
	     cantera_kinetics.f90 cantera_transport.f90 cantera_iface.f90              \
	     cantera_funcs.f90 cantera.f90

h_sources               = flib_defs.h

INC                     = -I./ -I$(top_builddir)/build/include/
AM_CPPFLAGS             = $(INC)
AM_CXXFLAGS             = $(AM_CPPFLAGS)
AM_FCFLAGS              = $(AM_CPPFLAGS)

LINK                    = -luser -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo
LINK                   += -lctnumerics -lctmath -ltpx -lctspectra -lconverters -lctbase -lcvode
LINK                   += -lctlapack -lctblas -lctf2c -lctcxx -lm  -lctf2c -lstdc++

AM_LDFLAGS              = -L$(top_builddir)/build/lib/

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libfct.la
library_includedir      = $(top_builddir)/build/include
library_include_HEADERS = $(h_sources)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libfct_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libfct_la_SOURCES      = $(cc_sources) $(fc_sources) $(h_sources)

__top_builddir__build_lib_libfct_include_HEADERS = $(h_sources)
__top_builddir__build_lib_libfct_includedir      = $(prefix)/include


#------------------------                                                                                
# install fo'tran modules                                                                                
#------------------------                                                                                

moduledir          = @libdir@
module_DATA        = cantera_funcs.mod cantera_iface.mod cantera_kinetics.mod \
                     cantera.mod cantera_thermo.mod cantera_transport.mod     \
                     cantera_xml.mod fct.mod fctxml.mod


# install-data-local: $(top_builddir)/build/lib/libfct.la
# 	@echo installing .mod files
# 	$(INSTALL_DATA) $(top_builddir)/Cantera/fortran/src/cantera.mod $(top_builddir)/build/include/cantera.mod
# 	$(INSTALL_DATA) $(top_builddir)/Cantera/fortran/src/cantera_kinetics.mod $(top_builddir)/build/include/cantera_kinetics.mod
# 	$(INSTALL_DATA) $(top_builddir)/Cantera/fortran/src/cantera_funcs.mod $(top_builddir)/build/include/cantera_funcs.mod
# 	$(INSTALL_DATA) $(top_builddir)/Cantera/fortran/src/cantera_iface.mod $(top_builddir)/build/include/cantera_iface.mod
# 	$(INSTALL_DATA) $(top_builddir)/Cantera/fortran/src/cantera_thermo.mod $(top_builddir)/build/include/cantera_thermo.mod
# 	$(INSTALL_DATA) $(top_builddir)/Cantera/fortran/src/cantera_transport.mod $(top_builddir)/build/include/cantera_transport.mod
# 	$(INSTALL_DATA) $(top_builddir)/Cantera/fortran/src/cantera_xml.mod $(top_builddir)/build/include/cantera_xml.mod
# 	$(INSTALL_DATA) $(top_builddir)/Cantera/fortran/src/fctxml.mod $(top_builddir)/build/include/fctxml.mod

CLEANFILES = *.o

# header file accumulation
all-local:
	@echo copying headers $<
	cp -p *.h $(top_builddir)/build/include/
