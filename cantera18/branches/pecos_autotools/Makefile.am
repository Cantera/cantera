include $(top_srcdir)/doxygen/aminclude.am

# Build in these directories:
SUBDIRS  = Cantera tools test_problems
SUBDIRS += doxygen examples build

# Distribute these directories:
DISTDIRS = Cantera examples build

# generated documentation will be included in distributed archive
EXTRA_DIST  = bootstrap $(DX_CONFIG) config.h.in cantera.pc
EXTRA_DIST += AUTHORS README

# Support for pkgconfig
pkgconfigdir   = $(libdir)/pkgconfig
pkgconfig_DATA = cantera.pc

