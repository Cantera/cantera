include $(top_srcdir)/doxygen/aminclude.am

# Build in these directories:
SUBDIRS  = ext src test_problems
SUBDIRS += doxygen

# Distribute these directories:
DISTDIRS = src

# generated documentation will be included in distributed archive
EXTRA_DIST  = bootstrap $(DX_CONFIG)
EXTRA_DIST += AUTHORS README

# Support for pkgconfig
pkgconfigdir   = $(libdir)/pkgconfig
pkgconfig_DATA = cantera.pc
