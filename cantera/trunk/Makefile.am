include $(top_srcdir)/doc/doxygen/aminclude.am

# Build in these directories:
SUBDIRS  = ext src test_problems
SUBDIRS += doc

# Distribute these directories:
DISTDIRS = ext src test_problems doc include

# generated documentation will be included in distributed archive
EXTRA_DIST  = bootstrap $(DX_CONFIG)
EXTRA_DIST += AUTHORS README

# Support for pkgconfig
pkgconfigdir   = $(libdir)/pkgconfig
pkgconfig_DATA = cantera.pc
