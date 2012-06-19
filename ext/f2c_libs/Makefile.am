SUBDIRS = arithchk

cc_sources              = abort_.c c_log.c d_atn2.c dfe.c d_prod.c \
                          due.c exit_.c getenv_.c hl_le.c i_dim.c  \
                          inquire.c lread.c pow_ii.c r_acos.c      \
                          r_dim.c r_nint.c r_tan.c                 \
                          uninit.c z_abs.c close.c                 \
                          d_cnjg.c d_imag.c d_sign.c ef1asc_.c     \
                          f77_aloc.c h_abs.c hl_lt.c i_dnnt.c      \
                          i_sign.c lwrite.c r_asin.c               \
                          rewind.c rsfe.c r_tanh.c util.c \
                          z_cos.c backspac.c c_sin.c d_cos.c       \
                          d_int.c d_sin.c ef1cmc_.c f77vers.c      \
                          h_dim.c h_mod.c i_indx.c lbitbits.c      \
                          pow_ri.c r_atan.c r_exp.c r_sign.c       \
                          s_cat.c s_rnge.c wref.c z_div.c c_abs.c  \
	                  c_sqrt.c d_cosh.c d_lg10.c d_sinh.c      \
                          endfile.c fmt.c h_dnnt.c h_nint.c iio.c  \
                          lbitshft.c open.c pow_zi.c r_atn2.c      \
                          r_imag.c r_sin.c s_cmp.c s_stop.c        \
                          wrtfmt.c z_exp.c cabs.c d_abs.c d_dim.c  \
                          d_log.c d_sqrt.c fmtlib.c h_indx.c\
                          h_sign.c i_len.c l_ge.c pow_ci.c pow_zz.c\
                          r_cnjg.c r_int.c r_sinh.c s_copy.c sue.c \
                          wsfe.c z_log.c c_cos.c d_acos.c \
                          d_mod.c d_tan.c                   \
                          h_len.c i77vers.c ilnw.c l_gt.c pow_dd.c \
                          r_cos.c r_lg10.c rsli.c sfe.c \
                          system_.c wsle.c z_sin.c c_div.c d_asin.c\
                          d_nint.c d_tanh.c err.c ftell_.c\
                          hl_ge.c i_abs.c i_mod.c l_le.c pow_di.c  \
                          r_cosh.c r_log.c rsne.c                  \
                          sig_die.c typesize.c wsne.c z_sqrt.c     \
                          c_exp.c d_atan.c d_exp.c dolio.c dtime_.c\
                          etime_.c hl_gt.c                         \
                          i_nint.c l_lt.c pow_hh.c r_abs.c rdfmt.c \
                          r_mod.c r_sqrt.c signal_.c uio.c xwsne.c 

h_sources               = 

AM_CPPFLAGS             = -I.
AM_CXXFLAGS             = $(AM_CPPFLAGS)
AM_FCFLAGS              = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libctf2c.la
library_includedir      = $(top_builddir)/include
library_include_HEADERS = $(h_sources)

#-----------------------
# Cantera Converters C/C++ library
#-----------------------

__top_builddir__build_lib_libctf2c_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libctf2c_la_SOURCES      = $(cc_sources) $(h_sources)

__top_builddir__build_lib_libctf2c_include_HEADERS = $(h_sources)
__top_builddir__build_lib_libctf2c_includedir = $(prefix)/include


# header file accumulation
all:
	@echo copying headers $<
	cp -p sysdep1.h0 sysdep1.h
	cp -p signal1.h0 signal1.h

CLEANFILES = *.o
