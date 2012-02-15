# Distribute these directories:
EXTRA_DIST = include

dist-hook:
	rm -rf `find $(distdir)/include -type d -name .svn`
	rm -rf `find $(distdir)/include -name *.h`