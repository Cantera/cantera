                             C A N T E R A
                   E X T E N D E D   T E S T   S U I T E 

                              release 1.8

                               03/2010

***************************************************************
***   Test suite including unit tests for the C++ interface
***   to Cantera.
***************************************************************


License information
===================

See the file "License.txt" for information on the terms & conditions
for usage, and a DISCLAIMER OF ALL WARRANTIES.

All trademarks referenced herein are property of their respective
holders.


To Use the Test Suite
=====================

The extended test suite assumes that there exists a working Cantera
installation and utilizes the installation directories of that installation
exclusively. 

The base location of the installation must be specified in the preconfig file
using the environmental variable CANTERA_INSTALL_DIR.

To configure the test suite, either edit the preconfig file. 
or alternatively environment variables may be set before calling the preconfig
file.
 
1)  Run autoconf to create the file configure from the file configure.in.

1)  Edit the preconfig file (or call preconfig with edited environmental 
    variables).

2)  Run the preconfig file which will run the file configure. This will
    create Makefiles in all of the test directories.

3)  Type 'make test' to run all of the tests.

Some of the tests still use asc diff file comparisons within the test.
For these tests, small changes in the asc output may indicate falsely
that the tests failed, while in actuality there are only smal differences
in numerical round-off.  Therefore, the output of these tests should
be manually examined to determine whether the tests failed significantly.


Web sites
==========

The current main Cantera web site is under construction. There exists three
web sites for communication amongst Cantera users:

1. The Cantera User's Group.
    http://groups.google.com/group/cantera-users
   This site has a message board, and some miscellaneous files and utilities.

2. The Cantera googleCode site.  Distribution of the Cantera source code is
   done using googleCode. The site is http://code.google.com/p/cantera/.

3. The Cantera developers site .
   http://groups.google.com/group/cantera-dev
   Limited access site where developers can discuss development ideas.

=============================================================

*  Old Cantera source code site. Distribution of the Cantera source code was
   done using SourceForge. The old site is located at
   http://sourceforge.net/projects/cantera.
   It still contains the project cvs history.

