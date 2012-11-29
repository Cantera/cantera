#!/bin/sh
###############################################################
#  $Author: hkmoffa $
#  $Date: 2010-03-19 18:26:02 -0600 (Fri, 19 Mar 2010) $
#  $Revision: 420 $
###############################################################
#

all: 
	cd StoichSolidKinetics;         gmake all
	cd surfkin_dupl_test;         gmake all

depends:
	cd StoichSolidKinetics ;         gmake depends
	cd surfkin_dupl_test;         gmake depends

clean:
	cd StoichSolidKinetics ;         gmake clean
	cd surfkin_dupl_test;         gmake clean

test:
	cd StoichSolidKinetics ;         gmake -s test
	cd surfkin_dupl_test;         gmake  -s test

rebaseline:
	cd StoichSolidKinetics ;         gmake rebaseline
	cd surfkin_dupl_test;         gmake rebaseline



