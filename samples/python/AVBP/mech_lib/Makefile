customkinetics.so: customkinetics.f90
	gfortran -ffixed-line-length-0 -c customkinetics.f90 -g -fPIC -o customkinetics.o
	gfortran -shared -o customkinetics.so customkinetics.o            
    