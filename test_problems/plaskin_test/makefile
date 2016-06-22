include $(HOME)/Codes/plasmacantera-0.0.1/mycantera/include/cantera/Cantera.mak

CXX    = g++
FLAG   = #-g -O0
INCS   = $(CANTERA_INCLUDES)
LIBS   = $(CANTERA_LIBS)

all: exec

exec: main.o
	$(CXX) -o exec main.o $(LIBS)

main.o: main.cpp
	$(CXX) -c $(FLAG) $< $(INCS)

clean:
	rm *.o

dist-clean:
	rm *~
