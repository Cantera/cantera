
SRCS = kinetics1.cpp

OBJS = $(SRCS:.cpp=.o)
EXES = $(SRCS:.cpp=.x)
MKS = $(SRCS:.cpp=.mak)

all: $(EXES)

%.mak:
	ctnew; sed 's/demo/$*/g' demo.mak > $*.mak

%.x:
	make $*.mak; make -f $*.mak; mv $* $*.x
	@echo type '$*.x' to run the program

clean:
	rm -f $(OBJS) $(EXES) $(MKS)


