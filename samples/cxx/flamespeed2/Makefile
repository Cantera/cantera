include /usr/local/include/cantera/Cantera.mak

CC=gcc
CXX=g++
RM=rm -f
CXXFLAGS=-pthread -O3 -Wno-inline -g -std=c++0x
CPPFLAGS=$(CANTERA_INCLUDES)
LDFLAGS=
LDLIBS=$(CANTERA_LIBS)

SRCS=flamespeed.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: flamespeed

clean:
	$(RM) $(OBJS) flamespeed

dist-clean: clean
	$(RM) *~ 
