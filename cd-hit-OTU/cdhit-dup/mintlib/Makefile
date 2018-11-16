
CC = g++

CFLAGS = -ggdb -O2 -Wall -Wno-unused -I.
LFLAGS = -fPIC


UNAME = $(shell uname)

ifeq ($(UNAME), Linux)
	CFLAGS += -DUNIX
endif

ifeq ($(UNAME), Darwin)
	CFLAGS += -DUNIX -DMAC_OSX
endif


OBJECTS = minString.o minMap.o test.o

first: all

all: test

.SUFFIXES: .c .obj .cpp .cc .cxx .C

.cxx.o:
	$(CC) -c $(CFLAGS) -o $@ $<

test: $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o test

clean:
	rm *.o
