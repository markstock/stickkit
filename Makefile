OPTS=-O2
#OPTS=-g -ggdb

EXE=stickkit
OBJS=stickkit.o svgout.o png.o bob.o
#CFLAGS=-I/usr/include/wxSVG
CFLAGS=-fPIC -I/usr/include/wxSVG `wx-config --cflags`
LDFLAGS=-fPIC `wx-config --libs`

all : $(EXE)

bob.o : stickkit.h bob.c
	cc $(OPTS) -std=c99 $(CFLAGS) -c bob.c

png.o : stickkit.h png.c
	cc $(OPTS) -std=c99 $(CFLAGS) -c png.c

svgout.o : stickkit.h svgout.cpp
	g++ $(OPTS) $(CFLAGS) -c svgout.cpp

stickkit.o : stickkit.h stickkit.c
	cc $(OPTS) -std=c99 $(CFLAGS) -c stickkit.c

$(EXE) : $(OBJS)
	g++ $(OPTS) -o $(EXE) $(OBJS) $(LDFLAGS) -lm -lwxsvg -lpng

clean :
	rm $(EXE) $(OBJS)
