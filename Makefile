OPTS=-O2
#OPTS=-g -ggdb

EXE=stickkit
OBJS=stickkit.o svgout.o png.o bob.o

all : $(EXE)

bob.o : stickkit.h bob.c
	cc $(OPTS) -std=c99 -c bob.c

png.o : stickkit.h png.c
	cc $(OPTS) -std=c99 -c png.c

svgout.o : stickkit.h svgout.cpp
	g++ $(OPTS) -I/usr/include/wxSVG `wx-config --cxxflags` -c svgout.cpp

stickkit.o : stickkit.h stickkit.c
	cc $(OPTS) -std=c99 -c stickkit.c

$(EXE) : $(OBJS)
	g++ $(OPTS) -o $(EXE) $(OBJS) `wx-config --libs` -lm -lwxsvg -lpng

clean :
	rm $(EXE) $(OBJS)
