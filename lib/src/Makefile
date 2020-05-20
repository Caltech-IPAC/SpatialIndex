all:        sptQueryLib.o libhtm

libhtm:
				(cd tinyhtm/src; make)

sptQueryLib.o:	sptQueryLib.c
				gcc -fPIC -m64 -c -Wall -I. -Itinyhtm/include sptQueryLib.c

install:
				install sptQueryLib.o ../objs
				install sptQuery.h ../include
				(cd tinyhtm/src; install *.o ../../objs)

clean:
				(cd tinyhtm/src; make clean)
				rm -f *.o
