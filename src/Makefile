all:        libs sptIndx sptQuery htmCell

libs:
				(cd ../lib; make; make install)

sptIndx:		sptIndx.o
				gcc -m64 -o sptIndx sptIndx.o ../lib/src/tinyhtm/src/libhtm.a -lm

sptIndx.o:	sptIndx.c
				gcc -m64 -c -Wall -I../lib/src/tinyhtm/include sptIndx.c

sptQuery:	sptQuery.o
				gcc -m64 -o sptQuery sptQuery.o ../lib/sptQueryLib.o ../lib/src/tinyhtm/src/libhtm.a -lm

sptQuery.o:	sptQuery.c
				gcc -m64 -c -Wall -I../lib/include -I../lib/src/tinyhtm/include sptQuery.c

htmCell:		htmCell.o
				gcc -m64 -o htmCell htmCell.o ../lib/src/tinyhtm/src/libhtm.a -lm

htmCell.o:	htmCell.c
				gcc -m64 -c -Wall -I../lib/src/tinyhtm/include htmCell.c

install:
	         mkdir -p ../bin
				cp sptIndx sptQuery htmCell ../bin

clean:
				rm -f sptIndx sptQuery htmCell *.o
