all:
		(cd src; make)

install:
		(cd src; install sptQueryLib.o ..)
		(cd src; install sptQuery.h ../include)
		(cd src/tinyhtm/src; install libhtm.a ../../..)
		(cd src/tinyhtm/include; install -d tinyhtm ../../../include/tinyhtm)
		(cd src/tinyhtm/include; install tinyhtm/* ../../../include/tinyhtm)
		(cd src/tinyhtm/include; install tinyhtm.h ../../../include)

clean:
		rm -f *.a *.o include/*.h
		(cd src; make clean)
