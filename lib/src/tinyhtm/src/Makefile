CC      =	gcc
CFLAGS  =	-fPIC -m64 -g -c -I../include -std=gnu99
LIBNAME =	libhtm.a

SRCS =	common.c     \
			geometry.c   \
			htm.c        \
			id_list.c    \
			select.c     \
			tree.c       \
			tree_count.c \
			tree_gen.c

OBJS    =	${SRCS:.c=.o}

libhtm.a:	${OBJS}
		rm -f ${LIBNAME}
		ar q ${LIBNAME} ${OBJS}

clean:
		rm -f libhtm.a *.o *.a a.out core
