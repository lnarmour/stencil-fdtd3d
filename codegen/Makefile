CFLAGS= -O3  -std=c99  -I/usr/include/malloc/
LIBRARIES=-lm
CC=gcc
OBJS = fdtd3d.o rampsig.o printme.o
all: plain

debug: CFLAGS =-DDEBUG -g -Wall -Wextra -std=c99 -I/usr/include/malloc/
debug: all
		
plain: $(OBJS)
	$(CC) fdtd3d-wrapper.c -o fdtd3d $(OBJS) $(CFLAGS) $(LIBRARIES)

fdtd3d.o : fdtd3d.c
	$(CC) fdtd3d.c -o fdtd3d.o $(CFLAGS) $(LIBRARIES) -c

rampsig.o : rampsig.c
	$(CC) rampsig.c -o rampsig.o $(CFLAGS) $(LIBRARIES) -c

printme.o : printme.c
	$(CC) printme.c -o printme.o $(CFLAGS) $(LIBRARIES) -c

clean:
	rm -rf *.o fdtd3d 
