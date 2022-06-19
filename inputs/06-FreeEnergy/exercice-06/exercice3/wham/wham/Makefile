# Linux flags
CC=gcc
CFLAGS= -O2 -Wall -W 
#Debugging flags
#CFLAGS= -g  -Wall -W
LFLAGS= -lm

#Linux flags for intel compiler
#CC=icc
#CFLAGS= -fast -Wall
#LFLAGS= -lm

# SGI
#CC=cc
#CFLAGS = -O2 -mips4
#LFLAGS= -lm

# Alpha
#CC=cc
#CFLAGS = -fast -arch host -tune host
#LFLAGS = -fast -non_shared -om -WL,-om_no_inst_sched -lm

# location of Numerical Recipes code
NR=../nr

SRCS=wham.c file_read.c histogram.c bootstrap.c $(NR)/ran2.c $(NR)/locate.c
OBJS=wham.o file_read.o histogram.o bootstrap.o $(NR)/ran2.o $(NR)/locate.o
HEADER_FILES= wham.h

#all: wham driver
all: wham

wham: $(OBJS)
	$(CC) $(CFLAGS) -o wham $(OBJS) $(LFLAGS)

#driver: bootstrap.o driver.o $(NR)/ran2.o $(NR)/locate.o
#	$(CC) $(CFLAGS) -o driver bootstrap.o driver.o ran2.o locate.o $(LFLAGS)

depend:
	makedepend -- $(CFLAGS) -- $(SRCS)

clean:
	rm -f wham driver *.o .*.swp *~

remake:
	make clean
	make depend
	make wham

tags: *.c
	ctags *.c *.h

dist:
	. dist_command

get:
	rsync -avuz  --exclude '.*.swp' --exclude 'data' --exclude 'test' /home/agrossf/gsa/wham/ .
put:
	rsync -avuz  --exclude '.*.swp' --exclude 'data' --exclude 'test' ./ /home/agrossf/gsa/wham/

sync: get put
# DO NOT DELETE
