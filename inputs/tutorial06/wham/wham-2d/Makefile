#Linux flags for intel compiler
#CC=icc
#CFLAGS= -fast -Wall
#LFLAGS= -lm

# Linux flags
CC=gcc
CFLAGS= -O2 -Wall -W
#Debugging flags
#CFLAGS= -g  -Wall -W
LFLAGS= -lm

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

SRCS=wham-2d.c file_read.c histogram.c $(NR)/ran2.c $(NR)/locate.c ../wham/bootstrap.c
OBJS=wham-2d.o file_read.o histogram.o $(NR)/ran2.o $(NR)/locate.o ../wham/bootstrap.o

wham-2d: $(OBJS)
	$(CC) $(CFLAGS) -o wham-2d $(OBJS) $(LFLAGS)

depend:
	makedepend -- $(CFLAGS) -- $(SRCS)

clean:
	rm -f wham-2d ${OBJS} .*.swp *~

remake:
	make clean
	make depend
	make

tags: ${SRCS} *.h
	ctags ${SRCS} *.h

get:
	rsync -avuz  --exclude '.*.swp' /home/agrossf/gsa/wham/wham-2d/ .
put:
	rsync -avuz  --exclude '.*.swp' ./ /home/agrossf/gsa/wham/wham-2d/

sync: get put
# DO NOT DELETE
