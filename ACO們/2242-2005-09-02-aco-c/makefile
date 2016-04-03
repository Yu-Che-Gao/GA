# makefile for compiling ACO-C source code
CC=gcc
LD=gcc
RM=rm -fv
CFLAGS=-g
OBJS:=$(patsubst %.c,%.o,$(wildcard *.c))
MAIN=aco
all:$(MAIN)
$(MAIN):$(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(MAIN) -lm
%.o: %.c dec.h func.h  par.h parameters.h
	$(CC) $(CFLAGS) -c $<
clean:
	$(RM) $(OBJS)

