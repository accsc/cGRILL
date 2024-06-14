# Makefile for cGRILL

CC = gcc
CFLAGS = -O3 -w -I .


PROGS = cGRILL
RM      = /bin/rm -f

OBJ=main.o

 
all:    $(PROGS)

cGRILL:	$(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -lm -o cGRILL.exe

clean:
	/bin/rm -f *.o cGRILL.exe
