TARGET=solve_Poisson

FLAGS=-O0 -g

SOURCE=\
	solve_Poisson2d_parallel.c \
	cg_parallel.c \
	bicg_parallel.c \
	splitting2d_parallel.c \
	mpi_debug.c

OBJS=\
	solve_Poisson2d_parallel.o \
	cg_parallel.o \
	bicg_parallel.o \
	splitting2d_parallel.o \
	mpi_debug.o

CC=mpicc

all: ${OBJS} comm.o
	${CC} -o ${TARGET} ${FLAGS} ${SOURCE} comm.o -lm


comm_dup: ${OBJS} comm_duplicate.o
	${CC} -o ${TARGET} ${FLAGS} ${SOURCE} comm_duplicate.o -lm

comm_nodup: ${OBJS} comm_no_duplicate.o
	${CC} -o ${TARGET} ${FLAGS} ${SOURCE} comm_no_duplicate.o -lm


clean:
	rm -rf *.o
