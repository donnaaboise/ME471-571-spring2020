CC = mpicc

SRC = ../../src
INCLUDE_DIR = -I../../src 

objs = ${SRC}/array_util.c

.PHONY: demo1
 
demo1: demo1.c ${objs}
	${CC} -o mpi_demo1  ${INCLUDE_DIR} demo1.c ${objs}  
