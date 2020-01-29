CC = mpicc

SRC_DIR = ../../src
INCLUDE_DIR = -I${SRC_DIR}

UTIL = ${SRC_DIR}/array_util.c

.PHONY: demo1 demo2
 
demo1: ${SRC} demo1.c
	${CC} -o mpi_demo1  ${INCLUDE_DIR} ${UTIL} demo1.c 

demo2: ${SRC} demo2.c
	${CC} -o mpi_demo2  ${INCLUDE_DIR} demo2.c
