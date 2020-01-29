CC = mpicc

SRC_DIR = ../../src
INCLUDE_DIR = -I${SRC_DIR}

UTIL = ${SRC_DIR}/array_util.c

.PHONY: demo0 demo1 demo2

demo0: ${SRC} demo0.c
	${CC} -o mpi_demo0  ${INCLUDE_DIR} ${UTIL} demo0.c 
 
demo1: ${SRC} demo1.c
	${CC} -o mpi_demo1  ${INCLUDE_DIR} ${UTIL} demo1.c 

demo2: ${SRC} demo2.c
	${CC} -o mpi_demo2  ${INCLUDE_DIR} ${UTIL} demo2.c
