CC = mpicc

SRC_DIR = ../src
INCLUDE_DIR = -I${SRC_DIR}

UTIL = ${SRC_DIR}/array_util.c

.PHONY: integral

integral: ${SRC} integral.c
	${CC} -o integral  ${INCLUDE_DIR} ${UTIL} integral.c -lm
