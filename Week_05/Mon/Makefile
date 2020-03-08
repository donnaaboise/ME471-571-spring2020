CC = mpicc

SRC_DIR = ../../src
INCLUDE_DIR = -I${SRC_DIR}

UTIL = ${SRC_DIR}/array_util.c

.PHONY: demo0 demo1 demo2 demo3 demo4

demo0: ${SRC} demo0.c
	${CC} -o demo0  ${INCLUDE_DIR} ${UTIL} demo0.c
 
demo1: ${SRC} demo1.c
	${CC} -o demo1  ${INCLUDE_DIR} ${UTIL} demo1.c 

demo2: ${SRC} demo2.c
	${CC} -o demo2  ${INCLUDE_DIR} ${UTIL} demo2.c

demo3: ${SRC} demo3.c
	${CC} -o demo3  ${INCLUDE_DIR} ${UTIL} demo3.c

demo4: ${SRC} demo4.c
	${CC} -o demo4  ${INCLUDE_DIR} ${UTIL} demo4.c

demo5: ${SRC} demo5.c
	${CC} -o demo5  ${INCLUDE_DIR} ${UTIL} demo5.c

integral: ${SRC} integral.c
	${CC} -o integral  ${INCLUDE_DIR} ${UTIL} integral.c -lm
