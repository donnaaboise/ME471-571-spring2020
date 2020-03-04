CC = mpicc


deriv: ${SRC} deriv_serial.c
	${CC} -o deriv  ${INCLUDE_DIR} ${UTIL} deriv_serial.c

derivp: ${SRC} deriv_parallel.c
	${CC} -o derivp  ${INCLUDE_DIR} ${UTIL} deriv_parallel.c
