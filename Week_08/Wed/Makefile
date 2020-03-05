CC = mpicc


deriv: ${SRC} deriv_serial.c
	${CC} -o deriv  ${INCLUDE_DIR} ${UTIL} deriv_serial.c

conn: ${SRC} connections.c
	${CC} -o conn  ${INCLUDE_DIR} ${UTIL} connections.c

derivp: ${SRC} deriv_parallel.c
	${CC} -o derivp  ${INCLUDE_DIR} ${UTIL} deriv_parallel.c
