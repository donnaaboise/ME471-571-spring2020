CC = gcc
NVC = nvcc

me471 = ${HOME}/ME471-571-master

SRC = ${me471}/src
INCLUDE_DIR = -I${me471}/src -I/usr/include

objs = ${SRC}/demo_util.c

.PHONY: hello_demo add_simple add count_devices parallel_simple2 pool

hello_demo: hello_demo.cu ${objs}
	${NVC} -o hello_demo  ${INCLUDE_DIR} hello_demo.cu ${objs} -lm

add_simple: add.cu ${objs}
	${NVC} -o add_simple  ${INCLUDE_DIR} add_simple.cu ${objs} -lm

add: add.cu ${objs}
	${NVC} -o add  ${INCLUDE_DIR} add.cu ${objs} -lm

count_devices: count_devices.cu ${objs}
	${NVC} -o count_devices  ${INCLUDE_DIR} count_devices.cu ${objs} -lm

parallel_simple: parallel_simple.cu ${objs}
	${NVC} -o parallel_simple  ${INCLUDE_DIR} parallel_simple.cu ${objs} -lm

parallel_simple2: parallel_simple2.cu ${objs}
	${NVC} -o parallel_simple2  ${INCLUDE_DIR} parallel_simple2.cu ${objs} -lm

pool: pool.cu ${objs}
	${NVC} -o pool  ${INCLUDE_DIR} pool.cu ${objs} -lm
