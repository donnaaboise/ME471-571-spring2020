CC = gcc

heat1d: heat1d_serial.c
	${CC} -o heat1d  heat1d_serial.c -lm

heat2d: heat2d_serial.c
	${CC} -o heat2d  heat2d_serial.c -lm

spiral: spiral_serial.c
	${CC} -o spiral  spiral_serial.c -lm


clean:
	rm -rf heat1d heat2d spiral *.out

