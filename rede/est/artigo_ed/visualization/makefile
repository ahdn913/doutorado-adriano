COMPILER = gcc
FLAGS = -Wall -O3 -mtune=native

LIB = -lgsl -lgslcblas -lm

triang:
	@${COMPILER} ${FLAGS} src/triang.c 		${LIB} -o rede-triang.out

honey:
	@${COMPILER} ${FLAGS} src/honeycomb.c ${LIB} -o rede-honeycomb.out

kagome:
	@${COMPILER} ${FLAGS} src/kagome.c 		${LIB} -o rede-kagome.out

vor:
	@${COMPILER} ${FLAGS} src/vor.c ${LIB} -o vor.out

clean:
	@rm plt/*.png
	@rm dat/*.dat
