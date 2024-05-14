COMPILER = gcc
FLAGS = -Wall -O3 -mtune=native

LIB = -lgsl -lm

sir:
	@${COMPILER} ${FLAGS} sir.c ${LIB} -o sir.out

clean:
	@rm -f *.out
	@rm -f dat/*.dat
	@rm -f plt/*.png
	@rm -f plt/*.pdf
	@rm -f plt/*.mp4
