rps:
	@gcc -O3 -Wall src/rps.c -lgsl -lgslcblas -lm -o rps.out

avg:
	@gcc -O3 -Wall src/avg-sd.c -lgsl -lgslcblas -lm -o avg-sd.out

fft:
	@gcc -O3 -Wall src/fft.c -lfftw3 -lm -o fft.out

f_p:
	@gcc -O3 -Wall src/f_p.c -lgsl -lgslcblas -lm -o f_p.out

clean:
	@rm -f *.out
	@rm -f dat*/*.dat

clean-fig:
	@rm -f plt*/*.png
	@rm -f plt*/*.pdf

