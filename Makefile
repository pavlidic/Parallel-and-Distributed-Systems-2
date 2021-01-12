all:
	mpicc -O3 kNN.c lib.c VPT.c insert.c -o kNN -lopenblas -lm

clean:
	rm kNN
