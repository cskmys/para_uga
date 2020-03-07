all: exo1 exo2 exo3

exo1: exo1.o
	gcc -fopenmp -o exo1 exo1.o

exo1.o: exo1.c
	gcc -c -O0 -fopenmp exo1.c

exo2: exo2.o
	gcc -fopenmp -o exo2 exo2.o

exo2.o: exo2.c
	gcc -c -O2 -fopenmp exo2.c

exo3: exo3.o
	gcc -fopenmp -o exo3 exo3.o

exo3.o: exo3.c
	gcc -c -O2 -fopenmp exo3.c

clean:
	rm -f exo1 exo2 exo3 exo1.o exo2.o exo3.o *~
