CC=g++
CFLAGS=-I.
main: main.o neuron.o
	$(CC) -o main main.o neuron.o -I.

clean:
	rm -f *.o