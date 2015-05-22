CC = clang
CFLAGS = -g -O0 -W -Wall -pedantic -std=c11
LDFLAGS = -lm -lpthread

build:
	clang -g -O0 -W -Wall -pedantic -std=c11   -c -o scaler.o scaler.c
	clang -g -O0 -W -Wall -pedantic -std=c11   -c -o vector.o vector.c
	clang -lm -lpthread  scaler.o vector.o   -o scalerCorrect

clean:
	-rm -f *.o
	-rm -f scaler vector
