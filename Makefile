###############################
### DO NOT MODIFY THIS FILE ###
###############################

CC = clang
# CFLAGS = -g -O0 -W -Wall -pedantic -std=c11
LDFLAGS = -lm -lpthread

.PHONY: all clean haste test rain update submission

all: scaler
	clang -g -O0 -W -Wall -pedantic -std=c11   -c -o vector.o vector.c
	clang -lm -lpthread  scaler.o vector.o   -o scaler

scaler: scaler.o vector.o

clean:
	-rm -f *.o
	-rm -f scaler vector

