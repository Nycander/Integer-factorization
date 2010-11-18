CC = gcc
CFLAGS = -O2 -std=gnu99 -lgmp -lm -Wall -g

#Make will build the first target in the file
main.exe: factoring.c factor.c pollard.c
	$(CC) factor.c pollard.c factoring.c $(CFLAGS) -o factoring.exe

