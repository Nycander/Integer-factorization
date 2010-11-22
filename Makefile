CC = gcc
CFLAGS = -O2 -std=gnu99 -lgmp -lm -Wall -g

#Make will build the first target in the file
main.exe: primes.c factor_list.c pollard.c trialdivision.c qs.c factoring.c
	$(CC) primes.c factor_list.c pollard.c trialdivision.c qs.c factoring.c $(CFLAGS) -o factoring.exe

generator: prime_generator.c
	$(CC) prime_generator.c $(CFLAGS) -o generator.exe
