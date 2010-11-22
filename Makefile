CC = gcc
CFLAGS = -O2 -std=gnu99 -lgmp -lm -Wall -g

#Make will build the first target in the file
<<<<<<< HEAD
main.exe: factor_list.c pollard.c trialdivision.c qs.c factoring.c
	$(CC) factor_list.c pollard.c trialdivision.c qs.c factoring.c $(CFLAGS) -o factoring.exe
=======

main.exe: primes.c factor_list.c pollard.c trialdivision.c qs.c factoring.c
	$(CC) primes.c factor_list.c pollard.c trialdivision.c qs.c factoring.c $(CFLAGS) -o factoring.exe
>>>>>>> 385ac5fc431695b59058c77c801ccec2ecf3937a
