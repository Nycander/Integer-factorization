main: factoring.c pollard.c
	gcc -o factoring.exe pollard.c factoring.c -O2 -std=gnu99 -lgmp -lm -Wall -g
