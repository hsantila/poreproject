CC=gcc
CFLAGS="-Wall" "-std=c99"

debug:clean
	$(CC) $(CFLAGS) -g -o pore main.c io.c mymath.c tools.c -lm
stable:clean
	$(CC) $(CFLAGS) -o pore main.c io.c mymath.c tools.c -lm
clean:
	rm -vfr *~ pore
