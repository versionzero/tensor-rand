CC=gcc
CFLAGS=-I.
LDFLAGS=-lm -ggdb
OBJ = debug.o gundersen.o main.o matrix.o tensor.o vector.o

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

tensor: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LDFLAGS)

clean:
	rm -f $(OBJ) tensor
