CC=gcc
CFLAGS=-I. -lm -ggdb
OBJ = tensor.o

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

tensor: $(OBJ)
	gcc -o $@ $^ $(CFLAGS)

clean:
	rm -f $(OBJ) tensor
