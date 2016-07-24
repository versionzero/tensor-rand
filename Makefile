CC=gcc
CFLAGS=-I.
OBJ = tensor.o

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

tensor: $(OBJ)
	gcc -o $@ $^ $(CFLAGS)

clean:
	rm -f $(OBJ)
