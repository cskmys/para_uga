CC = gcc
GIVEN_FLAGS = -O2
PEDANTIC_PARANOID_FREAK = -g3 -O0 -Wall -Wshadow -Wcast-align \
				-Waggregate-return -Wstrict-prototypes \
				-Wredundant-decls -Wnested-externs \
				-Wpointer-arith -Wwrite-strings -finline-functions -Werror

CFLAGS = -fopenmp $(PEDANTIC_PARANOID_FREAK)
LDFLAGS = -fopenmp

EXEC = 	triangular_matrix.run

all: $(EXEC)

%.run: %.o
	$(CC) $(LDFLAGS) -o $@ $^

%.o: %.c
	$(CC) -c $(CONFIG_FLAGS) $(CFLAGS) $< -o $@

clean:
	rm -f $(EXEC) *.o *~

.PHONY: clean
