CC = mpicc
GIVEN_FLAGS = -O2
PEDANTIC_PARANOID_FREAK = -g3 -O0 -Wall -Wshadow -Wcast-align \
				-Waggregate-return -Wstrict-prototypes \
				-Wredundant-decls -Wnested-externs \
				-Wpointer-arith -Wwrite-strings -finline-functions -Werror

CFLAGS = $(PEDANTIC_PARANOID_FREAK)
LDFLAGS =

EXEC = 	main.run

HEADER_FILES = $(wildcard *.h)

all: $(EXEC)

%.run: %.o ring.o
	$(CC) $(LDFLAGS) -o $@ $^

%.o: %.c $(HEADER_FILES)
	$(CC) -c $(CONFIG_FLAGS) $(CFLAGS) $< -o $@

clean:
	rm -f $(EXEC) *.o *~

.PHONY: clean