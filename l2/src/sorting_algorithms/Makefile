CC = gcc
GIVEN_FLAGS = -O2
PEDANTIC_PARANOID_FREAK = -g3 -O0 -Wall -Wshadow -Wcast-align \
				-Waggregate-return -Wstrict-prototypes \
				-Wredundant-decls -Wnested-externs \
				-Wpointer-arith -Wwrite-strings -finline-functions -Werror

CFLAGS = -fopenmp $(PEDANTIC_PARANOID_FREAK)
LDFLAGS = -fopenmp

EXEC = 	bubble.run	\
	qsort.run	\
	mergesort.run


HEADER_FILES = $(wildcard *.h)

RAND_INIT=0

ifdef p
$(info Test focus on parallel)
CONFIG_FLAGS += -DTEST_ONLY_PARA
endif

ifdef v
$(info Print intermediate results)
CONFIG_FLAGS += -DPRINT_INT
endif

ifdef r
$(info Initialization of the vector is random)
CONFIG_FLAGS += -DRINIT
else
RAND_INIT=0
ifeq ($(RAND_INIT), 1)
$(info Initialization of the vector is random)
CONFIG_FLAGS += -DRINIT
endif
endif


all: $(EXEC)

%.run: %.o utils.o
	$(CC) $(LDFLAGS) -o $@ $^

%.o: %.c $(HEADER_FILES)
	$(CC) -c $(CONFIG_FLAGS) $(CFLAGS) $< -o $@

clean:
	rm -f $(EXEC) *.o *~

.PHONY: clean
