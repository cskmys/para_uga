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

EXTRA_SRC = utils.o ring.o matrix_op.o


CHECK_CORRECT=1

EVAL_PERF=1

ifeq ($(EVAL_PERF), 0)
ifeq ($(CHECK_CORRECT), 0)
$(error Wrong combination of options. EVAL_PERF or CHECK_CORRECT need to be set)
endif
endif

ifeq ($(CHECK_CORRECT), 1)
$(info Computation if checked for correctness)
CONFIG_FLAGS += -DCHECK_CORRECTNESS
endif

ifeq ($(EVAL_PERF), 1)
$(info Performance is evaluated)
CONFIG_FLAGS += -DPERF_EVAL
endif


all: $(EXEC)

%.run: %.o $(EXTRA_SRC)
	$(CC) $(LDFLAGS) -o $@ $^

%.o: %.c $(HEADER_FILES)
	$(CC) -c $(CONFIG_FLAGS) $(CFLAGS) $< -o $@

clean:
	rm -f $(EXEC) *.o *~

.PHONY: clean