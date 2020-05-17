CC = mpicc
PEDANTIC_PARANOID_FREAK = -g3 -O0 -Wall -Wshadow -Wcast-align \
				-Waggregate-return -Wstrict-prototypes \
				-Wredundant-decls -Wnested-externs \
				-Wpointer-arith -Wwrite-strings -finline-functions -Werror

CFLAGS = $(PEDANTIC_PARANOID_FREAK)
LDFLAGS = -lm

EXEC = 	main.run

HEADER_FILES = $(wildcard *.h)

EXTRA_SRC = cart.o fox.o matrix_op.o mpi_comm.o ring.o utils.o


CHECK_CORRECT=1

EVAL_PERF=1

DBG_ON_FD=0

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

ifeq ($(DBG_ON_FD), 1)
$(info mpi_printfs on fd)
CONFIG_FLAGS += -DDBG_ON_FD
endif

all: $(EXEC)

%.run: %.o $(EXTRA_SRC)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c $(HEADER_FILES)
	$(CC) -c $(CONFIG_FLAGS) $(CFLAGS) $< -o $@

clean:
	rm -f $(EXEC) *.o *~

.PHONY: clean