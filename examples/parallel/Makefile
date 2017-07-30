CC=gcc
FLAGS= -I/home/sdyachen/usr/include -L/home/sdyachen/usr/lib -Ofast -flto
LFLAGS=-lmfft -lmpfr -lpthread -lm

SRCS_BENCH=bench.c
SRCS_THREAD=threads.c
OBJS_BENCH = $($SRCS_BENCH:.c=.o)
OBJS_THREAD = $($SRCS_THREAD:.c=.o)
BENCH=bench.x
THREAD=threads.x

bench: bench.c
	$(CC) $(FLAGS) -o $(BENCH) $(SRCS_BENCH) $(LFLAGS)

threads: threads.c
	$(CC) $(FLAGS) -o $(THREAD) $(SRCS_THREAD) $(LFLAGS)

clean:
	rm -f threads.x
	rm -f bench.x
