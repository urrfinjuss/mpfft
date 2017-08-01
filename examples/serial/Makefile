CC=gcc
FLAGS= -std=gnu99 -I/home/orange/usr/include -L/home/orange/usr/lib -Ofast -flto
LFLAGS=-lmfft -lmpfr -lpthread -lm

SRCS=main.c
OBJS=$($SRCS:.c=.o)
TARGET=bench.x

bench: main.c
	$(CC) $(FLAGS) -o $(TARGET) $(SRCS) $(LFLAGS)
