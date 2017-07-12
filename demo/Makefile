CC=gcc
FLAGS=-I/home/sdyachen/usr/include -L/home/sdyachen/usr/lib -Ofast -flto
LFLAGS=-lmfft -lmpfr

SRCS=main.c
OBJS=$($SRCS:.c=.o)
TARGET=bench.x

bench: main.c
	$(CC) $(FLAGS) -o $(TARGET) $(SRCS) $(LFLAGS)
