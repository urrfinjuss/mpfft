CC=gcc
FLAGS=-std=gnu99 -I/home/orange/usr/include -L/home/orange/usr/lib -Ofast -flto
LFLAGS=-lmpfft -lmpfr -lpthread -lfftw3 -lm

SRCS1=main.c
SRCS2=one_fft.c
OBJS1=$($SRCS1:.c=.o)
OBJS2=$($SRCS2:.c=.o)
TARGET1=bench.x
TARGET2=one_fft.x

.PHONY: all
all: $(TARGET1) $(TARGET2)

bench.x: main.c
	$(CC) $(FLAGS) -o $(TARGET1) $(SRCS1) $(LFLAGS)

one_fft.x: one_fft.c
	$(CC) $(FLAGS) -o $(TARGET2) $(SRCS2) $(LFLAGS)
clean:
	rm -f $(TARGET1) $(TARGET2)
