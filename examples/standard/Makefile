CC=gcc
FLAGS=-I/home/sdyachen/usr/include -L/home/sdyachen/usr/lib -Ofast -flto
LFLAGS=-lmfft -lmpfr -lpthread -lm

SRCS1=main.c
SRCS2=one_fft.c
OBJS1=$($SRCS1:.c=.o)
OBJS2=$($SRCS2:.c=.o)
TARGET1=bench.x
TARGET2=one_fft.x

bench: main.c
	$(CC) $(FLAGS) -o $(TARGET1) $(SRCS1) $(LFLAGS)

one_fft: one_fft.c
	$(CC) $(FLAGS) -o $(TARGET2) $(SRCS2) $(LFLAGS)

