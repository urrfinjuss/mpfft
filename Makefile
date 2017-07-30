CC=gcc
RM=rm -f
FLAGS=-std=gnu99 -march=native -Ofast -I/home/sdyachen/Work/mfft/include/
CFLAGS=-std=gnu99 -fPIC -Wall -Wextra -g -I/home/sdyachen/Work/mfft/include/
LDFLAGS=-shared

LIB_PATH=$(HOME)/usr/lib
INC_PATH=$(HOME)/usr/include
TARGET=libmfft.so

SOURCES=$(wildcard src/*.c)
HEADERS=$(wildcard src/*.h)
OBJECTS=$(SOURCES:.c=.o)

.PHONY: all
all: ${TARGET}

$(TARGET): $(OBJECTS)
	$(CC) ${FLAGS} ${LDFLAGS} -o $@ $^
	cp $(TARGET) $(LIB_PATH)
	cp $(HEADERS) $(INC_PATH)

$(SOURCES:.c=.d):%.d:%.c
	$(CC) $(CFLAGS) -MM $< >$@
include $(SOURCES:.c=.d)

.PHONY: clean
clean:
	-${RM} ${TARGET} ${OBJECTS} $(SOURCES:.c=.d)


