CC=gcc
RM=rm -f
LIB_PATH=$(HOME)/usr/lib
INC_PATH=$(HOME)/usr/include
FLAGS=-std=gnu99 -march=native -Ofast -I/home/orange/code/mpfft/src/ -I$(INC_PATH)
CFLAGS=-std=gnu99 -fPIC -Wall -Wextra -g -I/home/orange/code/mpfft/src/ -I$(INC_PATH)
LDFLAGS=-shared

TARGET=libmpfft.so
TARGET_SERIAL=libmpfft.so
TARGET_THREAD=libmpfft_threads.so

SOURCES=$(wildcard src/*.c)
SOURCES_COMMON=src/mpfft_bitrev.c src/mpfft_init.c src/mpfft_set.c src/dfft_bitrev.c src/dfft_algorithm.c src/mpfft_version.c src/mpfr_version.c 
SOURCES_SERIAL=$(SOURCES_COMMON) src/mpfft_serial.c
SOURCES_THREAD=$(SOURCES_COMMON) src/mpfft_pthread.c

HEADERS=$(wildcard src/*.h)
HEADERS_COMMON=src/header.h src/dfft_header.h src/mpfft_header.h
HEADERS_SERIAL=$(HEADERS_COMMON) src/mpfft_serial.h
HEADERS_THREAD=$(HEADERS_COMMON) src/mpfft_pthread.h

OBJECTS=$(SOURCES:.c=.o)
OBJECTS_SERIAL=$(SOURCES_SERIAL:.c=.o)
OBJECTS_THREAD=$(SOURCES_THREAD:.c=.o)

.PHONY: all
all: ${TARGET_SERIAL} ${TARGET_THREAD}

$(TARGET_SERIAL): $(OBJECTS_SERIAL)
	$(CC) ${FLAGS} ${LDFLAGS} -o $@ $^
	cp $(TARGET_SERIAL) $(LIB_PATH)
	cp $(HEADERS_SERIAL) $(INC_PATH)

$(TARGET_THREAD): $(OBJECTS_THREAD)
	$(CC) ${FLAGS} ${LDFLAGS} -o $@ $^
	cp $(TARGET_THREAD) $(LIB_PATH)
	cp $(HEADERS_THREAD) $(INC_PATH)

$(SOURCES:.c=.d):%.d:%.c
	$(CC) $(CFLAGS) -MM $< >$@
include $(SOURCES:.c=.d)

.PHONY: clean
clean:
	-${RM} ${TARGET_SERIAL} ${OBJECTS_SERIAL} $(SOURCES_SERIAL:.c=.d)
	-${RM} ${TARGET_THREAD} ${OBJECTS_THREAD} $(SOURCES_THREAD:.c=.d)


