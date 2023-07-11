CC := gfortran
CFLAGS := -O3 -fforce-addr -fcaller-saves -funroll-loops

SRCS := main.f90

all: main

main: main.f90
	$(CC) $(SRCS) -o main $(CFLAGS)
