CC=g++
CFLAGS=-std=c++11 -larmadillo -O3
DEPS = jacobi.h
OBJ = jacobi.o tests-project.o

%.o: %.cpp $(DEPS)
        $(CC) -c -o $@ $< $(CFLAGS)

test: jacobi.cpp
	$(CC) -o jacobi jacobi.cpp $(CFLAGS)

