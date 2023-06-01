
CC=g++

all: flexiplex 

flexiplex: flexiplex.c++ 
	${CC} $^ -Ofast -pthread -std=c++17 -o $@ edlib-1.2.7/edlib.cpp -I edlib-1.2.7/

clean:
	rm flexiplex 
