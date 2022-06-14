

CC=g++

#GPROF=/opt/homebrew/Cellar/gperftools/2.9.1_1

all: flexiplex 

flexiplex: flexiplex.c++ 
	${CC} $^ -O3 -std=c++11 -o $@ edlib-1.2.7/edlib.cpp -I edlib-1.2.7/

clean:
	rm flexiplex 
