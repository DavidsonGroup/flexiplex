CC=g++

all: flexiplex filter

flexiplex: flexiplex.c++ 
	${CC} $^ -Ofast -pthread -std=c++11 -o $@ edlib-1.2.7/edlib.cpp -I edlib-1.2.7/

clean:
	rm flexiplex 

install:
	cp flexiplex /usr/local/bin/flexiplex
	python3 -m pip install scripts/

uninstall:
	rm /usr/local/bin/flexiplex
	python3 -m pip list 2>/dev/null | grep -Fq flexiplex-filter && python3 -m pip uninstall -y flexiplex_filter
