CXX?=g++

# set DIR to /usr/local/bin if not given through env
DIR?=bin
PYTHON?=python3

all: flexiplex

flexiplex: flexiplex.c++ 
	${CXX} $^ -Ofast -pthread -std=c++11 -o $@ edlib-1.2.7/edlib.cpp -I edlib-1.2.7/ ${CFLAGS}

clean:
	rm flexiplex 

install:
	cp flexiplex ${DIR}
	${PYTHON} -m pip install scripts/

uninstall:
	rm ${DIR}/flexiplex
	${PYTHON} -m pip list 2>/dev/null | grep -Fq flexiplex-filter && python3 -m pip uninstall -y flexiplex_filter
