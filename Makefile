CXX?=g++

# set DIR to /usr/local/bin if not given through env
DIR?=bin
PYTHON?=python3

.PHONY: test clean install uninstall

all: flexiplex

flexiplex: flexiplex.c++ 
	${CXX} $^ -Ofast -pthread -std=c++11 -o $@ edlib-1.2.7/edlib.cpp -I edlib-1.2.7/ ${CFLAGS}

clean:
	rm flexiplex 

install:
	cp flexiplex ${DIR}
	${PYTHON} -m pip install scripts/

test: flexiplex
	@bash -c "if diff --color tests/output.fastq <(./flexiplex -k tests/barcodes.txt -b GGGG -x TTT tests/input.fastq 2>/dev/null); then echo 'Test success'; else echo 'Test fail'; fi"

uninstall:
	rm ${DIR}/flexiplex
	${PYTHON} -m pip list 2>/dev/null | grep -Fq flexiplex-filter && python3 -m pip uninstall -y flexiplex_filter
