


EDLIB=../edlib-1.2.7
#HD5=/opt/homebrew/Cellar/hdf5/1.12.1_1
#HD5_CC=h5c++
CC=g++
#DEP=ssw.h ssw.c ssw_cpp.cpp ssw_cpp.h ONT_barcode_extractor.c++ extract_passed_barcodes.c++
#SSW=./Complete-Striped-Smith-Waterman-Library/src
GPROF=/opt/homebrew/Cellar/gperftools/2.9.1_1

#all: ONT_barcode_extractor.c++ extract_passed_barcodes.c++
#	${CC} -O0 extract_passed_barcodes.c++ -o extract_passed_barcodes ${HD5}
#	${CC} -Ofast ONT_barcode_extractor.c++ -o ONT_barcode_extractor ${EDILIB} ${HD5}

all: flexiplex #extract_passed_barcodes

#extract_passed_barcodes: extract_passed_barcodes.c++
#	${HD5_CC} -Ofast -o $@ $^ -I${HD5}/include 

flexiplex: flexiplex.o #${SSW}/ssw_cpp.o ${SSW}/ssw.o
	${CC} $^ -o $@ -L${EDLIB}/meson-build -ledlib -L${GPROF}/lib/ -lprofiler 

flexiplex.o: flexiplex.c++
	${CC} -Ofast -c flexiplex.c++ -I${EDLIB}/edlib/include/ -I${GPROF}/include/ #-I${SSW}

clean:
	rm flexiplex.o flexiplex #extract_passed_barcodes
