# Contents

- [What does felxiplex do?](#what-does-felxiplex-do)
- [Installing flexiplex](#installing-flexiplex)
- [Usage](#usage)
- [Examples of use](#examples-of-use)
   - [Assigning single cell reads to 10x 3’ cellular barcodes (when barcodes are known)](#assigning-single-cell-reads-to-10x-3-cellular-barcodes-when-barcodes-are-known)
   - [Assigning single cell reads to 10x 3' cellular barcodes (when barcodes are unknown)](#assigning-single-cell-reads-to-10x-3-cellular-barcodes-when-barcodes-are-unknown)
   - [Demultiplexing other read data by barcode](#demultiplexing-other-read-data-by-barcode)
   - [Assigning genotype to cells](#assigning-genotype-to-cells)
   - [Simple search](#simple-search)
 - [Output](#output)
    - [New reads file](#new-reads-file)
    - [Table of barcodes found for each read](#table-of-barcodes-found-for-each-read)
    - [Table of barcode frequency](#table-of-barcode-frequency)
    - [Table of the number of barcode at each barcode frequency](#table-of-the-number-of-barcode-at-each-barcode-frequency)
 - [Support or Contact](#support-or-contact)


# What does flexiplex do?

Flexiplex is a light weight, flexible, error tolerant search and demultiplexing tool. Given a set of reads as either .fastq or .fasta it will demultiplex and/or identify target sequences, reporting matching reads and read-barcode assignment. It has been designed to demultiplex single cell long read RNA-Seq data, but can be used on any read data like an error tolerance "grep". Flexiplex is built with [edlib](https://github.com/Martinsos/edlib). 

Flexiplex first uses edlib to search for a left and right flanking sequence (primer and polyT) within each read (with barcode and UMI sequence left as a wildcard). For the best match with an edit distance of "f" or less it will trim to the barcode + UMI sequence +/- 5 bp either side, and search for the barcode against a known list. The best matching barcode with an edit distance of "e" or less will be reported. To identify second or more barcodes within a read, flexiplex repeats the search again with previously found primer to ployT sequence masked out. 

If the set of possible barcodes is unknown, flexiplex can be run in discovery mode (by leaving -k option off). In this mode, flexiplex will search for the primer and ployT sequence like usual, and take "b"bp after the primer sequence as a barcode. The frequency that barcodes are found in the data are reported for follow up analysis. For example, if 1000 cells were expected, the top 1000 most frequent barcodes can be used as the known list for a subsequent run of flexiplex.

The primer, ployT, list of barcodes and UMI length and maximum edit distances can all be set through user options (see [Usage](#usage)).

![Search sequence structure](/flexiplex/docs/assets/flexplex1.png)


# Installing flexiplex
Clone the [git repository](https://github.com/DavidsonGroup/flexiplex):
```
git clone https://github.com/DavidsonGroup/flexiplex.git
```

Change into the source directory and compile:
```
cd flexiplex ; make
```

You should now have a binary file called flexiplex which you can execute.
To see usage information, run 
```./flexiplex -h```


# Usage

```
FLEXIPLEX 0.9
usage: flexiplex [options] [reads_input]
  reads_input: a .fastq or .fasta file. Will read from stdin if empty.
  options:
     -k known_list   Either 1) a text file of expected barcodes in the first column,
                     one row per barcode, or 2) acomma separate string of barcodes.
                     Without this option, flexiplex will search and report possible barcodes.
                     The generated list can be used for known_list in subsequent runs.
     -r true/false   Replace read ID with barcodes+UMI, remove search strings
                     including flanking sequenence and split read if multiple
                     barcodes found (default: true).
     -s true/false   Sort reads into separate files by barcode (default: false)
     -p primer   Left flank sequence to search for (default: CTACACGACGCTCTTCCGATCT).
     -T polyT    Right flank sequence to search for (default: TTTTTTTTT).
     -n prefix   Prefix for output filenames.
     -b N   Barcode length (default: 16).
     -u N   UMI length (default: 12).
     -e N   Maximum edit distance to barcode (default 2).
     -f N   Maximum edit distance to primer+ployT (default 10).
     -h     Print this usage information.
```

# Examples of use

## Assigning single cell reads to 10x 3' cellular barcodes (when barcodes are known)
  
The default settings work for 10x 3' v3 chemistry. To demultiplex run:
```
flexiplex -k barcode_list.txt reads.fastq > new_reads.fastq
```

Or for older chemistry with 10bp UMIs use:
```
flexiplex -u 10 -k barcode_list.txt reads.fastq > new_reads.fastq
```

If dealing with large gzipped files you can pipe reads into flexiples. e.g.
```
gunzip -c read.fastq.gz | flexiplex -k barcode_list.txt | gzip > new_reads.fastq.gz
```
  
## Assigning single cell reads to 10x 3' cellular barcodes (when barcodes are unknown)

Flexiplex can be run in two passes: 1) to find the barcode sequences and 2) assign them to reads.
To find barcodes, set the flanking edit distance to 0 (a perfect match) as these are less likely to have errors in the barcodes.
```
flexiplex -f 0 reads.fastq
```

The table written to standard output can help to select the number of cells. For example, it can be captured and plotted in R or similar to make a knee plot.

Once the approximate number of cells is know, subset the list of barcodes by this number:
```
head -n <number of cells> flexiplex_barcodes_counts.txt > my_barcode_list.txt
```

Then use this list to assign barcodes to reads:
```
flexiplex -k my_barcode_list.txt reads.fastq > new_reads.fastq
```

## Demultiplexing other read data by barcode

To demultiplex with other flanking and barcodes sequences, set these 
```
flexiplex -p <left flank> -k "<barcode1>,<barcode2>" -T <right flank> -b <barcode length> -u 0 reads.fastq > new_reads.fastq
```

This assumes no UMI sequence is present. -e and -f which are the maximum barcode and flanking sequence edit distances may also need to be adjusted. As a guide we use -e 2 for 16bp barcodes and -f 12 for 32bp (left+right) flanking sequence.

If barcodes are expected at the start and end of reads, force flexiplex not to chop reads when mutiple barcodes are seen (-r false). Reads can also be separated into different files by barcodes (-s true).
```
flexiplex -r false -s true -p <left flank> -k "<barcode1>,<barcode2>" -T <right flank> -b <barcode length> -u 0
```

## Assigning genotype to cells

This is similar to [Demultiplexing other read data by barcode], but different alleles can be used in place of barcodes. e.g. to search for the KRAS variant c.34G>A run:

```
flexiplex -r false -p "GTATCGTCAAGGCACTCTTGCCTACGC" -k "CACTAGC,CACCAGC" -T "TCCAACTACCACAAGTTTATATTCAGT" -e 0 -f 15 -b 7 -u 0 reads.fasta > kras_var_reads.fasta
```

Where -k here lists the mutant and wild type variants (reverse complimented), with a few bp either side, and -p and -T are the adjacent sequence left and right of these respectively.

Multiple searches can be chained together. e.g. assign cellular barcodes then search for a specific mutation:
```
flexiplex -k barcode_list.txt reads.fasta | flexiplex -n barcode_mutation_mapping -r false -p "GTATCGTCAAGGCACTCTTGCCTACGC" -k "CACTAGC,CACCAGC" -T "TCCAACTACCACAAGTTTATATTCAGT" -e 0 -f 15 -b 7 -u 0 reads.fasta > kras_var_reads_with_barcodes.fasta
```

## Simple search

To perform a simple error tolerant grep-like search of a single sequence, split the sequence between the -p and -k (or -k and -T options):
```
flexiplex -r false -p "CACTCTTGCCTACGC" -k "CACTAGC" -f 3 -e 0 -b 7 reads.fasta
```

Matching reads will be printed to standard out. Edit distances (-e and -f ) can be adjusted as required.

# Output

## New reads file

Read with a matching barcode will be reported to standard output (or to individual files if the -s true option is provided).

If read chopping and ID replacement is used (-r true, default):
  - Read IDs will be replaced with the following format (similar to FLAMES): <barcode>_<UMI>#<original ID>_<+/-><N>of<M> where <+/-> indicates whether the barcode was found on the forward or reverse strand of the original read, M is the number of barcodes found in direction indicated by +/- and N is the 1st, 2nd etc. of those barcode.
  - Reads will be reverse complimented if the barcode was found in the reverse direction. For 10x 3' data this puts all reads in the reverse direction of the mRNA (3'->5')
  - If multiple barcodes are found in the same direction the read is split at the position of the second or subsequent primer, and multiple reads reported.
  - The primer+barcode+umi+polyT sequence is removed from the read.
If barcodes are found in both the forward and reverse directions on a read, the same read would be reported multiple time (once forward and once reverse). To overcome this duplication, data can be mapped as stranded

![Reads chopping](/flexiplex/docs/assets/flexplex2.png)
Schematic of default behaviour if multiple barcodes are identified in a read

  
 If read chopping and ID replacement is not used (-r false):
  - Any read with matching flank and barcode sequence will be reported with read ID appended with _<+/-> as above.
  - Read is reported only once even if multiple flank/barcode sequences found
  - If a flank/barcode is only found in the reverse direction, the read will be reverse complimented.

  
## Table of barcodes found for each read

This will be a file called flexiplex_reads_barcodes.txt or <prefix>_reads_barcodes.txt if the -n option was provided.
It is a tab delimited text file with a row for each barcode identified. e.g.
```
Read	CellBarcode	FlankEditDist	BarcodeEditDist	NextBestBarcodeEditDist
SRR12282458.3	CTCAGAAGTTTGCATG	4	0	100
SRR12282458.4	CGACTTCAGCTGTCTA	1	0	100
...
```
The NextBestBarcodeEditDist column is not currently filled.

## Table of barcode frequency

This will be a file called flexiplex_barcodes_counts.txt or <prefix>_reads_barcodes.txt if the -n option was provided. It will only be created if flexiplex is run in barcode discovery mode (no barcode list provided with -k).
 
It lists the number of reads that each barcode was found in, ordered by the more frequent. e.g. for a small dataset
```
GATCGATTCATCGCTC	15
AACCATGCACCTTGTC	12
ACGTCAAGTTTAAGCC	11
GTAACTGCATTGGGCC	11
AAGGCAGCATGGGAAC	10
ATGCGATTCACCTTAT	10
GATGCTAGTACAGCAG	10
GCAGCCATCATGCATG	10
GTTCATTCAATGTTGC	10
...
```
If the approximate number of cells is estimated, the top of this list can be used to generate a list of known barcodes e.g.
```
head -n <number of cell> flexiplex_barcodes_counts.txt > known_barcodes.txt
```
Then flexiplex run with -k known_barcodes.txt to more accurately assign the barcodes to reads (note that the second column in known_barcodes.txt gets ignored).
  
## Table of the number of barcode at each barcode frequency

This is printed to standard output when no barcodes are provided (ie. flexiplex is in barcode discovery mode and -k not provided).
e.g. for a small dataset:
```
Reads	Barcodes
15	1
14	0
13	0
12	1
11	2
10	5
9	6
8	13
7	22
6	17
5	29
4	62
3	86
2	124
1	1867
```
The first line can be interpreted as there was 1 barcode which was found in 15 reads etc. This data is provided to help work out how many cells were sequenced (e.g. by creating a knee plot).
  
# Support or Contact

To report issues, provide feedback or make a request please post a new [github issue](https://github.com/DavidsonGroup/flexiplex/issues). We are keen to hear if you are using flexiplex for a use case not described above or if you have suggestions for predefined settings/options (e.g. bulk sample demultiplexing).

   
