# What does flexiplex do?

Flexiplex is a light weight, flexible, error tolerant search and demultiplex tool. Given a set of reads as either .fastq or .fasta it will demultiplex and/or identify target sequences, reporting matching reads and read-sequence assignment. It has been designed to demultuplex single cell long read RNA-Seq data, but can be used on any read data like an error tolerance "grep". Flexiplex is built on [edlib](https://github.com/Martinsos/edlib). 


# Installing flexiplex
Clone the [git repository](https://github.com/DavidsonGroup/flexiplex):
```git clone https://github.com/DavidsonGroup/flexiplex.git```

change into the source directory and compile:
```cd flexiplex ; make```

You should now have a binary file called flexiplex which you can execute

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

# Output

### new reads file

### table of barcodes found for each read

### table of barcode frequency

### table of the number of barcode at each barcode frequency


# Examples of use

## Assigning single cell reads to 10x cellular barcodes
  
The default setting work for 10x 3' v3 chemistry. To demultiplex run:
```flexiplex -k barcode_list.txt reads.fastq > new_reads.fastq```

For 10bp UMIs use:
```flexiplex -u 10 reads.fastq > new_reads.fastq```

If dealing with large gzipped files you can pipe reads into flexiples. e.g.
```gunzip -c read.fastq.gz | flexiplex -k barcode_list.txt | gzip > new_reads.fastq.gz```
  
## Assigning single cell reads to 10x cellular barcodes (without knowing the barcodes)

Flexiplex can be run in two passes: 1) to find the barcode sequences and 2) assign them to reads.
To find barcodes, set the flanking edit distance to 0 (a perfect match) as these are less likely to have errors in the barcodes.
```flexiplex -f 0 reads.fastq```

The table written to standard output can help to select the number of cells. For example, it can be captured and plotted in R or similar to make a knee plot.

Once the approximate number of cells is know, subset the list of barcodes by this number:
```head -n <number of cells> flexiplex_barcodes_counts.txt > my_barcode_list.txt```

Then use this list to assign barcodes to reads:
```flexiplex -k my_barcode_list.txt reads.fastq > new_reads.fastq```

## Demultiplexing other read data by barcode

To demultiplex with other flanking and barcodes sequences, set these 
```flexiplex -p <left flank> -k "<barcode1>,<barcode2>" -T <right flank> -b <barcode length> -u 0 reads.fastq > new_reads.fastq```

This assumes no UMI sequence is present. -e and -f which are the maximum barcode and flanking sequence edit distances may also need to be adjusted. As a guide we use -e 2 for 16bp barcodes and -f 12 for 32bp (left+right) flanking sequence.

If barcodes are expected at the start and end of reads, force flexiplex not to chop reads when mutiple barcodes are seen (-r false). Reads can also be separated into different files by barcodes (-s true).
``` flexiplex -r false -s true -p <left flank> -k "<barcode1>,<barcode2>" -T <right flank> -b <barcode length> -u 0```

## Assigning karyotype to cells

This is similar to [Demultiplexing other read data by barcode], but different alleles can be used in place of barcodes. e.g. to search for the KRAS variant .. run:

```flexiplex -r false -p "GTATCGTCAAGGCACTCTTGCCTACGC" -k "CACTAGC,CACCAGC" -T "TCCAACTACCACAAGTTTATATTCAGT" -e 0 -f 15 -b 7 -u 0 reads.fasta > kras_var_reads.fasta```

Multiple searches can be strung together. e.g. assign cellular barcodes then search for mutation:
```flexiplex -k barcode_list.txt reads.fasta | flexiplex -n barcode_mutation_mapping -r false -p "GTATCGTCAAGGCACTCTTGCCTACGC" -k "CACTAGC,CACCAGC" -T "TCCAACTACCACAAGTTTATATTCAGT" -e 0 -f 15 -b 7 -u 0 reads.fasta > kras_var_reads_with_barcodes.fasta```

## Simple search

To perform a simple error tolerant grep-like search of a single sequence, split the sequence between the -p and -k (or -k and -T options):
```flexiplex -r false -p "CACTCTTGCCTACGC" -k "CACTAGC" -f 3 -e 0 -b 7 reads.fasta```

Matching reads will be printed to standard out. Edit distances (-e and -f ) can be adjusted as required.

### Support or Contact

To report issues or provide feedback please add a [github issue](https://github.com/DavidsonGroup/flexiplex/issues)
