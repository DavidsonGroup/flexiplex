<div align="center">
	<h1>Flexiplex</h1>
	<p>
		<b>A versatile demultiplexer and search tool for omics data</b>
	</p>
</div>

Flexiplex is a fast, multithreaded, and user-configurable demultiplexer. Given a set of reads as either FASTQ or FASTA, it will demultiplex and/or identify a sequence of interest, reporting matching reads and read-barcode assignment. Flexiplex works in two modes: (i) when one or more sequences of interest are known, such as barcodes, and (ii) discovery mode—when only the sequence which flanks the region of interest is known.

For a detailed description of how Flexiplex works and compares against other tools, please see our [paper](https://academic.oup.com/bioinformatics/article/40/3/btae102/7611801).

## Usage

Below are a few things that Flexiplex can do, but we recommend you [**check out the documentation for a more detailed guide**](https://davidsongroup.github.io/flexiplex/).

**Assign reads - short reads or single-cell long reads - to cellular barcodes**  
Presets for various chemistries, including 10x Chromium and Visium, are [available in the docs](https://davidsongroup.github.io/flexiplex/#assigning-long-reads-to-10x-barcodes-when-barcodes-are-known).
```sh
flexiplex -d 10x3v3 -k barcode_list.txt reads.fastq > new_reads.fastq
```

**Barcode discovery directly from long reads**  
Flexiplex can search for barcodes, filter against a whitelist, and remove empty droplets using an inflection point method.
```sh
flexiplex -f 0 reads.fastq
flexiplex-filter --whitelist 3M-february-2018.txt --outfile filtered.txt flexiplex_barcodes_counts.txt
```

**Error-tolerant search**  
Flexiplex can perform an error-tolerant grep-like search for a single sequence.
```sh
flexiplex -x "CACTCTTGCCTACGCCACTAGC" -d grep -f 3 reads.fasta
```

**And more!**  
For documentation visit [https://davidsongroup.github.io/flexiplex/](https://davidsongroup.github.io/flexiplex/).

For usage information type `flexiplex -h` and `flexiplex-filter -h`.


## Installation

Precompiled binaries for Flexiplex are located in the [Releases](https://github.com/DavidsonGroup/flexiplex/releases) section. You can also install Flexiplex using Anaconda: `conda install -c bioconda -c conda-forge flexiplex`

flexiplex-filter can be installed locally using Make, but we recommend using the [`uv` package manager](https://docs.astral.sh/uv/getting-started/installation/), which will automatically manage dependencies in a temporary virtualenv. Just invoke flexiplex-filter using

```sh
uvx --from git+https://github.com/davidsongroup/flexiplex.git#subdirectory=scripts \
    flexiplex-filter --help
```

### Compiling from source

To compile flexiplex, ensure that gcc is installed, then run:
`make`

To install flexiplex-filter and binaries into /usr/local/bin, ensure that `python3.9`, and `python-pip` are installed, then run:
`make install`

To uninstall, run:
`make uninstall`
