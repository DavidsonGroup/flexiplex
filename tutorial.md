---
title: Long-read single-cell quickstart tutorial
layout: default
permalink: /tutorial.html
---

# Long-read single-cell quickstart tutorial

This tutorial shows a **minimal**, workflow for long-read single-cell RNA-seq for the tools flexiplex and nailpolish:

1. Extract 10x barcodes from long-read single-cell data using **flexiplex**  
2. Use **flexiplex-filter** to find the knee point and shortlist barcodes  
3. Demultiplex reads with **flexiplex**  
4. Clean reads and deduplicate UMIs with **nailpolish**
5. Align to a reference using **minimap2**
6. Quantify in **oarfish** (single-cell mode) to get a count matrix  
7. Load counts in **R/Seurat** and plot a UMAP

---

## Prerequisites (install)

This tutorial assumes that you have already installed the required software:
 * [**flexiplex** and **flexiplex-filter**](index.html)
 * [**nailpolish**](https://davidsongroup.github.io/nailpolish/quickstart.html)
 * [**minimap2**](https://github.com/lh3/minimap2)
 * [**oarfish**](https://github.com/COMBINE-lab/oarfish)
 * [**R** and **seurat**](https://satijalab.org/seurat/articles/install_v5.html)

You will also need to download the [**demo data**]()
```bash
wget .....scmixology2_250k.fastq.gz
```

Reference data (hg38 **chr1**):
```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.transcripts.fa.gz
gunzip gencode.v48.transcripts.fa.gz
```



## 1. Barcode discovery
First we will need to find out which single-cell barcodes are present in the dataset, which can be done by running **flexiplex** in discovery mode (with the barcode list paramter, -k, missing):

```bash
gunzip -c scmixology2_250k.fastq.gz | flexiplex -d 10x3v3 -f 0 > 1_flexiplex.out
```
The output should look something like:
```
FLEXIPLEX 1.02.1
Using predefined settings for 10x3v3.
Adding flank sequence to search for: CTACACGACGCTCTTCCGATCT
Setting barcode to search for: ????????????????
Setting UMI to search for: ????????????
Adding flank sequence to search for: TTTTTTTTT
Setting max flanking sequence edit distance to 8
Setting max barcode edit distance to 2
Setting max flanking sequence edit distance to 0
For usage information type: flexiplex -h
No filename given... getting reads from stdin...
Searching for barcodes...
0.01 million reads processed..
0.02 million reads processed..
0.03 million reads processed..
0.04 million reads processed..
0.05 million reads processed..
0.06 million reads processed..
0.07 million reads processed..
0.08 million reads processed..
0.09 million reads processed..
0.1 million reads processed..
0.2 million reads processed..
Number of reads processed: 250000
Number of reads where a barcode was found: 112261
Number of reads where more than one barcode was found: 3975
All done!
```
You'll see that the number of reads where a barcode was found is less than half, but this is okay because we have selected just the best quality read (-f 0 means no sequencing errors in the flank around the barcode).

The key output from this step is a file called "flexiplex_barcodes_counts.txt" which lists the number of reads for each barcode, and will be used for generating a knee plot in the next step.
```
CTCCGATCATGGCCAC	1090
CATCGCTCAAGTAGTA	1050
GGGACCTTCTTGATTC	1025
TCATTTGTCACGGTCG	1009
GGAGGATTCTTCTAAC	946
AATCACGGTCCTCCTA	943
ATAGACCCACCGTCTT	923
AGCGCCACAATCCAGT	907
GCCCAGACAACACAAA	881
ATGCCTCGTCAAGCCC	867
...
```

## 2. Barcode filtering

Now, we are ready to refine the barcode list to a short-list of high quality barcodes. **flexiplex-filter** is a python script packaged with flexiplex that will automatically identify the inflection point of the barcode knee plot, and (optionally) remove barcode not seen in known list - such as the inclusion lists provided by Cell Ranger. Check the flexiplex documentaion carefully to make sure you've correctly installed it, as it requires an extra step beyond the installiton of **flexiplex**.

To run:
```bash
flexiplex-filter flexiplex_barcodes_counts.txt > flexiplex_barcodes_final.txt
```
Optionally, if you have an inclusion list from Cell Ranger (but this is not required):
```bash
flexiplex-filter -w [inclusion_list] flexiplex_barcodes_counts.txt > flexiplex_barcodes_final.txt
```

It's always a good idea to manually check that flexiplex-filter has choosen the correct inflection point. This can be done by generating an interactive plot of the barcode frequencies:
```bash
flexiplex-filter -g flexiplex_barcodes_counts.txt
```
In this case, is appear to have choosen correct.
[image]

In our example, flexiplex-filter picks the 180 more frequent barcodes, which is close to the number of cells known in this dataset. The process above will also work for samples with many more cells, and has been tested on datasets with tens of thousands of cells.

# 2. Barcode demultiplexing

Now we have a short-list of barcodes, we can run flexiplex a second time to actually assign read to cellular barcodes:
```bash
gunzip -c scmixology2_250k.fastq.gz | flexiplex -d 10x3v3 -k flexiplex_barcodes_final.txt > scmixology2_250k.demultiplexed.fastq
```


