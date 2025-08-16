---
title: Long-read single-cell quickstart tutorial
layout: default
permalink: /tutorial.html
---

# Long-read single-cell tutorial

This tutorial shows a simple workflow for pre-processing long-read single-cell RNA-seq using our tools flexiplex and nailpolish:

1. Extract 10x barcodes from long-read single-cell data using **flexiplex**  
2. Use **flexiplex-filter** to find the knee point and shortlist barcodes  
3. Demultiplex reads with **flexiplex**  
4. Clean reads and deduplicate UMIs with **nailpolish**
5. Align to a reference using **minimap2**
6. Quantify in **oarfish** (single-cell mode) to get a count matrix  
7. Load counts in **R/Seurat** and plot a UMAP

---

## 0. Prerequisites

This tutorial assumes that you have already installed the required software:
 * [**flexiplex** and **flexiplex-filter**](index.html)
 * [**nailpolish**](https://davidsongroup.github.io/nailpolish/quickstart.html)
 * [**minimap2**](https://github.com/lh3/minimap2)
 * [**oarfish**](https://github.com/COMBINE-lab/oarfish)
 * [**R** and **seurat**](https://satijalab.org/seurat/articles/install_v5.html)

You will also need to download the demo dataset:
```bash
wget .....scmixology2_250k.fastq.gz
```
This dataset comes from scmixology2 from []() and rebase called by [](). It contains sequencing single cell data from a mixture of 5
lung cancer cells lines sequenced using 10x 3' version 3 and Oxford Nanopore Technologies. 
We have subset down to just 250 thousand reads from the 250 most variable genes, to make the tutrial run in a reasonable time.

You will also need to download the human reference transcriptome:
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
FLEXIPLEX 1.02.5
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
Number of reads where more than one barcode was found: 1422
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
![flexiplex-filter graph](/flexiplex/docs/assets/tutorial.flexiplex-filter.png)

In our example, flexiplex-filter picks the 180 more frequent barcodes, which is close to the number of cells known in this dataset. The process above will also work for samples with many more cells, and has been tested on datasets with tens of thousands of cells.

## 3. Barcode demultiplexing

Now we have a short-list of barcodes, we can run flexiplex a second time to actually assign read to cellular barcodes:
```bash
gunzip -c scmixology2_250k.fastq.gz | flexiplex -d 10x3v3 -k flexiplex_barcodes_final.txt > scmixology2_250k.demultiplexed.fastq
```

The output of this look like:
```
FLEXIPLEX 1.02.5
Using predefined settings for 10x3v3.
Adding flank sequence to search for: CTACACGACGCTCTTCCGATCT
Setting barcode to search for: ????????????????
Setting UMI to search for: ????????????
Adding flank sequence to search for: TTTTTTTTT
Setting max flanking sequence edit distance to 8
Setting max barcode edit distance to 2
Setting known barcodes from flexiplex_barcodes_final.txt
Number of known barcodes: 180
For usage information type: flexiplex -h
No filename given... getting reads from stdin...
File flexiplex_reads_barcodes.txt already exists, overwriting.
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
Number of reads where a barcode was found: 213458
Number of reads where more than one barcode was found: 5923
All done!
```
And you see that we are able to find barcodes in 84% if reads. The read ID now contains the barcode and UMI information needed for downstream analysis:
scmixology2_250k.demultiplexed.fastq
```bash
@GCCCGAACAATACCCA_GCATAAATTGTA#82f7a317-c081-4d86-bbb3-60926dce0db9_+1of1	CB:Z:GCCCGAACAATACCCA	UB:Z:GCATAAATTGTA
TTTTTTTTTTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCATATTGGCCAGGCTGGTCTCGAACTCCTGACCTCATCATCCACCTGCCTTAGCCTCCCAAAATGCTGGGATTACAGGCATGAGCCACCGCACCCGGCCAGGAATTGCTTCTTATGGATGAACAAAGAAAATAGTTTCCTGAGATGGTATCTACTCTGTGAAGATGCACAACATTGTTGAAATGACAACAAAGGATCTGAATATTACATAAACTTGTTGATAAAAACAGTGGCAGGGTTTGAGATCATTGATTCAGTTTTGAAAGAAGTTCTACTGTGGGTAAAAATGCCATCAAACAGCATCACATGCTGCAGAGAAATCTTTCATGAAAGGAAGAATCAATCGATATGTATCCTCCATTTTGTCTTATTTTAAAAATTTGCCACAGCCACTCCAACCTTCAGTAACTACTACCTTGATCTGTCAGCAGCCATCAACATCAAGGCAACAGACCCTCCCACCAGCAAAAAAATTACAACTTGCTGAAGGCTCAGATGATCATTACCATTTTTTAACAATATTTTAAGATAGGTGTGTACATTGTTTTTAGACATAATCTATTGCATACTTAATAGACTTCTACAGTAATAGTGCAAAACATAACCTTTATGTACACTAGGGAACCAAAAATTTACACGTGACTCCTTTTATTTTGATGGTCTGGAACTGAGTCCACAGTATCTCCAAGGTATGCCTATATTGACTGATGAATAAACAAAATGTGGTATAGCCATACAATGGAATATTGTTCAGTCACAAAAAGAATAACATTCTGATACATGCTACAACATAGATGAACCTTGAAAATATTATGCTAAGTGAAAGAAGTCAGACACAAAAGACAAGTATTATATGATTCCATTTATATAAAATGTCCAGGATAGTCAAATCAACGAGAGACCAAAAAGTAGATTAATGGTTGACAGAGGCTGAAGGTAGGGAGGAAATGGGGAGTGACTTTAATGGGTACAGGGTTTTAGATAAAAGAAAATGTTCTGGAATTAGTGGTGATGCACAACCTTGTACAACCTTTTTGCAATATAGTGCAATAAGACAATATTGTACAATCTTGTTTATATACTAAAGCCATTAAATTACACACTTTAAAGGGTGCATTTTATGACATGGGAAGTGTACATTAAAAATAAATAACCCCATGTACTCTGCGTTGATACCACTGCTTGGCCATTAGGCCTGTAAACAATACGTAACTTG
...
```

## 4. Deduplication and polishing

As single-cell data often contains multiple reads for each unique barcode and UMI combination, we need to deduplicate reads to a single barcode/UMI combination. This is typically achieved after alingment and during quantification (ie. only counting a barcode and UMI combination once). However, **nailpolish** allows this deduplication to be performed earlier, and uses the duplicate reads to generate conconsensus sequence - hence improving the read accuracy whilst also deduplicating.

Nailpolish works in two steps. First the fastq file is indexed.
```bash
nailpolish index scmixology2_250k.demultiplexed.fastq
```
At this point you can use the output to examine the duplication rate within the data. 
```
nailpolish summary scmixology2_250k.demultiplexed.fastq
```

Loading scmixology2_250k.demultiplexed.summary.html in a browser shows that the duplicate rate is low. That's because this tutorial uses a small subset the data, but a range of 5-40% would be typical for long-read single-cell data.
[image]

Now we can perform the concensus calling:
```bash
nailpolish consensus scmixology2_250k.demultiplexed.fastq > scmixology2_250k.demultiplexed.deduplicated.fastq
```
The resulting data will be a mixture of original (singlet) and consensus (duplicated) reads, which each have a unique barcode and UMI. The barcode and UMI are encoded in the read ID which will be transfers into the mapped bam files by minimap2.

## 5. Read Mapping

We are now ready to align the reads. As we'll be using oarfish for quantification, mapping is done again the reference transcriptome. I this instance gencode (downloaded in step 0).
```bash
minimap2 -x map-ont -y gencode.v48.transcripts.fa scmixology2_250k.demultiplexed.deduplicated.fastq > scmixology2_250k.demultiplexed.deduplicated.sam
```
The "-y" flag here is important as it tell minimap2 to add the barcode and UMI tags into the bam (these will be used by oarfish):
```bash
head scmixology2_250k.demultiplexed.deduplicated.sam
processed_0_1	1256	221	484	-	ENST00000206423.8|ENSG00000091986.16|OTTHUMG00000159265.4|OTTHUMT00000354219.2|CCDC80-201|CCDC80|12301|protein_coding|	12301	4687	4953	64	266	4	tp:A:P	cm:i:6	s1:i:63	s2:i:55	dv:f:0.1056	rl:i:109	MI:Z:GCCCGAACAATACCCA_GCATAAATTGTA	nI:i:0	CB:Z:GCCCGAACAATACCCA	UB:Z:GCATAAATTGTA	nT:Z:simplex
processed_0_1	1256	221	530	+	ENST00000357401.3|ENSG00000197934.9|OTTHUMG00000078661.2|OTTHUMT00000171614.1|CYYR1-AS1-201|CYYR1-AS1|3412|lncRNA|	3412	1765	2077	56	313	0	tp:A:S	cm:i:4	s1:i:55	dv:f:0.1394	rl:i:109	MI:Z:GCCCGAACAATACCCA_GCATAAATTGTA	nI:i:0CB:Z:GCCCGAACAATACCCA	UB:Z:GCATAAATTGTA	nT:Z:simplex
processed_0_1	1256	144	484	+	ENST00000651069.1|ENSG00000180769.10|OTTHUMG00000141281.8|OTTHUMT00000502559.2|WDFY3-AS2-206|WDFY3-AS2|2179|lncRNA|	2179	1403	1748	51	346	0	tp:A:S	cm:i:5	s1:i:49	dv:f:0.1280	rl:i:109	MI:Z:GCCCGAACAATACCCA_GCATAAATTGTA	nI:i:0CB:Z:GCCCGAACAATACCCA	UB:Z:GCATAAATTGTA	nT:Z:simplex
```

Before proceeding we also need to sort the



## 6. Transcript Quantification

Next we quantify the read 

## 7. Count analysis

At this point your count matrix can be analysed is whatever way makes the most sense for you and your experiment. Here we show a simple example using R and seurat to generate a UMAP of the cells. 



