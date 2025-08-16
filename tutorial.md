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
```wget .....
gunzip 
```

Reference data (hg38 **chr1**):
```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.transcripts.fa.gz
gunzip gencode.v48.transcripts.fa.gz
``



## 1. Barcode discovery
First we will need to find out which single-cell barcodes are present in the dataset. 







