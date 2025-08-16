---
title: Long-read single-cell quickstart
layout: default
permalink: /tutorial.html
---

# Long-read single-cell quickstart

This walkthrough shows a **minimal**, end-to-end pipeline for long-read single-cell RNA-seq:

1. Extract 10x barcodes from long-read sc data  
2. Use `flexiplex-filter` to find the knee point and shortlist barcodes  
3. Demultiplex reads with `flexiplex`  
4. Clean & deduplicate with **nailpolish** (barcode + UMI consensus)  
5. Quantify in **oarfish** (single-cell mode) to get a count matrix  
6. Load counts in **R/Seurat** and plot a UMAP

> Tested on Linux/macOS. Commands are intentionally short and copy-pasteable.

---

## Prerequisites (quick install)

We’ll install the core tools with `conda` and use the latest **nailpolish** binary.

```bash
# (Recommended) new env
mamba create -n lrsc -c conda-forge -c bioconda \
  flexiplex oarfish minimap2 samtools coreutils pigz -y
mamba activate lrsc

# flexiplex-filter (installed with flexiplex via conda) OR run via uv:
# uvx --from git+https://github.com/davidsongroup/flexiplex.git#subdirectory=scripts flexiplex-filter --help

# nailpolish (download latest nightly or release for your platform)
curl -L "https://github.com/DavidsonGroup/nailpolish/releases/download/nightly_develop/nailpolish" -o nailpolish
chmod +x nailpolish
