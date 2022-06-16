# What does flexiplex do

Flexiplex is a light weight, flexible, error tolerant search and demultiplex tool. Given a set of reads as either .fastq or .fasta it will demultiplex and/or identify target sequences, reporting matching reads and a table of read-sequence assignments. It has been designed to demultuplex single cell long read RNA-Seq data, but can be used on any read data like an error tolerance "grep".

<pic>
  
<pic>


# Installing flexiplex
Download and untar the latex package:

If using linux or mac run you can run the code directly as
...
and 
... 
respectively
otherwise, compile the code with
make flexiplex.c++

# Examples of use

## Assigning single cell reads to 10x cellular barcodes
  
## Assigning single cell reads to 10x cellular barcodes (without knowing the barcodes)

## Demultiplexing bulk read data by barcode

## Assigning karyotype to cells




You can use the [editor on GitHub](https://github.com/DavidsonGroup/flexiplex/edit/gh-pages/index.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [Basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/DavidsonGroup/flexiplex/settings/pages). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.
