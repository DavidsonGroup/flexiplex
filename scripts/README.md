# scripts

This directory contains a Python tool, `flexiplex-filter`, which can be used to:
* Find an inflection point automatically
* Graph the knee plot of the barcode ranks/counts
* Filter using a whitelist of known barcodes

It is recommended that you run `flexiplex-filter` using the [`uv` package manager](https://docs.astral.sh/uv/getting-started/installation/):

```sh
uvx --from git+https://github.com/davidsongroup/flexiplex.git#subdirectory=scripts \
    flexiplex-filter --help
```

More information, including installation instructions, is detailed in the [usage file](usage.md).
