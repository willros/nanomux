# nanomux
Demultiplex your dual barcodeded 16s Nanopore reads. Uses either `greedy` or `fuzzy` search to identify the barcodes.

## Install:
```console
$ git clone https://github.com/willros/nanomux.git
$ cd nanomux
$ pip install -e .
```

## Quick Start:
```console
$ nanomux -h
```
## Barcode csv
The `.csv` file with barcodes must have the **following** shape:
```csv
name,fwd_barcode,rvs_barcode
sample1,ATCGTAGCGA,GTCAGCGTTG
sample2,AGTTCGATTG,GATGCGATTT
```
