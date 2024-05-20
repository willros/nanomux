# nanomux
Demultiplex your nanopore reads

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
