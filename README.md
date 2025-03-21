# This repo is no longer maintained
Go to my other project: [nanomux_c](https://github.com/willros/nanomux_c), which is implemented in c and is much faster compared to this.  

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

## Usage
```console
usage: nanomux [-h] -f FASTX -o OUTPUT -b BARCODES -m {greedy,fuzzy} [-mm {1,2,3,4}] [-bc_start BARCODE_START] [-bc_end BARCODE_END] [-len_min_start READ_LEN_MIN_START]
               [-len_min_after READ_LEN_MIN_AFTER] [-len_max READ_LEN_MAX] [-min_reads MINIMUM_READS] [-par | --parquet | --no-parquet] [-t | --trim | --no-trim]
               [-nin | --number_in_name | --no-number_in_name]

Demultiplex your Nanopore reads!

options:
  -h, --help            show this help message and exit
  -f FASTX, --fastx FASTX
                        Fastx file
  -o OUTPUT, --output OUTPUT
                        Output folder
  -b BARCODES, --barcodes BARCODES
                        .csv files with barcodes
  -m {greedy,fuzzy}, --mode {greedy,fuzzy}
                        Mode to find barcodes. Either 'greedy' or 'fuzzy'
  -mm {1,2,3,4}, --mismatch {1,2,3,4}
                        If fuzzy, how many mismatches are allowed? [DEFAULT]: 1
  -bc_start BARCODE_START, --barcode_start BARCODE_START
                        Where do the barcode start in the read?
  -bc_end BARCODE_END, --barcode_end BARCODE_END
                        Where do the barcode end in the read?
  -len_min_start READ_LEN_MIN_START, --read_len_min_start READ_LEN_MIN_START
                        Minimum length of read before searching for barcodes
  -len_min_after READ_LEN_MIN_AFTER, --read_len_min_after READ_LEN_MIN_AFTER
                        Minimum length of read after finding barcodes
  -len_max READ_LEN_MAX, --read_len_max READ_LEN_MAX
                        Maximum length of read
  -min_reads MINIMUM_READS, --minimum_reads MINIMUM_READS
                        Minimum amount of reads to keep the barcode
  -par, --parquet, --no-parquet
                        Save all demuxed files as parquet file? (default: True)
  -t, --trim, --no-trim
                        Trim the reads? (default: False)
  -nin, --number_in_name, --no-number_in_name
                        Number of reads in the fastq name? (default: False)
```

## Example
```console
nanomux -f tests/test.fastq -o 240524_TEST_greedy -b tests/barcodes_with_flanking_sequence.csv -m greedy -min_reads 30
```


## Changelog

### 240524
* Trim option
* Capped mismatch to 1
* Chunked search for barcodes to limit memory consumtion


## 240603
* Minimun length before finding barcodes and after trimming found barcodes
* Reports barcodes with no reads assigned to it in the summary csv-file
* Option to keep number of reads in fastq name or not
