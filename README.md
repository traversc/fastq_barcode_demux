# Fastq Barcode Demux

[![CI](https://github.com/traversc/fastq_barcode_demux/actions/workflows/ci.yml/badge.svg)](https://github.com/traversc/fastq_barcode_demux/actions/workflows/ci.yml)


- A flexible tool for demultiplexing FASTQ files containing barcodes
- Splits FASTQ files into separate output files based on barcode sequence
- Supports both `zstd` and `gzip` compression formats
- Written in Rust using `rust-bio` and `rayon`

## Usage

```
demux_barcodes \
    sample_R1.fastq.gz \
    sample_R2.fastq.gz \
    barcode_config.txt
```

### Additional Options
- **--output-path:** Specify the directory where the output files will be written
- **--nthreads:** Set the number of threads to use
- **--compress-algo:** Choose the output compression algorithm (gz/gzip or zstd)
- **--compress-level:** Define the compression level (e.g., 0 means default; valid ranges: 0-9 for gzip, 0-22 for zstd).
- **--verbose:** Enable verbose output (0 or 1).

## Installation

Rust must be installed:
- via official shell script (https://www.rust-lang.org/tools/install)
- or via conda `conda install -c conda-forge rust`

Compile executable:
```
cargo build --release
```

Executable `demux_barcodes` is located in `targets` folder and can be moved anywhere. 

## The barcode config file

- The barcode config file controls how barcodes are searched in each read pair
- There is a simple and advanced config format
- The configuration files can be tab or comma delimited

**Simple configuration**

Two columns, no column header: `barcode_name` and `barcode_sequence`. This will search for barcode sequences anywhere within the read pair, with a maximum edit distance of 1. If found, reads will be demultiplexed into files named based on `barcode_name`.

Don't use any weird symbols for `barcode_name` as it will be part of the output file path.

**Advanced configuration**

Five columns, no column header:

`barcode_name`, `barcode_sequence`, `max_edit_distance`, `Read1_config`, `Read2_config` (See next section)

## Read Config: Barcode Direction / Match Position Configuration

For each read in the pair, this parameter is specified as a coupled value with two components separated by a slash. The first component is the Barcode Direction and the second is the Match Position.

### Barcode Direction Allowed Values:

    Forward barcode sequence only: `F`, `Fwd`, `Forward`
    Reverse Complement barcode sequence only: `RC`, `RevComp`, `Reverse Complement`
    Allow matches to either direction: `E`, `Either`

### Match Position Allowed Values:

    Match only at the start of the read: `S`, `Start`
    Match only at the end of the read: `E`, `End`
    Match anywhere in the read: `A`, `Any`, `Anywhere`

### Usage Example
Specify as BarcodeDirection/MatchPosition. E.g., `F/A` for searching for the forward barcode sequence with a match allowed anywhere. Leave the value blank or use `NA` if not allowing any match for that read in the pair.
