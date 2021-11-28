
# fastkmers

A simple program for getting k-mer counts from a fastq file, written in Rust.

## Description

This command line program takes a fastq file as input (can be `.gz` also) and outputs the counts of [k-mers](https://en.wikipedia.org/wiki/K-mer) of a specified length. It is implemented using hash table and a simple algortihm but is still reasonably fast. The maximum supported k-mer size is 21.

## Install

I do not provide precompiled binaries here, but it is simple to compile and run:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
git clone https://github.com/angelovangel/fastkmers.git

cd fastkmers
cargo build --release

```

## Usage


```bash

# the binary is now under ./target/release/, run it like this:
./target/release/fastkmers -k 4 /path/to/fastq/file.fastq.gz

# to get 4-mer counts and a summary
fastkmers -k 4 -s file.fastq.gz

# output json
fastkmers -k 4 -j file.fastq.gz

```

The k-mer counts are printed to `stdout` as a tab-separated table or as `json`.
