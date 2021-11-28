
# fastkmers

A simple program for getting k-mer counts from a fastq file, written in Rust.

## Description

This command line program takes a fastq file as input (can be `.gz` also) and outputs the counts of [k-mers](https://en.wikipedia.org/wiki/K-mer) of a specified length. It is implemented using hash table and a simple algortihm but is still reasonably fast.

## Install

I do not provide precompiled binaries here, but it is simple to compile and run:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
git clone https://github.com/angelovangel/faster.git

cd fastkmers
cargo build --release

# the binary is now under ./target/release/, run it like this:
./target/release/fastkmers /path/to/fastq/file.fastq.gz

```

## Usage


```bash
# to get 4-mer counts
fastkmers -k 4 file.fastq.gz

```

The k-mer counts are printed to `stdout` as a tab-separated table.
## Reference

`fastkmers` uses the excellent Rust-Bio library:

[Köster, J. (2016). Rust-Bio: a fast and safe bioinformatics library. Bioinformatics, 32(3), 444-446.](https://academic.oup.com/bioinformatics/article/32/3/444/1743419)
