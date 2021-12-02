use bio::io::{fastq, fasta};
use flate2::bufread;
use std::{fs, str, io::BufReader, collections::HashMap, process};
use serde::{Deserialize, Serialize};

extern crate clap;
use clap::{App, Arg};


// own functions
// fastq reader, file as arg, decide based on extension
fn get_fastq_reader(path: &String) -> Box<dyn (::std::io::Read)> {
    if path.ends_with(".gz") {
        let f = fs::File::open(path).unwrap();
        Box::new(bufread::MultiGzDecoder::new(BufReader::new(f)))
    } else {
        let f = fs::File::open(path).unwrap();
        Box::new(BufReader::new(f))
    }
}

// define struct used for serializing the HashMap to json
#[derive(Serialize, Deserialize)]
struct Kmer {
    bases: String,
    count: i32,
}

fn main() {
    let matches = App::new("fastkmers")
        .author("https://github.com/angelovangel")
        .about("get k-mer counts and multiplicity frequency from a fastq file")

        .arg(Arg::with_name("kmer")
        .required(true)
        .help("k-mer size, maximal value is 21 (required)")
        .takes_value(true)
        .long("kmer_size")
        .short("k"))
        
        .arg(Arg::with_name("summary")
        .required(false)
        .long("summary")
        .short("s")
        .takes_value(false)
        .help("display fastq file summary information at the end of the output (optional)")).
        
        arg(Arg::with_name("fasta")
        .required(false)
        .long("fasta")
        .short("a")
        .takes_value(false)
        .help("input is fasta file (default is fastq)"))
        
        .arg(Arg::with_name("json")
        .conflicts_with("summary")
        .required(false)
        .long("json")
        .short("j")
        .takes_value(false)
        .help("output data in json format (optional)"))
        
        .arg(Arg::with_name("freq")
        .long("freq")
        .short("f")
        .required(false)
        .takes_value(false)
        .help("Output a frequency table (multiplicity versus occurence), for building k-mer spectra, see https://en.wikipedia.org/wiki/K-mer (optional)"))
        
        .arg(Arg::with_name("INPUT")
        .help("Path to a fastq/fasta file")
        .required(true)
        .index(1)).get_matches();

    let infile = matches.value_of("INPUT").unwrap().to_string();
    // fastq or fasta
    //let mut record = fastq::Record::new();
    let mut kmer_counts: HashMap<String,i32> = HashMap::new();
    
    let k = matches.value_of("kmer").unwrap().trim().parse::<usize>().expect("k-mer length argument not valid!");
    
    if k >= 22 {
        println!("use k-mer size below the limit...");
        process::exit(0);
    }
    
    let mut reads: i64 = 0;
    let mut kmers: i64 = 0;
    
    if matches.is_present("fasta") {
        let reader = fasta::Reader::new(get_fastq_reader(&infile));

        for result in reader.records(){
            reads += 1;

            let record = result.expect("Error");
            let seq = record.seq();
            let seq_str = str::from_utf8(seq).unwrap().to_string();
            //println!("{:?}", seq);
            for c in 0..seq_str.len() - k + 1 {
                let subseq = &seq_str[c..c + k];
            
                kmers += subseq.len() as i64;
        
                *kmer_counts.entry(subseq.to_string() ).or_insert(0) += 1;
            }
        }
    } else {
        let reader = fastq::Reader::new(get_fastq_reader(&infile));

        for result in reader.records(){
            reads += 1;

            let record = result.expect("Error");
            let seq = record.seq();
            let seq_str = str::from_utf8(seq).unwrap().to_string();
            //println!("{:?}", seq);
            for c in 0..seq_str.len() - k + 1 {
                let subseq = &seq_str[c..c + k];
            
                kmers += subseq.len() as i64;
        
                *kmer_counts.entry(subseq.to_string() ).or_insert(0) += 1;
            }
        }
    }

    if matches.is_present("json") {
        let j = serde_json::to_string(&kmer_counts).unwrap();
        println!("{}", j);
        process::exit(0);
    }

    if matches.is_present("freq") {
        
        let values: Vec<i32> = kmer_counts.values().cloned().collect(); // get counts from main hash table
        let mut freq_hash: HashMap<i32, i32> = HashMap::new();
        println!("occ\tcount");

        for i in values {
            
            *freq_hash.entry(i).or_insert(0) += 1;
        }
        
        for (key, value) in &freq_hash {
            println!("{}\t{}", key, value)
        }
        process::exit(0);
    }

    println!("kmer\tcount");
    for (key, value) in &kmer_counts {
        //let kmer_name = str::from_utf8(&key).unwrap();
        println!("{}\t{}", key, value)
    }

    if matches.is_present("summary") {
        let unique_kmers = kmer_counts.len();
        println!("---");
        println!("reads\ttotal_k-mers\tunique_k-mers");
        println!("{}\t{}\t{}", reads, kmers, unique_kmers);
        println!("---");
    }
    
}