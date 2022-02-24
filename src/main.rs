use std::{collections::HashMap, process};
use regex::Regex;
use kseq::parse_path;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

extern crate clap;
use clap::{App, Arg};

fn main() {
    let matches = App::new("fastkmers")
        .version("0.1.3")
        .author("https://github.com/angelovangel")
        .about("get k-mer counts and multiplicity frequency from a fastx file")

        .arg(Arg::with_name("cycle")
        .required(false)
        .takes_value(false)
        .long("cycle")
        .short("c")
        //.conflicts_with("kmer")
        .help("Print bases per cycle (optional). Number of cycles must be given with the -k argument Use only with Illumina data."))

        .arg(Arg::with_name("kmer")
        //.default_value("0")
        .required(true)
        .help("k-mer size")
        .takes_value(true)
        .long("kmer_size")
        .short("k"))
        
        .arg(Arg::with_name("summary")
        .required(false)
        .long("summary")
        .short("s")
        .takes_value(false)
        .help("display fastq file summary information at the end of the output (optional)"))

        .arg(Arg::with_name("valid")
        .long("valid")
        .short("v")
        .required(false)
        .takes_value(false)
        .help("Remove k-mers with ambigous bases, e.g. only the ones containing ATGC are retained"))
        
        .arg(Arg::with_name("json")
        .conflicts_with("summary")
        .required(false)
        .long("json")
        .short("j")
        .takes_value(false)
        .help("output data in json format (optional)"))
        
        .arg(Arg::with_name("query")
        .long("query")
        .short("q")
        .conflicts_with("summary")
        .conflicts_with("freq")
        .required(false)
        .takes_value(true)
        .help("Output counts for a k-mer sequence passed as query on the command line. Works with regex too!"))
        
        .arg(Arg::with_name("freq")
        .long("freq")
        .short("f")
        .required(false)
        .takes_value(false)
        .help("Output a histogram of k-mer occurence, for building k-mer spectra, see https://en.wikipedia.org/wiki/K-mer (optional)"))
        
        .arg(Arg::with_name("INPUT")
        .help("Path to a fastq/fasta file")
        .required(true)
        .index(1))
        
        .get_matches();

    let infile = matches.value_of("INPUT").unwrap().to_string();
    let mut records = parse_path(infile).unwrap();

    let mut kmer_counts: HashMap<String,i32> = HashMap::new();
    
    
    let k = matches.value_of("kmer")
        .unwrap()
        .trim()
        .parse::<usize>()
        .expect("k-mer length argument not valid!");
    
    let mut reads: i64 = 0;
    let mut kmers: i64 = 0;
    
    
    while let Some(record) = records.iter_record().unwrap() {
        reads += 1;
        let seq_str = record.seq();

        for c in 0..seq_str.len() - k + 1 {
            let subseq = &seq_str[c..c + k];

            if matches.is_present("valid") {
                if subseq.chars().all(|x| matches!(x, 'A'|'T'|'G'|'C'|'a'|'t'|'g'|'c') ) {

                    kmers += subseq.len() as i64;
                    *kmer_counts.entry(subseq.to_string() ).or_insert(0) += 1;
                }
            } else {

            kmers += subseq.len() as i64;
            *kmer_counts.entry(subseq.to_string() ).or_insert(0) += 1;
            
            }
        }
    }
// get per cycle content by first running the main function with k = cycles (to make sure the reads are aligned) 
// then calculate base contents at each position
    if matches.is_present("cycle") {
        
        //let mut b: Vec<char> = Vec::with_capacity(k);
        println!("cycle\ta\tt\tg\tc\tn");
        
       // for i in 0..k {
       //     let b: Vec<char> = kmer_counts
       //         .keys()
       //         .map(|x| x.chars().nth(i).unwrap())
       //         .collect();
        (0..k).into_par_iter().for_each(|i| {
            let b: Vec<char> = kmer_counts
                .keys()
                .map(|x| x.chars().nth(i).unwrap())
                .collect();

                let a: f32 = b.iter().filter(|&x| matches!(x, 'A'|'a')).count() as f32 / b.len() as f32;
                let t: f32 = b.iter().filter(|&x| matches!(x, 'T'|'t')).count() as f32 / b.len() as f32;
                let g: f32 = b.iter().filter(|&x| matches!(x, 'G'|'g')).count() as f32 / b.len() as f32;
                let c: f32 = b.iter().filter(|&x| matches!(x, 'C'|'c')).count() as f32 / b.len() as f32;
                let n: f32 = b.iter().filter(|&x| matches!(x, 'N'|'n')).count() as f32 / b.len() as f32;
            
                println!("{}\t{}\t{}\t{}\t{}\t{}", i+1, a, t, g, c, n);
        });

           // let a: f32 = b.iter().filter(|&x| matches!(x, 'A'|'a')).count() as f32 / b.len() as f32;
           // let t: f32 = b.iter().filter(|&x| matches!(x, 'T'|'t')).count() as f32 / b.len() as f32;
           // let g: f32 = b.iter().filter(|&x| matches!(x, 'G'|'g')).count() as f32 / b.len() as f32;
           // let c: f32 = b.iter().filter(|&x| matches!(x, 'C'|'c')).count() as f32 / b.len() as f32;
           // let n: f32 = b.iter().filter(|&x| matches!(x, 'N'|'n')).count() as f32 / b.len() as f32;
           // 
           // println!("{}\t{}\t{}\t{}\t{}\t{}", i+1, a, t, g, c, n);
            //println!("{:?}", b.len());
        
        
       
        //modules::cycle_vec(records);
        process::exit(0);
    }
    
    if matches.is_present("json") {
        let j = serde_json::to_string(&kmer_counts).unwrap();
        println!("{}", j);
        process::exit(0);
    }

    if matches.is_present("query") {
        let string = matches.value_of("query").unwrap().trim();
        let re = Regex::new(string).expect("Failed to construct regex from string!");
        println!("kmer\tcount");
        
        for (key, value) in kmer_counts {
            if re.is_match(&key) {
                println!("{}\t{}", key, value)
            }
        }
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
        println!("reads\ttotal_kmers\tunique_kmers");
        println!("{}\t{}\t{}", reads, kmers, unique_kmers);
        println!("---");
    }
    
}