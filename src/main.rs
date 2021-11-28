use bio::io::fastq;
use flate2::bufread;
use std::{fs, str, io::BufReader, collections::HashMap, process};

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

fn main() {
    let matches = App::new("fastkmers")
        .author("https://github.com/angelovangel")
        .about("get k-mer counts from a fastq file")

        .arg(Arg::with_name("kmer")
        .required(true)
        .help("k-mer size, maximal value is 31")
        .takes_value(true)
        .long("kmer_size")
        .short("k"))
        
        .arg(Arg::with_name("summary")
        .required(false)
        .long("summary")
        .short("s")
        .takes_value(false)
        .help("display fastq file summary information at the end of the output"))
        
        .arg(Arg::with_name("INPUT")
        .help("Path to a fastq file")
        .required(true)
        .index(1)).get_matches();

    let infile = matches.value_of("INPUT").unwrap().to_string();
    let reader = fastq::Reader::new(get_fastq_reader(&infile));
    //let mut record = fastq::Record::new();
    let mut kmer_counts: HashMap<String,i32> = HashMap::new();
    
    let k = matches.value_of("kmer").unwrap().trim().parse::<usize>().expect("k-mer length argument not valid!");
    
    if k >= 32 {
        println!("use k-mer size below the limit...");
        process::exit(0);
    }
    
    let mut reads: i64 = 0;
    let mut kmers: i64 = 0;

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
    for (key, value) in &kmer_counts {
        //let kmer_name = str::from_utf8(&key).unwrap();
        println!("{} \t {}", key, value)
    }

    if matches.is_present("summary") {
        let unique_kmers = kmer_counts.len();
        println!("---");
        println!("reads \t total_k-mers \t unique_k-mers");
        println!("{} \t {} \t {}", reads, kmers, unique_kmers);
        println!("---");
    }


    
}