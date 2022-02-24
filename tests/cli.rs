use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::process::Command; // Run programs


#[test]
fn test_kmer_fastq() -> Result<(), Box<dyn std::error::Error>> {

    let mut cmd = Command::cargo_bin("fastkmers")?;
    cmd.arg("-k 3").arg("tests/test.fastq");

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("ATC\t10"));

    Ok(())
}

#[test]
fn test_kmer_fasta() -> Result<(), Box<dyn std::error::Error>> {

    let mut cmd = Command::cargo_bin("fastkmers")?;
    cmd.arg("-k 3").arg("tests/test.fasta");

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("ATC\t10"));

    Ok(())
}


#[test]
fn test_freq() -> Result<(), Box<dyn std::error::Error>> {

    let mut cmd = Command::cargo_bin("fastkmers")?;
    
    cmd.arg("-k 3").arg("-f").arg("tests/test.fastq");
    cmd.assert().success().stdout(predicate::str::contains("10\t5"));

    Ok(())
}

#[test]
fn test_regex() -> Result<(), Box<dyn std::error::Error>> {

    let mut cmd = Command::cargo_bin("fastkmers")?;
    
    cmd.arg("-k 5").arg("-q").arg("A[T|G]A$").arg("tests/test.fastq");
    cmd.assert().success().stdout(predicate::str::contains("GAAGA\t3"));

    Ok(())

}

#[test]
fn test_cycle() -> Result<(), Box<dyn std::error::Error>> {
    
    let mut cmd = Command::cargo_bin("fastkmers")?;
    cmd.arg("-k 126").arg("-c").arg("tests/test.fastq");
    cmd.assert().success().stdout(predicate::str::contains("1\t0.5\t0.25\t0\t0.25\t0"));

    Ok(())
}