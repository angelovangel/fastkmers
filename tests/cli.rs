use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::process::Command; // Run programs


#[test]

fn find_content_in_file() -> Result<(), Box<dyn std::error::Error>> {

    let mut cmd = Command::cargo_bin("fastkmers")?;
    cmd.arg("-k 3").arg("tests/test.fastq");

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("CGT \t 6"));

    Ok(())
}
