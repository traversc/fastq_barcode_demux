use std::process::Command;
use tempfile::{TempDir};
use itertools::iproduct;
use std::path::{Path, PathBuf};
use std::env;
use std::io::{BufWriter, Write, BufReader, BufRead, Read};
use std::fs::File;
use flate2::read::MultiGzDecoder;

// temporary dir, but you can prevent it from deletion by setting
// the environment variable KEEP_TEMP_DIR
pub struct MaybeTempDir {
    inner: Option<TempDir>,
    path: PathBuf,
}
impl MaybeTempDir {
    pub fn new() -> Self {
        let temp_dir = tempfile::tempdir().expect("Failed to create temporary directory");
        let path = temp_dir.path().to_path_buf();
        let inner = if env::var("KEEP_TEMPDIR").is_ok() {
            // Persist the directory (prevent auto deletion)
            // into_path consumes the tempDir, preventing deletion when temp_dir goes out of scope
            let _ = temp_dir.into_path();
            None
        } else {
            Some(temp_dir)
        };
        MaybeTempDir { inner, path }
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn into_path(mut self) -> PathBuf {
        if let Some(temp_dir) = self.inner.take() {
            temp_dir.into_path()
        } else {
            self.path
        }
    }
}

#[test]
fn integration_test_generate_and_demux() {
    // Create a temporary directory for test output.
    let temp_dir = MaybeTempDir::new();
    let temp_path = temp_dir.path().to_str().unwrap();

    // Run generate_fastq
    // test_barcodes are 19 bp long and at least an edit distance of 6 to each other AND their reverse complement
    println!("Running generate_fastq");
    let gen_output = Command::new("cargo")
        .args(&["run", "--bin", "generate_fastq", "--", "tests/test_barcodes.txt", temp_path])
        .output()
        .expect("Failed to run generate_fastq");
        println!("generate_fastq output: {}", String::from_utf8_lossy(&gen_output.stdout));
    assert!(
        gen_output.status.success(),
        "generate_fastq failed: {}",
        String::from_utf8_lossy(&gen_output.stderr)
    );

    // Check that the two FASTQ files were created.
    let r1_path = temp_dir.path().join("test_data_R1.fastq.gz");
    let r2_path = temp_dir.path().join("test_data_R2.fastq.gz");
    assert!(
        r1_path.exists(),
        "R1 FASTQ file not found at: {}",
        r1_path.display()
    );
    assert!(
        r2_path.exists(),
        "R2 FASTQ file not found at: {}",
        r2_path.display()
    );

    let max_edits = vec!["0", "1", "2"];
    // Match Direction/Match Position s.t. 
    // Match Direction in {Forward, Reverse Complement, Either} and 
    // Match Position in {Start of read, End of read, Anywhere}
    let match_config = vec![
        "F/S", "F/E", "F/A",
        "RC/S", "RC/E", "RC/A",
        "E/S", "E/E", "E/A",
    ];

    // read in line by line
    // test/test_barcodes.txt is a list of barcode sequence
    let barcodes = std::fs::read_to_string("tests/test_barcodes.txt").unwrap();
    let barcodes = barcodes.lines(); // iterator
    let barcode_names: Vec<_> = barcodes.clone().enumerate().map(|(i, _)| format!("bc{}", i)).collect();
    let configs = iproduct!(max_edits, match_config);
    for (max_edit, mconf) in configs {
        println!("Testing max_edit: {}, mconf: {}", max_edit, mconf);
        // create a temporary barcode config CSV file

        // for each barcode, column 1 is barcode name bc0,bc1, etc.
        // column 2 is the barcode sequence
        // column 3 is the max_edit
        // column 4 is match config for R1
        // column 5 is match config for R2
        let barcode_config = temp_dir.path().join("barcode_config.csv");
        let mut writer = BufWriter::new(File::create(&barcode_config).unwrap());
        for (i, barcode) in barcodes.clone().enumerate() {
            writeln!(writer, "bc{},{},{},{},{}", i, barcode, max_edit, mconf, mconf).unwrap();
        }
        writer.flush().unwrap();

        let r1_fastq = temp_dir.path().join("test_data_R1.fastq.gz");
        let r2_fastq = temp_dir.path().join("test_data_R2.fastq.gz");
        // Run demux_barcodes single-threaded
        // demux_barcodes <R1 FASTQ> <R2 FASTQ> <barcodes> --output <output directory> --nthreads 1
        println!("Running single-threaded demux_barcodes");
        let st_output_dir = MaybeTempDir::new();
        let status = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--bin",
            "demux_barcodes",
            "--",
            r1_fastq.to_str().unwrap(),
            r2_fastq.to_str().unwrap(),
            barcode_config.to_str().unwrap(),
            "--output-path",
            st_output_dir.path().to_str().unwrap(),
            "--nthreads",
            "1",
        ])
        .status()
        .unwrap();
        assert!(status.success());

        println!("Running multi-threaded demux_barcodes");
        let mt_output_dir = MaybeTempDir::new();
        let status = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--bin",
            "demux_barcodes",
            "--",
            r1_fastq.to_str().unwrap(),
            r2_fastq.to_str().unwrap(),
            barcode_config.to_str().unwrap(),
            "--output-path",
            mt_output_dir.path().to_str().unwrap(),
            "--nthreads",
            "4",
        ])
        .status()
        .unwrap();
        assert!(status.success());

        // File outputs are gzip fastq
        // Check that
        // 1) all the file names are the same
        // 2) All the file contents are the same
        let st_files: Vec<_> = std::fs::read_dir(st_output_dir.path()).unwrap().collect();
        let mt_files: Vec<_> = std::fs::read_dir(mt_output_dir.path()).unwrap().collect();
        assert_eq!(st_files.len(), mt_files.len(), "Output directories have different number of files for {} {}", max_edit, mconf);
        for (st_entry, mt_entry) in st_files.into_iter().zip(mt_files.into_iter()) {
            let st_path = st_entry.unwrap().path();
            let mt_path = mt_entry.unwrap().path();
            assert_eq!(st_path.file_name(), mt_path.file_name(), "File names do not match");
            
            let st_file = File::open(&st_path).unwrap();
            let mt_file = File::open(&mt_path).unwrap();
            
            let mut st_decoder = MultiGzDecoder::new(st_file);
            let mut mt_decoder = MultiGzDecoder::new(mt_file);
            
            let mut st_contents = String::new();
            let mut mt_contents = String::new();
            st_decoder.read_to_string(&mut st_contents).unwrap();
            mt_decoder.read_to_string(&mut mt_contents).unwrap();
            
            assert_eq!(st_contents, mt_contents, "File contents do not match for file {:?}", st_path.file_name());
        }
        
        // The output files are <sample_name>_R1,2__<barcode_name>.fastq.gz
        // Read in R1 and R2 files at the same time
        // Check that the headers (first element before space) are the same
        let st_dir = st_output_dir.path();
        let mut r1_files: Vec<_> = std::fs::read_dir(st_dir)
            .unwrap()
            .filter_map(|e| {
                let p = e.unwrap().path();
                let name = p.file_name()?.to_string_lossy();
                if name.contains("test_data_R1") { Some(p) } else { None }
            })
            .collect();
        let mut r2_files: Vec<_> = std::fs::read_dir(st_dir)
            .unwrap()
            .filter_map(|e| {
                let p = e.unwrap().path();
                let name = p.file_name()?.to_string_lossy();
                if name.contains("test_data_R2") { Some(p) } else { None }
            })
            .collect();
        
        r1_files.sort();
        r2_files.sort();
        assert_eq!(r1_files.len(), r2_files.len());
        
        for (r1_path, r2_path) in r1_files.iter().zip(r2_files.iter()) {
            let file_name = r1_path.file_name().unwrap().to_string_lossy(); // lossy here means drop any non-utf8 info
            let parts: Vec<&str> = file_name.split("__").collect();
            let barcode_name = parts[1].strip_suffix(".fastq.gz").unwrap();
            // get the corresponding barcode (sequence)
            let barcode = if barcode_name == "MULTIMAP" || barcode_name == "NOMAP" {
                barcode_name
            } else {
                barcodes.clone().nth(barcode_names.iter().position(|x| x == barcode_name).unwrap()).unwrap()
            };
            println!("Checking barcode: {} {}", barcode_name, barcode);
            let r1_headers: Vec<String> = BufReader::new(MultiGzDecoder::new(File::open(r1_path).unwrap()))
                .lines()
                .enumerate()
                .filter_map(|(i, line)| if i % 4 == 0 { Some(line.unwrap()) } else { None })
                .collect();
        
            let r2_headers: Vec<String> = BufReader::new(MultiGzDecoder::new(File::open(r2_path).unwrap()))
                .lines()
                .enumerate()
                .filter_map(|(i, line)| if i % 4 == 0 { Some(line.unwrap()) } else { None })
                .collect();
            assert_eq!(r1_headers.len(), r2_headers.len(), "Header count mismatch between R1 and R2");

            for (hdr1, hdr2) in r1_headers.iter().zip(r2_headers.iter()) {
                let hdr1_parts: Vec<&str> = hdr1.split_whitespace().collect();
                let hdr2_parts: Vec<&str> = hdr2.split_whitespace().collect();
                // Check that the header (before the space) is identical for both files.
                assert_eq!(hdr1_parts[0], hdr2_parts[0], "Header info mismatch between R1 and R2");
                // Check that the read indicator (after the space) is different.
                assert_ne!(hdr1_parts[1], hdr2_parts[1], "Read indicator should differ between R1 and R2");
            
                if barcode == "MULTIMAP" {
                    // Every read pair in this file should have a header that starts with "MULTIMAP"
                    assert!(hdr1_parts[0].starts_with("@MULTIMAP"), "Expected header to start with MULTIMAP, got {}", hdr1_parts[0]);
                } else if barcode == "NOMAP" {
                    if hdr1_parts[0].starts_with("@MULTIMAP") { // skip for now, assessment is complicated
                        continue;
                    }
                    // Extract variables from the header (expected format: @<barcode>_<which_read>_<position>_<direction>_<sim_edit>_<index>)
                    let tokens: Vec<&str> = hdr1_parts[0]
                        .trim_start_matches('@')
                        .split('_')
                        .collect();
                    assert_eq!(tokens.len(), 6, "Expected 6 tokens in header, got {}", tokens.len());
                    // let extracted_barcode = tokens[0];
                    // let which_read = tokens[1];
                    let direction = tokens[2];
                    let position = tokens[3];
                    let sim_edit = tokens[4];

                    // expect a read pair in this file if any of the following conditions are met:
                    // sim_edit > max_edit
                    // If match_direction == F and direction == revcomp
                    // If match_direction == RC and direction == forward
                    // If match position == S and position == end or middle
                    // If match position == E and position == start or middle
                    let parts: Vec<&str> = mconf.split("/").collect();
                    let match_direction = parts[0];
                    let match_position = parts[1];

                    // these are pairs that did NOT map to barcodes
                    let expect_read_pair = 
                    sim_edit.parse::<u32>().unwrap() > max_edit.parse::<u32>().unwrap() ||
                        (match_direction == "F" && direction == "revcomp") ||
                        (match_direction == "RC" && direction == "forward") ||
                        (match_position == "E" && (position == "start" || position == "middle")) ||
                        (match_position == "S" && (position == "end" || position == "middle"));
                    if !expect_read_pair {
                        panic!("Unexpected read pair in NOMAP file: match_direction {} match_position {}, header {} {} {} {}",
                            match_direction, match_position,
                            hdr1, direction, position, sim_edit);
                    }
                } else {
                    if hdr1_parts[0].starts_with("@MULTIMAP") {
                        continue;
                    }
                    // Extract variables from the header (expected format: @<barcode>_<which_read>_<position>_<direction>_<sim_edit>_<index>)
                    let tokens: Vec<&str> = hdr1_parts[0]
                        .trim_start_matches('@')
                        .split('_')
                        .collect();
                    assert_eq!(tokens.len(), 6, "Expected 6 tokens in header, got {}", tokens.len());
                    let extracted_barcode = tokens[0];
                    // let which_read = tokens[1];
                    let direction = tokens[2];
                    let position = tokens[3];
                    let sim_edit = tokens[4];
                    assert_eq!(extracted_barcode, barcode, "Extracted barcode does not match expected barcode");

                    // The logic here is a little complicated because edit distance calcs depend on alignment mode
                    // There are two cases:
                    // Case 1: match_position == A
                    // (match_direction == E) OR (match_direction == F and direction == forward) OR (match_direction == RC and direction == revcomp)
                    // max_edit may be less than sim_edit
                    // i.e. the number of allowed can be less than how many we generated, since semi-global is not forced to match the end or start
                    // Case 2: match_position (search parameter) is the same as position (simulation parameter)
                    // the number of simulated edits is strictly less than how many we allow: sim_edit <= max_edit
                    // (match_direction == E) OR (match_direction == F and direction == forward) OR (match_direction == RC and direction == revcomp)
                    let parts: Vec<&str> = mconf.split("/").collect();
                    let match_direction = parts[0];
                    let match_position = parts[1];

                    let expect_pair_case1 = 
                        match_position == "A" &&
                        ((match_direction == "E") ||
                         (match_direction == "F" && direction == "forward") ||
                         (match_direction == "RC" && direction == "revcomp"));
                    let expect_pair_case2 = 
                        ( (match_position == "S" && position == "start") || (match_position == "E" && position == "end") ) &&
                        sim_edit.parse::<u32>().unwrap() <= max_edit.parse::<u32>().unwrap() &&
                        ((match_direction == "E") ||
                         (match_direction == "F" && direction == "forward") ||
                         (match_direction == "RC" && direction == "revcomp"));
                    if !(expect_pair_case1 || expect_pair_case2) {
                        panic!("Unexpected read pair in file: match_direction {} match_position {}, header {} {} {} {}",
                            match_direction, match_position,
                            hdr1, direction, position, sim_edit);
                    }
                }
            }
        }
    }
}
