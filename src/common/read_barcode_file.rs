use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use crate::common::barcode_matcher::{ReadDir, MatchPos};

/// Parses a configuration string of the form "BarcodeDirection/MatchPosition".
///
/// Allowed BarcodeDirection values:
///   - Forward match only: "F", "Fwd", "Forward"
///   - Reverse Complement match only: "RC", "RevComp", "Reverse Compliment"
///   - Allow match to either direction: "E", "Either"
///
/// Allowed MatchPosition values:
///   - Match only at the start: "S", "Start"
///   - Match only at the end: "E", "End"
///   - Match anywhere: "A", "Any", "Anywhere"
///
/// If the string is blank or "NA", returns (ReadDir::NoMatchAllowed, MatchPos::NoMatchAllowed).
fn parse_config(config: &str) -> Result<(ReadDir, MatchPos), Box<dyn Error>> {
    let config = config.trim();
    if config.is_empty() || config.eq_ignore_ascii_case("NA") {
        return Ok((ReadDir::NoMatchAllowed, MatchPos::NoMatchAllowed));
    }
    let parts: Vec<&str> = config.split('/').collect();
    if parts.len() != 2 {
        panic!("Invalid configuration format: '{}'. Expected format 'BarcodeDirection/MatchPosition'", config);
    }
    let dir_str = parts[0].trim().to_ascii_lowercase();
    let pos_str = parts[1].trim().to_ascii_lowercase();

    let dir = match dir_str.as_str() {
        "f" | "fwd" | "forward" => ReadDir::Forward,
        "rc" | "revcomp" | "reverse compliment" => ReadDir::ReverseComplement,
        "e" | "either" => ReadDir::Either,
        _ => panic!("Invalid barcode direction: '{}'", parts[0]),
    };

    let pos = match pos_str.as_str() {
        "s" | "start" => MatchPos::Start,
        "e" | "end" => MatchPos::End,
        "a" | "any" | "anywhere" => MatchPos::Anywhere,
        _ => panic!("Invalid match position: '{}'", parts[1]),
    };

    Ok((dir, pos))
}

/// Reads a barcode file and returns 6 vectors, one per column:
/// - Names (String)
/// - Barcodes (String)
/// - Maximum barcode edit distance (u8)
/// - Read1 configuration (String, see readme)
/// - Read2 configuration (String, see readme)
///
/// Each line in the barcode file must have either 2 or 5 tab-separated columns.
/// For lines with 2 columns, max edit distance is 1 and Read 1 and 2 config are E/A (see readme)
/// Output: (names, sequences, max_distance, read1 barcode dir, read1 match position, read2 barcode dir, read2 match position)
pub fn read_barcode_file(
    path: &str,
) -> Result<
    (
        Vec<String>,
        Vec<String>,
        Vec<u8>,
        Vec<ReadDir>, // read 1
        Vec<MatchPos>, // read 1
        Vec<ReadDir>, // read 2
        Vec<MatchPos>, // read 2
    ),
    Box<dyn Error>,
> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut names = Vec::new();
    let mut sequences = Vec::new();
    let mut max_distances = Vec::new();
    let mut read1_dirs = Vec::new();
    let mut read1_matchpos = Vec::new();
    let mut read2_dirs = Vec::new();
    let mut read2_matchpos = Vec::new();

    for (line_no, line_result) in reader.lines().enumerate() {
        let line = line_result?;
        // Skip empty lines and comments.
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let delimiter = if line.contains('\t') {
            '\t'
        } else if line.contains(',') {
            ','
        } else {
            panic!("Line does not contain a tab or comma delimiter");
        };
        let parts: Vec<&str> = line.split(delimiter).collect();
        match parts.len() {
            2 => {
                // Two-column line: [name, sequence]
                names.push(parts[0].to_string());
                sequences.push(parts[1].to_string());
                // Default configuration is "E/A" for both reads.
                let (default_dir, default_pos) = (ReadDir::Either, MatchPos::Anywhere);
                max_distances.push(1);
                read1_dirs.push(default_dir);
                read1_matchpos.push(default_pos);
                read2_dirs.push(default_dir);
                read2_matchpos.push(default_pos);
            }
            5 => {
                // Five-column line: [name, sequence, read1 config, read2 config, max edit distance]
                names.push(parts[0].to_string());
                sequences.push(parts[1].to_string());
                let max_distance: u8 = parts[2].trim().parse().unwrap_or_else(|e| {
                    panic!("Error parsing max edit distance on line {}: {}", line_no + 1, e)
                });          
                let (r1_dir, r1_pos) = parse_config(parts[3])?;
                let (r2_dir, r2_pos) = parse_config(parts[4])?;
                max_distances.push(max_distance);
                read1_dirs.push(r1_dir);
                read1_matchpos.push(r1_pos);
                read2_dirs.push(r2_dir);
                read2_matchpos.push(r2_pos);
            }
            _ => {
                panic!("Invalid number of columns on line {}: expected 2 or 5, found {}", line_no + 1, parts.len());
            }
        }
    }

    Ok((names, sequences, max_distances, read1_dirs, read1_matchpos, read2_dirs, read2_matchpos))
}