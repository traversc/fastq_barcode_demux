use rand::Rng;
use bio::pattern_matching::myers::Myers;
use std::fs::File;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::{Write, BufReader, BufRead};

const READ_LEN: usize = 80;
const READS_PER_BARCODE: usize = 3141;

const MAX_EDIT_DISTANCE: u32 = 2;
// Barcodes and reverse complement are at least 6 bp apart in test_barcodes.txt
// To-do: should be calculated dynamically
const BARCODE_DISTANCE: u32 = 6;

fn read_barcodes(path: &str) -> std::io::Result<Vec<String>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut barcodes = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if !trimmed.is_empty() {
            barcodes.push(trimmed.to_string());
        }
    }
    Ok(barcodes)
}

fn revcomp(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _   => c,
        })
        .collect()
}

// we need to make sure all possible alignments have the same edit_distance
// alignments may be less e.g. if a deletion cancels out an insertion
fn validate_mut_barcode(barcode: &str, candidate: &str, edit_distance: u32, position: &str) -> bool {
    if position == "start" {
        // reverse complement candidate and barcode and do find_all_end
        let barcode_rc = revcomp(barcode);
        let candidate_rc = revcomp(candidate);
        let pattern = Myers::<u64>::new(barcode_rc.as_bytes());
        let mut occ = pattern.find_all_end(candidate_rc.as_bytes(), edit_distance as u8);
        return ! occ.any(|(end, dist)| {
            // end is position not end of range, so we need to add +1
            (dist as u32) + ( ( candidate.len() - (end+1) ) as u32 ) < edit_distance
        });
    } else if position == "end" {
        let pattern = Myers::<u64>::new(barcode.as_bytes());
        let mut occ = pattern.find_all_end(candidate.as_bytes(), edit_distance as u8);
        return ! occ.any(|(end, dist)| {
            // end is position not end of range, so we need to add +1
            (dist as u32) + ( ( candidate.len() - (end+1) ) as u32 ) < edit_distance
        });
    } else { // position == middle
        let pattern = Myers::<u64>::new(barcode.as_bytes());
        let mut occ = pattern.find_all_end(candidate.as_bytes(), edit_distance as u8);
        return ! occ.any(|(_, dist)| {
            (dist as u32) < edit_distance
        });
    }
}

fn generate_mut_barcode(barcode: &str, edit_distance: u32, read_dir: &str, position: &str, rng: &mut impl Rng) -> String {
    let mut barcode = barcode.to_string();
    if read_dir.eq_ignore_ascii_case("revcomp") {
        barcode = revcomp(&barcode);
    }
    // Hack-y but should always work given edit_distance << barcode.len()
    // Generate a random barcode with the specified edit distance and validate it
    loop {
        let mut candidate = barcode.clone();
        for _ in 0..edit_distance {
            match rng.random_range(0..3) {
                // Insertion: insert "I" at a random position (including the end)
                0 => {
                    let pos = rng.random_range(0..=candidate.len());
                    candidate.insert(pos, 'I');
                },
                // Deletion: remove a nucleotide from a random position.
                1 => {
                    let pos = rng.random_range(0..candidate.len());
                    candidate.remove(pos);
                },
                // Mutation: change a nucleotide at a random position to "M".
                2 => {
                    let pos = rng.random_range(0..candidate.len());
                    candidate.replace_range(pos..pos+1, "M");
                },
                _ => unreachable!(),
            }
        }
        if validate_mut_barcode(&barcode, &candidate, edit_distance, position) {
            return candidate;
        }
    }
}

fn insert_barcode(base: &str, barcode: &str, pos: &str) -> String {
    let len = base.len();
    match pos.to_ascii_lowercase().as_str() {
        "start" => format!("{}{}", barcode, &base[barcode.len().min(len)..]),
        "end" => format!("{}{}", &base[..len.saturating_sub(barcode.len())], barcode),
        "middle" => {
            let mid = len / 2;
            let start = mid.saturating_sub(barcode.len() / 2);
            let end = start + barcode.len();
            format!("{}{}{}", &base[..start], barcode, &base[end..])
        },
        _ => base.to_string(),
    }
}

fn generate_fastq_read_pair(
    barcode: &str,
    edit_distance: u32,
    which_read: &str,  // "r1", "r2", or "both"
    position: &str,    // "start", "end", or "middle"
    read_dir: &str,    // "forward" or "revcomp"
    index: usize,
    rng: &mut impl Rng,
) -> (String, String) {
    let base_read = "X".repeat(READ_LEN);
    let err_barcode = generate_mut_barcode(barcode, edit_distance, read_dir, position, rng);
    let seq_r1 = if which_read.eq_ignore_ascii_case("r1") || which_read.eq_ignore_ascii_case("both") {
        insert_barcode(&base_read, &err_barcode, position)
    } else {
        base_read.clone()
    };
    let seq_r2 = if which_read.eq_ignore_ascii_case("r2") || which_read.eq_ignore_ascii_case("both") {
        insert_barcode(&base_read, &err_barcode, position)
    } else {
        base_read
    };
    let header_r1 = format!(
        "@{}_{}_{}_{}_{}_{:06} /1 {}",
        barcode, which_read, read_dir, position, edit_distance, index, err_barcode
    );
    let header_r2 = format!(
        "@{}_{}_{}_{}_{}_{:06} /2 {}",
        barcode, which_read, read_dir, position, edit_distance, index, err_barcode
    );
    let qual_r1 = "I".repeat(seq_r1.len());
    let qual_r2 = "I".repeat(seq_r2.len());
    let fastq_r1 = format!("{}\n{}\n+\n{}", header_r1, seq_r1, qual_r1);
    let fastq_r2 = format!("{}\n{}\n+\n{}", header_r2, seq_r2, qual_r2);
    (fastq_r1, fastq_r2)
}

fn main() -> std::io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: {} <barcode_file> <output_dir>", args[0]);
        std::process::exit(1);
    }
    let barcode_file = &args[1];
    let barcodes = read_barcodes(barcode_file)?;
    let output_dir = &args[2];
    std::fs::create_dir_all(output_dir)?;

    let mut rng = rand::rng();

    let r1_path = format!("{}/test_data_R1.fastq.gz", output_dir);
    let r2_path = format!("{}/test_data_R2.fastq.gz", output_dir);

    let r1_file = File::create(&r1_path)?;
    let r2_file = File::create(&r2_path)?;
    let mut writer_r1 = GzEncoder::new(r1_file, Compression::default());
    let mut writer_r2 = GzEncoder::new(r2_file, Compression::default());

    let mut index: usize = 0;
    for barcode in barcodes.iter() {
        for _ in 0..READS_PER_BARCODE {
            let edit_distance : u32 = rng.random_range(0..=MAX_EDIT_DISTANCE);
            let which_read = if rng.random_bool(0.5) { "r1" } else { "r2" };
            let position = match rng.random_range(0..3) {
                0 => "start",
                1 => "end",
                _ => "middle",
            };
            let read_dir = if rng.random_bool(0.5) { "forward" } else { "revcomp" };
            let (fastq_r1, fastq_r2) = generate_fastq_read_pair(
                barcode,
                edit_distance,
                which_read,
                position,
                read_dir,
                index,
                &mut rng
            );
            writeln!(writer_r1, "{}", fastq_r1)?;
            writeln!(writer_r2, "{}", fastq_r2)?;
            index += 1;
        }
    }
    // Generate multimap read pairs (r1 and r2 have different barcodes).
    // set edit_distances to zero so we can always identify these cases
    for _ in 0..READS_PER_BARCODE {
        let edit_distance_r1 = 0;
        let edit_distance_r2 = 0;
        // Choose two distinct barcodes.
        let b1 = &barcodes[rng.random_range(0..barcodes.len())];
        let mut b2 = &barcodes[rng.random_range(0..barcodes.len())];
        while b2 == b1 {
            b2 = &barcodes[rng.random_range(0..barcodes.len())];
        }
        let position = match rng.random_range(0..3) {
            0 => "start",
            1 => "end",
            _ => "middle",
        };
        let read_dir_r1 = if rng.random_bool(0.5) { "forward" } else { "revcomp" };
        let read_dir_r2 = if rng.random_bool(0.5) { "forward" } else { "revcomp" };

        // Call generate_fastq_read_pair separately for b1 and b2.
        // We use "r1" for the first call and "r2" for the second call.
        let (mut fastq_r1, _) = generate_fastq_read_pair(b1, edit_distance_r1, "r1", position, read_dir_r1, index, &mut rng);
        let (_, mut fastq_r2) = generate_fastq_read_pair(b2, edit_distance_r2, "r2", position, read_dir_r2, index, &mut rng);

        // Modify the headers so they begin with "@MULTIMAP" instead of the barcode.
        let mut lines: Vec<String> = fastq_r1.lines().map(|s| s.to_string()).collect();
        lines[0] = format!("@MULTIMAP_{:06} r1", index);        
        fastq_r1 = lines.join("\n");

        let mut lines: Vec<String> = fastq_r2.lines().map(|s| s.to_string()).collect();
        lines[0] = format!("@MULTIMAP_{:06} r2", index);
        fastq_r2 = lines.join("\n");

        writeln!(writer_r1, "{}", fastq_r1)?;
        writeln!(writer_r2, "{}", fastq_r2)?;
        index += 1;
    }

    // NOMAP
    // create reads with edit distance 3, 4 or 5
    // Since input barcodes are >= 6 edit distance apart, these reads should not map to any barcode during testing.
    for _ in 0..READS_PER_BARCODE {
        // choose a random barcode
        let barcode = &barcodes[rng.random_range(0..barcodes.len())];
        let edit_distance : u32 = rng.random_range((MAX_EDIT_DISTANCE+1)..BARCODE_DISTANCE);
        let which_read = if rng.random_bool(0.5) { "r1" } else { "r2" };
        let position = match rng.random_range(0..3) {
            0 => "start",
            1 => "end",
            _ => "middle",
        };
        let read_dir = if rng.random_bool(0.5) { "forward" } else { "revcomp" };
        let (fastq_r1, fastq_r2) = generate_fastq_read_pair(
            barcode,
            edit_distance,
            which_read,
            position,
            read_dir,
            index,
            &mut rng
        );
        writeln!(writer_r1, "{}", fastq_r1)?;
        writeln!(writer_r2, "{}", fastq_r2)?;
        index += 1;
    }

    writer_r1.finish()?;
    writer_r2.finish()?;
    Ok(())
}
