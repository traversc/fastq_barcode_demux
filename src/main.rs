mod common;

use crate::common::barcode_matcher::{BarcodeMatcher, UniqueMatch};
use crate::common::read_barcode_file::read_barcode_file;
use crate::common::concurrent_barcode_writer::ConcurrentBarcodeWriter;
use crate::common::input_decompressor_traits::{open_input_decompressor};
use crate::common::output_compressor_traits::{OutputCompressor, Finish, open_output_compressor};

use std::io::{self, BufRead, Write};

use clap::Parser;

#[derive(clap::Parser, Debug)]
#[clap(author, version, about)]
struct Args {
    reads1_path: String,
    reads2_path: String,
    barcodes_path: String,
    #[clap(long, default_value = "output", help = "Output directory")]
    output_path: String,
    #[clap(long, default_value_t = 1, help = "Number of threads to use")]
    nthreads: usize,
    #[clap(
        long,
        default_value = "gz",
        value_parser = clap::builder::PossibleValuesParser::new(["gz", "gzip", "zstd"]),
        help = "Output compression algorithm"
    )]
    compress_algo: String,
    #[clap(
        long,
        default_value_t = 0,
        help = "Compression level: 0-9 for gzip, 0-22 for zstd, 0 means default (gzip: 6, zstd: 3)"
    )]
    compress_level: u32,
    #[clap(
        long,
        default_value_t = 0,
        help = "Verbosity level (0 or 1)"
    )]
    verbose: usize,
}

// hard coded multithreading parameters
// Rayon's approach of join/waiting has overhead, so we chunk reads into batches and process 
// multiple batches at a time in order to avoid the overhead of spawning tasks.
// release
#[cfg(not(debug_assertions))]
const MT_BATCH_SIZE: usize = 128; // how many reads to chunk into a batch and push onto the queue for writing
#[cfg(not(debug_assertions))]
const MT_EPOCH_SIZE: usize = 512; // do this many batches at a time before writing to disk and clearing queues

// dev/debug, decrease values to better test multithreading
#[cfg(debug_assertions)]
const MT_BATCH_SIZE: usize = 17;
#[cfg(debug_assertions)]
const MT_EPOCH_SIZE: usize = 7;

fn main() -> io::Result<()> {
    let Args { reads1_path, reads2_path, barcodes_path, output_path, nthreads, compress_algo, compress_level, verbose } = Args::parse();

    // compress_level range depends on the compression algorithm
    let compress_suffix = match compress_algo.as_str() {
        "gz" | "gzip" => {
            if compress_level > 9 {
                panic!("For gzip compression, compress_level must be between 0 and 9");
            }
            "gz"
        },
        "zstd" => {
            if compress_level > 22 {
                panic!("For zstd compression, compress_level must be between 0 and 22");
            }
            "zst"
        },
        _ => unreachable!("Unexpected compression algorithm"),
    }.to_string();

    // check that verbose is 0 or 1
    if verbose > 1 {
        panic!("Verbosity level must be 0 or 1");
    }

    // Create the output directory if it doesn't exist.
    std::fs::create_dir_all(&output_path)
        .unwrap_or_else(|e| panic!("Failed to create output directory '{}': {}", output_path, e));

    // Open the input FASTQ files with proper decompression inferred from file extension.
    let records1 = records(open_input_decompressor(&reads1_path)?.lines());
    let records2 = records(open_input_decompressor(&reads2_path)?.lines());

    let sample_name1 = get_sample_name(&reads1_path)
        .unwrap_or_else(|e| panic!("Failed to extract sample name from '{}': {}", reads1_path, e));
    let sample_name2 = get_sample_name(&reads2_path)
        .unwrap_or_else(|e| panic!("Failed to extract sample name from '{}': {}", reads2_path, e));
    
    // Read barcodes from the file (one per line).
    let (barcode_names, barcode_sequences, max_distance, read1_dir, read1_matchpos, read2_dir, read2_matchpos) = read_barcode_file(&barcodes_path)
        .unwrap_or_else(|e| panic!("Failed to read barcode file '{}': {}", barcodes_path, e));
    let barcode_matcher = BarcodeMatcher::new(barcode_sequences, max_distance, read1_dir, read1_matchpos, read2_dir, read2_matchpos);

    // Pre-create output writers, one writer per barcode.
    let writer_map = get_writer_map(&barcode_names, &sample_name1, &sample_name2, &output_path, &compress_suffix, compress_level)
        .unwrap_or_else(|e| panic!("Failed to initialize writer map: {}", e));

    // Create writers for MULTIMAP and NOMAP.
    let multimap_writer = (
        open_output_compressor(&format!("{}/{}__MULTIMAP.fastq.{}", output_path, sample_name1, compress_suffix), &compress_algo, compress_level),
        open_output_compressor(&format!("{}/{}__MULTIMAP.fastq.{}", output_path, sample_name2, compress_suffix), &compress_algo, compress_level),
    );
    let nomap_writer = (
        open_output_compressor(&format!("{}/{}__NOMAP.fastq.{}", output_path, sample_name1, compress_suffix), &compress_algo, compress_level),
        open_output_compressor(&format!("{}/{}__NOMAP.fastq.{}", output_path, sample_name2, compress_suffix), &compress_algo, compress_level),
    );
    if nthreads == 1 {
        barcode_finder_single_thread(
            records1,
            records2,
            barcode_matcher,
            writer_map,
            multimap_writer,
            nomap_writer,
            verbose
        )
        .unwrap_or_else(|e| panic!("Error in single-threaded barcode finder: {}", e));
    } else {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(nthreads)
            .build()
            .unwrap_or_else(|e| panic!("Failed to build thread pool with {} threads: {}", nthreads, e));
        pool.install(move || {
            barcode_finder_multi_thread(
                records1,
                records2,
                barcode_matcher,
                writer_map,
                multimap_writer,
                nomap_writer,
                verbose
            )
            .unwrap_or_else(|e| panic!("Error in multi-threaded barcode finder: {}", e));
        });
    }
    Ok(())
}

// Helper function: convert an iterator of lines into an iterator of FASTQ records (4 lines per record).
fn records<I>(
    mut lines: I,
) -> impl Iterator<Item = io::Result<(String, String, String, String)>> + Send
where
    I: Iterator<Item = io::Result<String>> + Send,
{
    std::iter::from_fn(move || {
        let header = match lines.next()? {
            Ok(line) => line,
            Err(e) => return Some(Err(e)),
        };
        let seq = match lines.next()? {
            Ok(line) => line,
            Err(e) => return Some(Err(e)),
        };
        let plus = match lines.next()? {
            Ok(line) => line,
            Err(e) => return Some(Err(e)),
        };
        let qual = match lines.next()? {
            Ok(line) => line,
            Err(e) => return Some(Err(e)),
        };
        Some(Ok((header, seq, plus, qual)))
    })
}

fn get_sample_name(reads_path: &str) -> std::io::Result<String> {
    use std::path::Path;
    let sample_path = Path::new(reads_path);
    let filename = sample_path
        .file_name()
        .unwrap_or_else(|| panic!("Failed to extract file name from path '{}'", reads_path))
        .to_str()
        .unwrap_or_else(|| panic!("Failed to convert file name to UTF-8 from path '{}'", reads_path))
        .to_owned();
    let sample_name = if filename.ends_with(".fastq.gz") {
        filename.trim_end_matches(".fastq.gz").to_string()
    } else if filename.ends_with(".fq.gz") {
        filename.trim_end_matches(".fq.gz").to_string()
    } else if filename.ends_with(".fastq.zst") {
        filename.trim_end_matches(".fastq.zst").to_string()
    } else if filename.ends_with(".fq.zst") {
        filename.trim_end_matches(".fq.zst").to_string()
    } else if filename.ends_with(".fastq") {
        filename.trim_end_matches(".fastq").to_string()
    } else if filename.ends_with(".fq") {
        filename.trim_end_matches(".fq").to_string()
    } else {
        filename
    };
    Ok(sample_name)
}

fn get_writer_map(barcode_names: &Vec<String>, sample_name1: &str, sample_name2: &str, output_path: &str, compress_algo: &str, compress_level: u32) -> io::Result<Vec<(OutputCompressor, OutputCompressor)>> {
    let compress_suffix = match compress_algo {
        "gz" | "gzip" => "gz",
        "zstd" => "zst",
        _ => unreachable!("Unexpected compression algorithm"),
    };
    let mut writer_map = Vec::with_capacity(barcode_names.len());
    for name in barcode_names {
        let out_filename1 = format!("{}/{}__{}.fastq.{}", output_path, sample_name1, name, compress_suffix);
        let comp1 = open_output_compressor(&out_filename1, compress_algo, compress_level);
        let out_filename2 = format!("{}/{}__{}.fastq.{}", output_path, sample_name2, name, compress_suffix);
        let comp2 = open_output_compressor(&out_filename2, compress_algo, compress_level);
        writer_map.push((comp1, comp2));
    }
    Ok(writer_map)
}

// Helper function to write a FASTQ record.
fn write_record(
    writer: &mut (OutputCompressor, OutputCompressor),
    record1: &(String, String, String, String),
    record2: &(String, String, String, String),
) -> io::Result<()> {
    let (header, seq, plus, qual) = record1;
    writeln!(writer.0, "{}", header)
        .unwrap_or_else(|e| panic!("Failed to write header to output file (writer0): {}", e));
    writeln!(writer.0, "{}", seq)
        .unwrap_or_else(|e| panic!("Failed to write sequence to output file (writer0): {}", e));
    writeln!(writer.0, "{}", plus)
        .unwrap_or_else(|e| panic!("Failed to write '+' line to output file (writer0): {}", e));
    writeln!(writer.0, "{}", qual)
        .unwrap_or_else(|e| panic!("Failed to write quality string to output file (writer0): {}", e));
    let (header, seq, plus, qual) = record2;
    writeln!(writer.1, "{}", header)
        .unwrap_or_else(|e| panic!("Failed to write header to output file (writer1): {}", e));
    writeln!(writer.1, "{}", seq)
        .unwrap_or_else(|e| panic!("Failed to write sequence to output file (writer1): {}", e));
    writeln!(writer.1, "{}", plus)
        .unwrap_or_else(|e| panic!("Failed to write '+' line to output file (writer1): {}", e));
    writeln!(writer.1, "{}", qual)
        .unwrap_or_else(|e| panic!("Failed to write quality string to output file (writer1): {}", e));
    Ok(())
}

fn par_next_batch<I, J>(
    records1: &mut I,
    records2: &mut J,
) -> io::Result<Vec<((String, String, String, String), (String, String, String, String))>>
where
    I: Iterator<Item = io::Result<(String, String, String, String)>> + Send,
    J: Iterator<Item = io::Result<(String, String, String, String)>> + Send,
{
    use rayon::join;
    let (batch1_result, batch2_result) = join(
        || {
            let mut batch = Vec::with_capacity(MT_BATCH_SIZE);
            for _ in 0..MT_BATCH_SIZE {
                match records1.next() {
                    Some(Ok(record)) => batch.push(record),
                    Some(Err(e)) => return Err(e),
                    None => break,
                }
            }
            Ok(batch)
        },
        || {
            let mut batch = Vec::with_capacity(MT_BATCH_SIZE);
            for _ in 0..MT_BATCH_SIZE {
                match records2.next() {
                    Some(Ok(record)) => batch.push(record),
                    Some(Err(e)) => return Err(e),
                    None => break,
                }
            }
            Ok(batch)
        }
    );
    let batch1 = batch1_result?;
    let batch2 = batch2_result?;
    if batch1.len() != batch2.len() {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "Mismatched batch lengths",
        ));
    }
    Ok(batch1.into_iter().zip(batch2.into_iter()).collect())
}

fn barcode_finder_multi_thread(
    mut records1: impl Iterator<Item = io::Result<(String, String, String, String)>> + Send,
    mut records2: impl Iterator<Item = io::Result<(String, String, String, String)>> + Send,
    barcode_matcher: BarcodeMatcher,
    writer_map: Vec<(OutputCompressor, OutputCompressor)>,
    multimap_writer: (OutputCompressor, OutputCompressor),
    nomap_writer: (OutputCompressor, OutputCompressor),
    verbose: usize
) -> io::Result<()> {
    let mut barcode_writer = ConcurrentBarcodeWriter::<{MT_EPOCH_SIZE * MT_BATCH_SIZE}, _>::new(writer_map, nomap_writer, multimap_writer);
    let mut records_processed = 0;
    'outer: loop {
        let barcode_matcher_ref = &barcode_matcher;
        let barcode_writer_ref = &barcode_writer; // interior mutability
        let mut break_inner = false;
        let records_processed_ref = &mut records_processed;
        rayon::scope(|s| {
            for batch_number in 0..MT_EPOCH_SIZE {
                let batch = par_next_batch(&mut records1, &mut records2)
                    .unwrap_or_else(|e| panic!("Failed to get next batch of records: {}", e));
                if batch.is_empty() {
                    break_inner = true;
                    break;
                }
                *records_processed_ref += batch.len();
                if verbose > 0 {
                    println!("Processed {} records", *records_processed_ref);
                }
                s.spawn(move |_| {
                    let mut idx = batch_number * MT_BATCH_SIZE;
                    for (record1, record2) in batch.into_iter() {
                        match barcode_matcher_ref.find_unique_match(&record1.1, &record2.1) {
                            UniqueMatch::Unique(bc) => {
                                barcode_writer_ref.insert(bc, idx, (record1, record2));
                            }
                            UniqueMatch::Multiple => {
                                barcode_writer_ref.insert_multimap(idx, (record1, record2));
                            }
                            UniqueMatch::None => {
                                barcode_writer_ref.insert_nomap(idx, (record1, record2));
                            }
                        }
                        idx += 1;
                    }
                    
                });
            }
        });
        barcode_writer.par_write_data();
        if break_inner {
            break 'outer;
        }
    }
    barcode_writer.par_finish();
    Ok(())
}

fn barcode_finder_single_thread(
    mut records1: impl Iterator<Item = io::Result<(String, String, String, String)>>,
    mut records2: impl Iterator<Item = io::Result<(String, String, String, String)>>,
    barcode_matcher: BarcodeMatcher,
    mut writer_map: Vec<(OutputCompressor, OutputCompressor)>,
    mut multimap_writer: (OutputCompressor, OutputCompressor),
    mut nomap_writer: (OutputCompressor, OutputCompressor),
    verbose: usize
) -> io::Result<()> {
    let mut records_processed = 0;
    for (rec1, rec2) in (&mut records1).zip(&mut records2) {
        let record1 = rec1.unwrap_or_else(|e| panic!("Error reading record from file1 at record {}: {}", records_processed, e));
        let record2 = rec2.unwrap_or_else(|e| panic!("Error reading record from file2 at record {}: {}", records_processed, e));
        records_processed += 1;
        if verbose > 0 && records_processed % 1_000 == 0 {
            println!("Processed {} records", records_processed);
        }
        // Splitting by space and matching the first element should work for the vast majority of read header formats
        // https://en.wikipedia.org/wiki/FASTQ_format
        if record1.0.split_whitespace().next()
            .unwrap_or_else(|| panic!("Missing header in record from file1: '{}'", record1.0))
            != record2.0.split_whitespace().next()
            .unwrap_or_else(|| panic!("Missing header in record from file2: '{}'", record2.0))
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Headers do not match: {} != {}", record1.0, record2.0),
            ));
        }
        match barcode_matcher.find_unique_match(&record1.1, &record2.1) {
            UniqueMatch::Unique(bc) => {
                write_record(&mut writer_map[bc], &record1, &record2)
                    .unwrap_or_else(|e| panic!("Failed to write record for barcode index {}: {}", bc, e));
            }
            UniqueMatch::Multiple => {
                write_record(&mut multimap_writer, &record1, &record2)
                    .unwrap_or_else(|e| panic!("Failed to write record to MULTIMAP writer: {}", e));
            }
            UniqueMatch::None => {
                write_record(&mut nomap_writer, &record1, &record2)
                    .unwrap_or_else(|e| panic!("Failed to write record to NOMAP writer: {}", e));
            }
        }        
    }
    for (writer1, writer2) in writer_map.into_iter() {
        writer1.finish().unwrap_or_else(|e| panic!("Failed to finish barcode writer (writer1): {}", e));
        writer2.finish().unwrap_or_else(|e| panic!("Failed to finish barcode writer (writer2): {}", e));
    }
    multimap_writer.0.finish().unwrap_or_else(|e| panic!("Failed to finish multimap writer (writer0): {}", e));
    multimap_writer.1.finish().unwrap_or_else(|e| panic!("Failed to finish multimap writer (writer1): {}", e));
    nomap_writer.0.finish().unwrap_or_else(|e| panic!("Failed to finish nomap writer (writer0): {}", e));
    nomap_writer.1.finish().unwrap_or_else(|e| panic!("Failed to finish nomap writer (writer1): {}", e));
    if records1.next().is_some() {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "File 1 has more records than file 2",
        ));
    }
    if records2.next().is_some() {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "File 2 has more records than file 1",
        ));
    }
    Ok(())
}
