// extern crate flate2;
// extern crate zstd;

use std::fs::File;
use std::io::{self, BufReader, Read, BufRead};

use flate2::read::MultiGzDecoder;

/// Enum for input decompression with static dispatch.
/// The Gzip variant wraps a MultiGzDecoder<File> inside a BufReader,
// NB: Zstd implementation is kind of funky
// It specifies BufReader<File> as the inner stream type even though it takes in JUST FILE as a parameter
// It also requires a lifetime trait, even though BufReader<File> is owned by the Decoder
pub enum InputDecompressor {
    Gzip(BufReader<MultiGzDecoder<File>>),
    Zstd(BufReader<zstd::stream::Decoder<'static, BufReader<File>>>),
}

impl Read for InputDecompressor {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            InputDecompressor::Gzip(reader) => reader.read(buf),
            InputDecompressor::Zstd(decoder) => decoder.read(buf),
        }
    }
}

impl BufRead for InputDecompressor {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        match self {
            InputDecompressor::Gzip(reader) => reader.fill_buf(),
            InputDecompressor::Zstd(decoder) => decoder.fill_buf(),
        }
    }
    fn consume(&mut self, amt: usize) {
        match self {
            InputDecompressor::Gzip(reader) => reader.consume(amt),
            InputDecompressor::Zstd(decoder) => decoder.consume(amt),
        }
    }
}

/// Open an input file and return a decompressor inferred strictly from the file extension.
pub fn open_input_decompressor(path: &str) -> io::Result<InputDecompressor> {
    if path.ends_with(".zst") {
        let file = File::open(path)
            .unwrap_or_else(|e| panic!("Failed to open input file '{}': {}", path, e));
        let decoder = zstd::stream::Decoder::new(file)
            .unwrap_or_else(|e| panic!("Failed to create zstd decoder for file '{}': {}", path, e));
        Ok(InputDecompressor::Zstd(BufReader::new(decoder)))
    } else if path.ends_with(".gz") {
        let file = File::open(path)
            .unwrap_or_else(|e| panic!("Failed to open input file '{}': {}", path, e));
        let decoder = MultiGzDecoder::new(file);
        Ok(InputDecompressor::Gzip(BufReader::new(decoder)))
    } else {
        panic!("Unsupported file format for input: {}", path);
    }
}
