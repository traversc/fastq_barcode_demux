// extern statements not necessary in rust 2021
// extern crate flate2;
// extern crate zstd;

use std::fs::File;
use std::io::{self, Write, BufWriter};

use flate2::write::GzEncoder;
use flate2::Compression;


/// A trait to abstract finishing the compression.
pub trait Finish {
    type Output;
    fn finish(self) -> io::Result<Self::Output>;
}

/// Enum for output compression with static dispatch.
/// The Gzip variant wraps a GzEncoder<File> inside a BufReader,
/// Zstd variant wraps a zstd::stream::write::Encoder<File> inside a BufReader.
// Zstd requires lifetime parameter
pub enum OutputCompressor {
    Gzip(GzEncoder<BufWriter<File>>),
    Zstd(zstd::stream::write::Encoder<'static, BufWriter<File>>),
}

impl Write for OutputCompressor {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
         match self {
             OutputCompressor::Gzip(encoder) => encoder.write(buf),
             OutputCompressor::Zstd(encoder) => encoder.write(buf),
         }
    }
    fn flush(&mut self) -> io::Result<()> {
         match self {
             OutputCompressor::Gzip(encoder) => encoder.flush(),
             OutputCompressor::Zstd(encoder) => encoder.flush(),
         }
    }
}

impl Finish for OutputCompressor {
    type Output = ();
    fn finish(self) -> io::Result<Self::Output> {
        match self {
            OutputCompressor::Gzip(encoder) => encoder.finish().map(|_| ()),
            OutputCompressor::Zstd(encoder) => encoder.finish().map(|_| ()),
        }
    }
}

/// Open an output file using the provided compression algorithm.
/// The compressor type is determined solely by the `compress_algo` argument.
/// Valid values are "gz", "gzip", and "zstd".
/// For zstd, the compress_level is an i32, and for gzip it's a u32.
pub fn open_output_compressor(path: &str, compress_algo: &str, compress_level: u32) -> OutputCompressor {
    let file = File::create(path)
        .unwrap_or_else(|e| panic!("Failed to create output file '{}': {}", path, e));
    let buf_writer = BufWriter::new(file);
    match compress_algo {
        "zstd" => {
            let encoder = zstd::stream::write::Encoder::new(buf_writer, compress_level as i32) // 0 defaults to 3
                .unwrap_or_else(|e| panic!("Failed to create zstd encoder for file '{}': {}", path, e));
            OutputCompressor::Zstd(encoder)
        },
        "gz" | "gzip" | _ => {
            let compress_level_gzip = match compress_level {
                0 => Compression::default(),
                1..=9 => Compression::new(compress_level),
                _ => panic!("For gzip compression, compress_level must be between 0 and 9"),
            };
            let encoder = GzEncoder::new(buf_writer, compress_level_gzip);
            OutputCompressor::Gzip(encoder)
        },
    }
}
