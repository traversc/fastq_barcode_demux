use rayon::prelude::*;
use rayon::iter::{Either, ParallelIterator};
use std::cell::UnsafeCell;
use std::io::{self, Write};

use crate::common::output_compressor_traits::Finish;

type RECORD = (String, String, String, String);

pub struct ConcurrentBarcodeWriter<const MAXP: usize, WriterType> {
    pub barcode_queue: Vec<Vec<UnsafeCell<(RECORD, RECORD)>>>,
    pub nomap_queue: Vec<UnsafeCell<(RECORD, RECORD)>>,
    pub multimap_queue: Vec<UnsafeCell<(RECORD, RECORD)>>,
    pub barcode_writer: Vec<(WriterType, WriterType)>,
    pub nomap_writer: (WriterType, WriterType),
    pub multimap_writer: (WriterType, WriterType),
    pub output_nomap: bool,
}

impl<const MAXP: usize, WriterType> ConcurrentBarcodeWriter<MAXP, WriterType>
where
    WriterType: Write + Send,
{
    fn write_record(
        writer: &mut (WriterType, WriterType),
        record1: &RECORD,
        record2: &RECORD,
    ) -> io::Result<()> {
        let (header, seq, plus, qual) = record1;
        writeln!(writer.0, "{}", header)?;
        writeln!(writer.0, "{}", seq)?;
        writeln!(writer.0, "{}", plus)?;
        writeln!(writer.0, "{}", qual)?;
        let (header, seq, plus, qual) = record2;
        writeln!(writer.1, "{}", header)?;
        writeln!(writer.1, "{}", seq)?;
        writeln!(writer.1, "{}", plus)?;
        writeln!(writer.1, "{}", qual)?;
        Ok(())
    }

    pub fn new(
        barcode_writer: Vec<(WriterType, WriterType)>,
        nomap_writer: (WriterType, WriterType),
        multimap_writer: (WriterType, WriterType),
        output_nomap: bool,
    ) -> Self {
        let barcode_queue: Vec<Vec<UnsafeCell<(RECORD, RECORD)>>> = barcode_writer
            .iter()
            .map(|_| {
                std::iter::repeat_with(|| {
                    UnsafeCell::new((RECORD::default(), RECORD::default()))
                })
                .take(MAXP)
                .collect()
            })
            .collect();

        let nomap_queue: Vec<UnsafeCell<(RECORD, RECORD)>> = std::iter::repeat_with(|| {
            UnsafeCell::new((RECORD::default(), RECORD::default()))
        })
        .take(MAXP)
        .collect();

        let multimap_queue: Vec<UnsafeCell<(RECORD, RECORD)>> = std::iter::repeat_with(|| {
            UnsafeCell::new((RECORD::default(), RECORD::default()))
        })
        .take(MAXP)
        .collect();

        Self {
            barcode_queue,
            nomap_queue,
            multimap_queue,
            barcode_writer,
            nomap_writer,
            multimap_writer,
            output_nomap,
        }
    }

    pub fn insert(&self, i: usize, j: usize, x: (RECORD, RECORD)) {
        unsafe {
            *self.barcode_queue.get_unchecked(i).get_unchecked(j).get() = x;
        }
    }

    pub fn insert_nomap(&self, j: usize, x: (RECORD, RECORD)) {
        unsafe {
            *self.nomap_queue.get_unchecked(j).get() = x;
        }
    }

    pub fn insert_multimap(&self, j: usize, x: (RECORD, RECORD)) {
        unsafe {
            *self.multimap_queue.get_unchecked(j).get() = x;
        }
    }

    pub fn par_write_data(&mut self) {
        // RECORD and QueueItem are fine because they don't refer to WriterType
        type RECORD = (String, String, String, String);
        type QueueItem = std::cell::UnsafeCell<(RECORD, RECORD)>;
    
        // Build the *optional* nomap‐slice iterator
        let nomap_queue_iter = if self.output_nomap {
            Either::Left(rayon::iter::once(self.nomap_queue.as_mut_slice()))
        } else {
            Either::Right(rayon::iter::empty::<&mut [QueueItem]>())
        };
    
        let queue_iter = self
            .barcode_queue
            .par_iter_mut()
            .map(|v| v.as_mut_slice())
            .chain(rayon::iter::once(self.multimap_queue.as_mut_slice()))
            .chain(nomap_queue_iter);
    
        // Build the *optional* nomap‐writer iterator, using the generic WriterType directly
        let nomap_writer_iter = if self.output_nomap {
            Either::Left(rayon::iter::once(&mut self.nomap_writer))
        } else {
            Either::Right(
                rayon::iter::empty::<&mut (WriterType, WriterType)>(),
            )
        };
    
        let writer_iter = self
            .barcode_writer
            .par_iter_mut()
            .chain(rayon::iter::once(&mut self.multimap_writer))
            .chain(nomap_writer_iter);
    
        queue_iter.zip(writer_iter).for_each(|(queue, writer)| {
            for cell in queue.iter() {
                let (r1, r2) = unsafe { &*cell.get() };
                if !r1.0.is_empty() {
                    Self::write_record(writer, r1, r2).unwrap();
                }
                unsafe {
                    *cell.get() = (RECORD::default(), RECORD::default());
                }
            }
        });
    }    
    
    
    pub fn par_finish(self)
    where
        WriterType: Finish + Send,
    {
        self.barcode_writer
            .into_par_iter()
            .chain(rayon::iter::once(self.nomap_writer))
            .chain(rayon::iter::once(self.multimap_writer))
            .for_each(|(w1, w2)| {
                let _ = w1.finish().unwrap();
                let _ = w2.finish().unwrap();
            });
    }
}

unsafe impl<const MAXP: usize, WriterType> Sync for ConcurrentBarcodeWriter<MAXP, WriterType> {}
