use bio::pattern_matching::myers::Myers;

/// An enum representing the allowed read directions.
// Defined individually for Read1 or Read2
#[derive(Debug, Clone, Copy)]
pub enum ReadDir {
    Forward,
    ReverseComplement,
    Either,
    NoMatchAllowed,
}

/// Allowed match positions.
#[derive(Debug, Clone, Copy)]
pub enum MatchPos {
    Start,
    End,
    Anywhere,
    NoMatchAllowed,
}

/// Returns the reverse complement of the given DNA string.
fn reverse_complement(dna: &str) -> String {
    dna.chars()
       .rev()
       .map(|nucleotide| match nucleotide {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            other => other,
       })
       .collect()
}

/// An enum representing the match status for a barcode.
#[derive(Debug)]
pub enum UniqueMatch {
    /// Exactly one barcode matched, with the given index.
    Unique(usize),
    /// No barcode matched.
    None,
    /// More than one barcode matched.
    Multiple,
}

/// BarcodeMatcher holds the parsed barcode data (without names):
/// - fw_patterns: Myers patterns built from the forward barcode sequences,
/// - rc_patterns: Myers patterns built from the reverse complement sequences,
/// - max_distances: maximum allowed edit distance for each barcode,
/// - read1_dirs: allowed read direction for read1,
/// - read1_matchpos: allowed match position for read1,
/// - read2_dirs: allowed read direction for read2.
/// - read2_matchpos: allowed match position for read2.
pub struct BarcodeMatcher {
    pub fw_patterns: Vec<Myers<u64>>,
    pub rc_patterns: Vec<Myers<u64>>,
    pub max_distances: Vec<u8>,
    pub read1_dirs: Vec<ReadDir>,
    pub read1_matchpos: Vec<MatchPos>,
    pub read2_dirs: Vec<ReadDir>,
    pub read2_matchpos: Vec<MatchPos>,
}

impl BarcodeMatcher {
    /// Constructs a BarcodeMatcher from four vectors.
    /// `sequences` is used to build the forward patterns and, via reverse complement, the reverse patterns.
    pub fn new(
        sequences: Vec<String>,
        max_distances: Vec<u8>,
        read1_dirs: Vec<ReadDir>,
        read1_matchpos: Vec<MatchPos>,
        read2_dirs: Vec<ReadDir>,
        read2_matchpos: Vec<MatchPos>,

    ) -> Self {
        let rc_sequences: Vec<String> = sequences.iter().map(|s| reverse_complement(s)).collect();
        let fw_patterns: Vec<Myers<u64>> = sequences
            .iter()
            .map(|s| Myers::<u64>::new(s.as_bytes()))
            .collect();
        let rc_patterns: Vec<Myers<u64>> = rc_sequences
            .iter()
            .map(|s| Myers::<u64>::new(s.as_bytes()))
            .collect();
        Self {
            fw_patterns,
            rc_patterns,
            max_distances,
            read1_dirs,
            read1_matchpos,
            read2_dirs,
            read2_matchpos,
        }
    }

    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.fw_patterns.len()
    }

    /// Given two reads, returns the match status as follows:
    ///
    /// - If exactly one barcode matches (based on its stored max_distance and allowed directions),
    ///   returns `UniqueMatch::Unique(index)`.
    /// - If no barcode matches, returns `UniqueMatch::None`.
    /// - If more than one barcode matches, returns `UniqueMatch::Multiple`.
    ///
    /// For a barcode at index `i`:
    /// - If `read1_dirs[i]` is Fwd or Both, read1 is matched against the forward pattern.
    /// - If `read1_dirs[i]` is Rev or Both, read1 is matched against the reverse pattern.
    /// - Similarly for read2 with `read2_dirs[i]`.
    /// Matching is performed with bio::pattern_matching::myers::Myers. If MatchPos is
    /// Start, the barcode is matched at the start of the read. Since Myers does not anchor
    /// patterns, the match is considered valid if the edit distance plus the start position
    /// is less than or equal to the max_distance. This is valid as Myers outputs all possible 
    /// alignments, so it is not favoring alignments that are farther away from the start. 
    /// Similar logic if MatchPos is End. If MatchPos is Anywhere, the barcode can be
    /// anywhere in the read.
    /// Implementation details: for finding matches at the start of the read, we use reverse_complement 
    /// of the read and then use find_all_end. This is because Myers::find_all requires 
    /// mutability, in order to perform an internal traceback, therefore cannot be called from
    /// a multithreaded context. Also Myers::find_all_end should be more efficient than 
    /// Myers::find_all as it does not require traceback. 
    pub fn find_unique_match(&self, read1: &str, read2: &str) -> UniqueMatch {
        let match_in_read = |read: &str,
                             fw_pat: &Myers<u64>,
                             rc_pat: &Myers<u64>,
                             max_distance: u8,
                             dir: &ReadDir,
                             bp: &MatchPos|
             -> bool {
                match bp {
                    MatchPos::Start => {
                        // Reverse complement the read so that the original beginning becomes
                        // the end of the reversed read.
                        let rev_read = reverse_complement(read);
                        // In the reversed context, the roles swap:
                        // - If the allowed direction is Fwd (or Both) in the original,
                        //   use the reverse-complement pattern (rc_pat).
                        if matches!(dir, ReadDir::Forward) || matches!(dir, ReadDir::Either) {
                            if rc_pat.find_all_end(rev_read.as_bytes(), max_distance)
                                     .any(|(end, dist)| (rev_read.len() + dist as usize) <= (max_distance as usize + end + 1)) {
                                return true;
                            }
                        }
                        // - If the allowed direction is Rev (or Both) in the original,
                        //   use the forward pattern (fw_pat).
                        if matches!(dir, ReadDir::ReverseComplement) || matches!(dir, ReadDir::Either) {
                            if fw_pat.find_all_end(rev_read.as_bytes(), max_distance)
                                     .any(|(end, dist)| (rev_read.len() + dist as usize) <= (max_distance as usize + end + 1)) {
                                return true;
                            }
                        }
                        false // No match found
                    },
                    MatchPos::End => {
                        if matches!(dir, ReadDir::Forward) || matches!(dir, ReadDir::Either) {
                            if fw_pat.find_all_end(read.as_bytes(), max_distance)
                                     .any(|(end, dist)| (read.len() + dist as usize) <= (max_distance as usize + end + 1)) {
                                return true;
                            }
                        }
                        if matches!(dir, ReadDir::ReverseComplement) || matches!(dir, ReadDir::Either) {
                            if rc_pat.find_all_end(read.as_bytes(), max_distance)
                                     .any(|(end, dist)| (read.len() + dist as usize) <= (max_distance as usize + end + 1)) {
                                return true;
                            }
                        }
                        false // No match found
                    },
                    MatchPos::Anywhere => {
                        // For non-anchored matches, use find_all_end (immutable) to check if any match exists.
                        if matches!(dir, ReadDir::Forward) || matches!(dir, ReadDir::Either) {
                            if fw_pat.find_all_end(read.as_bytes(), max_distance).next().is_some() {
                                return true;
                            }
                        }
                        if matches!(dir, ReadDir::ReverseComplement) || matches!(dir, ReadDir::Either) {
                            if rc_pat.find_all_end(read.as_bytes(), max_distance).next().is_some() {
                                return true;
                            }
                        }
                        false // No match found
                    },
                    MatchPos::NoMatchAllowed => false,
                }
            };
    
        let mut candidate: Option<usize> = None;
        for i in 0..self.fw_patterns.len() {
            let fw_pat = &self.fw_patterns[i];
            let rc_pat = &self.rc_patterns[i];
            let maxdis = self.max_distances[i];
            let rd1 = &self.read1_dirs[i];
            let mp1 = &self.read1_matchpos[i];
            let rd2 = &self.read2_dirs[i];
            let mp2 = &self.read2_matchpos[i];

            // |read: &str,
            // fw_pat: &Myers<u64>,
            // rc_pat: &Myers<u64>,
            // max_distance: u8,
            // dir: &ReadDir,
            // bp: &MatchPos|
            let is_match = match_in_read(read1, fw_pat, rc_pat, maxdis, rd1, mp1)
                        || match_in_read(read2, fw_pat, rc_pat, maxdis, rd2, mp2);
    
            if is_match {
                if candidate.is_some() {
                    return UniqueMatch::Multiple; // if a match was already found to a different barcode
                } else {
                    candidate = Some(i);
                }
            }
        }
        match candidate {
            Some(i) => UniqueMatch::Unique(i),
            None => UniqueMatch::None,
        }
    }
}
