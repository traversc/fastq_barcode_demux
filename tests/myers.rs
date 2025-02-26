/// This is an example of the difference in treatment between Start and End positions 
/// Using Myers alignment, as of rust-bio 2.0.3. End positions are not interchangeable 
/// with start positions after reverse complement.

use bio::pattern_matching::myers::Myers;

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

#[test]
fn main() -> std::io::Result<()> {
    let bc = "TGATTGTCTTTTAAACGTA".to_string();
    let mutbc = "TIGATTGTCTTTTAAACGTAI".to_string();
    let mut pattern = Myers::<u64>::new(bc.as_bytes());
    let occ: Vec<_> = pattern.find_all(mutbc.as_bytes(), 2).collect();
    println!("find_all");
    for (start, end, dist) in occ {
        println!("Start {} End {} Dist {} Mut len {}", start, end, dist, mutbc.len());
    }
    println!("find_all_end");
    let occ: Vec<_> = pattern.find_all_end(mutbc.as_bytes(), 2).collect();
    for (end, dist) in occ {
        println!("End {} Dist {} Mut len {}", end, dist, mutbc.len());
    }

    let bc_rc = revcomp(&bc);
    let mutbc_rc = revcomp(&mutbc);
    let mut pattern = Myers::<u64>::new(bc_rc.as_bytes());
    let occ: Vec<_> = pattern.find_all(mutbc_rc.as_bytes(), 2).collect();
    println!("find_all_rc");
    for (start, end, dist) in occ {
        println!("Start {} End {} Dist {} Mut len {}", start, end, dist, mutbc.len());
    }
    println!("find_all_end_rc");
    let occ: Vec<_> = pattern.find_all_end(mutbc_rc.as_bytes(), 2).collect();
    for (end, dist) in occ {
        println!("End {} Dist {} Mut len {}", end, dist, mutbc.len());
    }
    Ok(())
}