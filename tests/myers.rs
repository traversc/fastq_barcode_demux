use bio::pattern_matching::myers::Myers;

#[test]
fn main() -> std::io::Result<()> {
    let bc = "ACTGACTGACTGACTACGT".to_string();
    let mutbc = "ACTGACTGACTGACTACGTI".to_string();
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
    Ok(())
}