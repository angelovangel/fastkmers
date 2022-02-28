
pub fn prob_am(q: &[u8]) -> f64 {
    let mut prob_arithmetic_mean = 0.0;
    for &item in q {
        let phred = *&item as f64 - 33.0;
        let prob = 10.0_f64.powf(-phred / 10.0);
        prob_arithmetic_mean += prob;
    }
    prob_arithmetic_mean / q.len() as f64
}

// use log method to avoid overflow
// https://www.geeksforgeeks.org/geometric-mean-two-methods/

pub fn phred_gm(q: &[u8]) -> f64 {
    let mut phred_product = 0.;
    
    for &item in q {                       
    let phred = *&item as f64 - 33.0;

            phred_product += phred.ln();

    }
    let result = phred_product / q.len() as f64;
    result.exp()
}