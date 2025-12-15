//! Module for utilities

/// Check if the num is a power of 2 or not
///
/// # Arguments
///
/// * `n` - Number of positive integer
///
pub fn is_power_of_two(n: usize) -> bool {
    n > 0 && (n & (n - 1)) == 0
}

/// Calculate standard deviation
///
/// # Arguments
///
/// * `data` - data sequence
///
pub fn std_dev(data: &[f64]) -> Option<f64> {
    if data.is_empty() {
        return None;
    }

    // mean value
    let sum: f64 = data.iter().sum();
    let mean = sum / data.len() as f64;

    // variance
    let variance: f64 = data
        .iter()
        .map(|&x| {
            let diff = x - mean;
            diff * diff
        })
        .sum::<f64>()
        / data.len() as f64;

    Some(variance.sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::approx_eq;

    #[test]
    fn std_dev_success() {
        let data: Vec<f64> = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        assert!(approx_eq!(
            f64,
            std_dev(&data).unwrap(),
            2.0_f64.sqrt(),
            epsilon = 1.0e-12
        ));
    }
}
