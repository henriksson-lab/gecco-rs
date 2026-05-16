//! Shared utility helpers, mirroring the small grab-bag in
//! `gecco/_meta.py`.

use std::ops::Range;

/// Iterate over a sequence of length `length` with a sliding window of
/// width `window` and stride `step`.
///
/// # Panics
///
/// Panics if `window` is zero, or if `step` is zero or larger than `window`.
pub fn sliding_window(length: usize, window: usize, step: usize) -> Vec<Range<usize>> {
    assert!(window > 0, "Window size must be strictly positive");
    assert!(
        step > 0 && step <= window,
        "Window step must be strictly positive and <= window_size"
    );
    let mut windows = Vec::new();
    let mut i = 0;
    while i + window <= length {
        windows.push(i..i + window);
        i += step;
    }
    windows
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sliding_window() {
        let wins = sliding_window(10, 3, 1);
        assert_eq!(wins.len(), 8);
        assert_eq!(wins[0], 0..3);
        assert_eq!(wins[7], 7..10);
    }

    #[test]
    fn test_sliding_window_step2() {
        let wins = sliding_window(10, 3, 2);
        assert_eq!(wins.len(), 4);
        assert_eq!(wins[0], 0..3);
        assert_eq!(wins[3], 6..9);
    }
}
