//! # Chart Trends
//!
//! The `chart_trends` module provides utilities for detecting, analyzing, and breaking down trends in price charts.
//! These functions help identify overall direction, peaks, valleys, and trend segments in a time series.
//!
//! ## When to Use
//! Use chart trend indicators when you want to:
//! - Decompose a price series into upward/downward trends
//! - Find peaks and valleys for support/resistance analysis
//! - Quantify the overall or local trend direction of an asset
//!
//! ## Structure
//! Unlike other modules, `chart_trends` does not have `single` or `bulk` submodules.
//! All functions operate over slices and return either trend breakdowns or locations of key points.
//!
//! ## Included Functions
//! - [`break_down_trends`]: Segments the chart into distinct up/down trends
//! - [`overall_trend`]: Returns the overall trend (slope) for all price points
//! - [`peak_trend`]: Calculates the trend based on local peaks
//! - [`peaks`]: Finds all local maxima (peaks) in the series
//! - [`valley_trend`]: Calculates the trend based on local valleys
//! - [`valleys`]: Finds all local minima (valleys) in the series
//!
//! ## API Details
//! - All functions work on slices of `f64` prices (or equivalent).
//! - Returns are typically vectors of trend segments or indices/values of peaks/valleys.
//! - See each function's documentation for examples, panics, and usage tips.
//!
//! ---

use crate::basic_indicators::single::{max, mean, min};

/// Calculates all peaks over a given period
///
/// # Arguments
///
/// * `prices` - Slice of prices
/// * `period` - Period over which to find the peak
/// * `closest_neighbor` - Minimum distance between peaks
///
/// # Panics
///
/// Panics if:
///     * `period` == 0
///     * `period` > `prices.len()`
///
/// # Examples
///
/// ```rust
/// let highs = vec![103.0, 102.0, 107.0, 104.0, 100.0];
/// let period: usize = 3;
/// let closest_neighbor: usize = 1;
/// let peaks = rust_ti::chart_trends::peaks(&highs, period, closest_neighbor);
/// assert_eq!(vec![(107.0, 2)], peaks);
///
/// let highs = vec![103.0, 102.0, 107.0, 104.0, 100.0, 109.0];
/// let period: usize = 3;
/// let peaks = rust_ti::chart_trends::peaks(&highs, period, closest_neighbor);
/// assert_eq!(vec![(107.0, 2), (109.0, 5)], peaks);
///
/// let highs = vec![103.0, 102.0, 107.0, 104.0, 100.0, 109.0];
/// let period: usize = 6;
/// let peaks = rust_ti::chart_trends::peaks(&highs, period, closest_neighbor);
/// assert_eq!(vec![(109.0, 5)], peaks);
///
/// let highs = vec![103.0, 102.0, 107.0, 104.0, 100.0, 107.0];
/// let period: usize = 3;
/// let peaks = rust_ti::chart_trends::peaks(&highs, period, closest_neighbor);
/// assert_eq!(vec![(107.0, 2), (107.0, 5)], peaks);
///
/// // If there are 2 peaks it will take the most recent one
/// let highs = vec![103.0, 102.0, 107.0, 104.0, 100.0, 107.0];
/// let period: usize = 6;
/// let peaks = rust_ti::chart_trends::peaks(&highs, period, closest_neighbor);
/// assert_eq!(vec![(107.0, 5)], peaks);
/// ```
pub fn peaks(prices: &[f64], period: usize, closest_neighbor: usize) -> Vec<(f64, usize)> {
    if period == 0 {
        panic!("Period ({}) must be greater than 0", period)
    };
    let length = prices.len();
    if period > length {
        panic!(
            "Period ({}) cannot be longer than length of prices ({})",
            period, length
        )
    };

    let mut peaks: Vec<(f64, usize)> = Vec::new();
    let mut last_peak_idx: usize = 0;
    let mut last_peak: f64 = 0.0;

    for i in 0..=length - period {
        let window = &prices[i..i + period];
        let peak = max(window);
        let local_idx = window.iter().rposition(|&x| x == peak).unwrap();
        let idx = i + local_idx;

        if last_peak_idx != 0 {
            if idx <= last_peak_idx + closest_neighbor {
                if peak < last_peak {
                    last_peak_idx = idx;
                } else if peak > last_peak {
                    peaks.pop();
                    peaks.push((peak, idx));
                    last_peak_idx = idx;
                    last_peak = peak;
                }
            } else if !peaks.contains(&(peak, idx)) {
                peaks.push((peak, idx));
                last_peak_idx = idx;
                last_peak = peak;
            }
        } else {
            peaks.push((peak, idx));
            last_peak_idx = idx;
            last_peak = peak;
        }
    }
    peaks
}

/// Calculates all valleys for a given period.
///
/// # Arguments
///
/// * `prices` - Slice of prices
/// * `period` - Period over which to find the valley
/// * `closest_neighbor` - Minimum distance between valleys
///
/// # Panics
///
/// Panics if:
///     * `period` == 0
///     * `period` > `prices.len()`
///
/// # Examples
///
/// ```rust
/// let lows = vec![98.0, 101.0, 95.0, 100.0, 97.0];
/// let period: usize = 3;
/// let closest_neighbor: usize = 1;
/// let valleys = rust_ti::chart_trends::valleys(&lows, period, closest_neighbor);
/// assert_eq!(vec![(95.0, 2)], valleys);
///
/// let lows = vec![98.0, 101.0, 95.0, 100.0, 97.0, 93.0];
/// let period: usize = 3;
/// let valleys = rust_ti::chart_trends::valleys(&lows, period, closest_neighbor);
/// assert_eq!(vec![(95.0, 2), (93.0, 5)], valleys);
///
/// let lows = vec![98.0, 101.0, 95.0, 100.0, 97.0, 93.0];
/// let period: usize = 6;
/// let valleys = rust_ti::chart_trends::valleys(&lows, period, closest_neighbor);
/// assert_eq!(vec![(93.0, 5)], valleys);
///
/// let lows = vec![98.0, 101.0, 95.0, 100.0, 97.0, 95.0];
/// let period: usize = 3;
/// let valleys = rust_ti::chart_trends::valleys(&lows, period, closest_neighbor);
/// assert_eq!(vec![(95.0, 2), (95.0, 5)], valleys);
///
/// let lows = vec![98.0, 101.0, 95.0, 100.0, 97.0, 95.0];
/// let period: usize = 6;
/// let valleys = rust_ti::chart_trends::valleys(&lows, period, closest_neighbor);
/// assert_eq!(vec![(95.0, 5)], valleys);
/// ```
pub fn valleys(prices: &[f64], period: usize, closest_neighbor: usize) -> Vec<(f64, usize)> {
    if period == 0 {
        panic!("Period ({}) must be greater than 0", period)
    };
    let length = prices.len();
    if period > length {
        panic!(
            "Period ({}) cannot be longer than length of prices ({})",
            period, length
        )
    };
    let mut peaks: Vec<(f64, usize)> = Vec::new();
    let mut last_peak_idx: usize = 0;
    let mut last_peak: f64 = 0.0;

    for i in 0..=length - period {
        let window = &prices[i..i + period];
        let peak = max(window);
        let local_idx = window.iter().rposition(|&x| x == peak).unwrap();
        let idx = i + local_idx;

        if last_peak_idx != 0 {
            if idx <= last_peak_idx + closest_neighbor {
                if peak < last_peak {
                    last_peak_idx = idx;
                } else if peak > last_peak {
                    peaks.pop();
                    peaks.push((peak, idx));
                    last_peak_idx = idx;
                    last_peak = peak;
                }
            } else if !peaks.contains(&(peak, idx)) {
                peaks.push((peak, idx));
                last_peak_idx = idx;
                last_peak = peak;
            }
        } else {
            peaks.push((peak, idx));
            last_peak_idx = idx;
            last_peak = peak;
        }
    }
    let mut valleys: Vec<(f64, usize)> = Vec::new();
    let mut last_valley_idx: usize = 0;
    let mut last_valley: f64 = 0.0;

    for i in 0..=length - period {
        let window = &prices[i..i + period];
        let valley = min(window);
        let local_idx = window.iter().rposition(|&x| x == valley).unwrap();
        let idx = i + local_idx;

        if last_valley_idx != 0 {
            if idx <= last_valley_idx + closest_neighbor {
                if valley > last_valley {
                    last_valley_idx = idx;
                } else if valley < last_valley {
                    valleys.pop();
                    valleys.push((valley, idx));
                    last_valley_idx = idx;
                    last_valley = valley;
                }
            } else if !valleys.contains(&(valley, idx)) {
                valleys.push((valley, idx));
                last_valley_idx = idx;
                last_valley = valley;
            }
        } else {
            valleys.push((valley, idx));
            last_valley_idx = idx;
            last_valley = valley;
        }
    }
    valleys
}

/// OLS simple linear regression function
fn get_trend_line(p: &[(f64, usize)]) -> (f64, f64) {
    let length = p.len() as f64;
    let mean_x = p.iter().map(|&(_, x)| x as f64).sum::<f64>() / length;
    let mean_y = p.iter().map(|&(y, _)| y).sum::<f64>() / length;

    let (num, den) = p.iter().fold((0.0, 0.0), |(num, den), &(y, x)| {
        let x = x as f64;
        let dx = x - mean_x;
        (num + dx * (y - mean_y), den + dx * dx)
    });

    let slope = num / den;
    let intercept = mean_y - (slope * mean_x);
    (slope, intercept)
}

/// Returns the slope and intercept of the trend line fitted to peaks.
///
/// # Arguments
///
/// * `prices` - Slice of prices
/// * `period` - Period over which to calculate the peaks
///
/// # Examples
///
/// ```rust
/// let highs = vec![103.0, 102.0, 107.0, 104.0, 100.0, 109.0];
/// let period: usize = 3;
/// let peak_trend = rust_ti::chart_trends::peak_trend(&highs, period);
/// assert_eq!((0.6666666666666666, 105.66666666666667), peak_trend);
/// ```
pub fn peak_trend(prices: &[f64], period: usize) -> (f64, f64) {
    let peaks = peaks(prices, period, 1);
    get_trend_line(&peaks)
}

/// Calculates the slope and intercept of the trend line fitted to valleys.
///
/// # Arguments
///
/// * `prices` - Slice of prices
/// * `period` - Period over which to calculate the valleys
///
/// # Examples
///
/// ```rust
/// let lows = vec![98.0, 101.0, 95.0, 100.0, 97.0, 93.0];
/// let period: usize = 3;
/// let valley_trend = rust_ti::chart_trends::valley_trend(&lows, period);
/// assert_eq!((-0.6666666666666666, 96.33333333333333), valley_trend);
/// ```
pub fn valley_trend(prices: &[f64], period: usize) -> (f64, f64) {
    let valleys = valleys(prices, period, 1);
    get_trend_line(&valleys)
}

/// Calculates the slope and intercept of the trend line fitted to all prices.
///
/// # Arguments
///
/// * `prices` - Slice of prices
///
/// # Examples
///
/// ```rust
/// let prices = vec![100.0, 102.0, 103.0, 101.0, 100.0];
/// let overall_trend = rust_ti::chart_trends::overall_trend(&prices);
/// assert_eq!((-0.1, 101.4), overall_trend);
/// ```
pub fn overall_trend(prices: &[f64]) -> (f64, f64) {
    let indexed_prices: Vec<(f64, usize)> =
        prices.iter().enumerate().map(|(i, &y)| (y, i)).collect();
    get_trend_line(&indexed_prices)
}

/// Calculates price trends and their slopes and intercepts.
///
/// # Arguments
///
/// * `prices` - Slice of prices
/// * `max_outliers` - Allowed consecutive trend-breaks before splitting
/// * `soft_r_squared_minimum` - Soft minimum value for r squared
/// * `soft_r_squared_maximum` - Soft maximum value for r squared
/// * `hard_r_squared_minimum` - Hard minimum value for r squared
/// * `hard_r_squared_maximum` - Hard maximum value for r squared
/// * `soft_standard_error_multiplier` - Soft standard error multiplier
/// * `hard_standard_error_multiplier` - Hard multiplier for the standard error
/// * `soft_reduced_chi_squared_multiplier` - Soft chi squared multiplier
/// * `hard_reduced_chi_squared_multiplier` - Hard chi squared multiplier
///
/// # Panics
///
/// Panics if `prices.is_empty()`
///
///
/// # Examples
///
/// ```rust
/// let prices = vec![100.0, 102.0, 103.0, 101.0, 99.0, 99.0, 102.0, 103.0, 106.0, 107.0, 105.0,
/// 104.0, 101.0, 97.0, 100.0];
/// let max_outliers = 1;
/// let soft_r_squared_minimum = 0.75;
/// let soft_r_squared_maximum = 1.0;
/// let hard_r_squared_minimum = 0.5;
/// let hard_r_squared_maximum = 1.5;
/// let soft_standard_error_multiplier = 2.0;
/// let hard_standard_error_multiplier = 3.0;
/// let soft_reduced_chi_squared_multiplier = 2.0;
/// let hard_reduced_chi_squared_multiplier = 3.0;
/// let trend_break_down = rust_ti::chart_trends::break_down_trends(
///     &prices,
///     max_outliers,
///     soft_r_squared_minimum,
///     soft_r_squared_maximum,
///     hard_r_squared_minimum,
///     hard_r_squared_maximum,
///     soft_standard_error_multiplier,
///     hard_standard_error_multiplier,
///     soft_reduced_chi_squared_multiplier,
///     hard_reduced_chi_squared_multiplier
/// );
/// assert_eq!(
///     vec![
///         (0, 2, 1.5, 100.16666666666667),
///         (2, 4, -2.0, 107.0),
///         (4, 9, 1.7714285714285714, 91.15238095238095),
///         (9, 11, -1.5, 120.33333333333333),
///         (11, 13, -3.5, 142.66666666666669)
///     ], trend_break_down);
/// ```
pub fn break_down_trends(
    prices: &[f64],
    max_outliers: usize,
    soft_r_squared_minimum: f64,
    soft_r_squared_maximum: f64,
    hard_r_squared_minimum: f64,
    hard_r_squared_maximum: f64,
    soft_standard_error_multiplier: f64,
    hard_standard_error_multiplier: f64,
    soft_reduced_chi_squared_multiplier: f64,
    hard_reduced_chi_squared_multiplier: f64,
) -> Vec<(usize, usize, f64, f64)> {
    if prices.is_empty() {
        panic!("Prices cannot be empty");
    };

    let mut outliers: Vec<usize> = Vec::new();
    let mut trends: Vec<(usize, usize, f64, f64)> = Vec::new();
    let mut current_slope = 0.0;
    let mut current_intercept = 0.0;
    let mut start_index: usize = 0;
    let mut end_index: usize = 1;
    let mut indexed_points: Vec<(f64, usize)> = Vec::new();
    let mut previous_standard_error = 10000.0;
    let mut previous_reduced_chi_squared = 10000.0;

    for (index, &price) in prices.iter().enumerate() {
        indexed_points.push((price, index));

        if index == 0 {
            continue;
        }
        if index > end_index {
            let current_trend = get_trend_line(&indexed_points);
            let (standard_error, r_squared, reduced_chi_squared) =
                goodness_of_fit(&indexed_points, &current_trend);

            let soft_break = standard_error
                > soft_standard_error_multiplier * previous_standard_error
                && (r_squared < soft_r_squared_minimum || r_squared > soft_r_squared_maximum)
                && reduced_chi_squared
                    > soft_reduced_chi_squared_multiplier * previous_reduced_chi_squared;

            let hard_break = r_squared < hard_r_squared_minimum
                || r_squared > hard_r_squared_maximum
                || standard_error > hard_standard_error_multiplier * previous_standard_error
                || reduced_chi_squared
                    > hard_reduced_chi_squared_multiplier * previous_reduced_chi_squared;

            if soft_break || hard_break {
                if outliers.len() < max_outliers {
                    outliers.push(index);
                    indexed_points.pop();
                    continue;
                };
                trends.push((start_index, end_index, current_slope, current_intercept));
                start_index = end_index;
                end_index = index;
                indexed_points = (start_index..=index).map(|x| (prices[x], x)).collect();
                let current_trend = get_trend_line(&indexed_points);
                current_slope = current_trend.0;
                current_intercept = current_trend.1;
                // if list bigger than 2
                if indexed_points.len() > 2 {
                    (previous_standard_error, _, previous_reduced_chi_squared) =
                        goodness_of_fit(&indexed_points, &current_trend);
                } else {
                    previous_standard_error = 10000.0;
                    previous_reduced_chi_squared = 10000.0;
                };
                outliers.clear();
            } else {
                previous_standard_error = standard_error;
                previous_reduced_chi_squared = reduced_chi_squared;
                current_slope = current_trend.0;
                current_intercept = current_trend.1;
                outliers.clear();
            }
        }
        end_index = index;
    }
    trends.push((start_index, end_index, current_slope, current_intercept));
    trends
}

fn goodness_of_fit(indexed_points: &[(f64, usize)], trend: &(f64, f64)) -> (f64, f64, f64) {
    let trend_line: Vec<f64> = indexed_points
        .iter()
        .map(|&(_, x)| trend.1 + trend.0 * x as f64)
        .collect();
    let observed_prices: Vec<f64> = indexed_points.iter().map(|&(y, _)| y).collect();

    let trend_length = trend_line.len();
    if trend_length != observed_prices.len() {
        panic!(
            "trend line length ({}) and prices length ({}) must be equal!",
            trend_length,
            trend_line.len()
        );
    };

    let observed_mean = mean(&observed_prices);
    let (sum_sq_residuals, total_squares): (f64, f64) =
        (0..trend_length).fold((0.0, 0.0), |(ssr, tss), i| {
            let resid = observed_prices[i] - trend_line[i];
            let total = observed_prices[i] - observed_mean;
            (ssr + resid.powi(2), tss + total.powi(2))
        });

    let standard_error = (sum_sq_residuals / trend_length as f64).sqrt();
    let r_squared = 1.0 - (sum_sq_residuals / total_squares);
    let reduced_chi_squared = sum_sq_residuals / trend_length as f64;
    (standard_error, r_squared, reduced_chi_squared)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn peaks_single_peak() {
        let highs = vec![101.26, 102.57, 102.32, 100.69];
        assert_eq!(vec![(102.57, 1)], peaks(&highs, 4_usize, 1usize));
    }

    #[test]
    fn peaks_multiple_peaks() {
        let highs = vec![101.26, 102.57, 102.32, 100.69, 100.83, 101.73, 102.01];
        assert_eq!(
            vec![(102.57, 1), (102.01, 6)],
            peaks(&highs, 4_usize, 1usize)
        );
    }

    #[test]
    fn peaks_multiple_peaks_same_period() {
        let highs = vec![101.26, 102.57, 102.57, 100.69, 100.83, 101.73, 102.01];
        assert_eq!(
            vec![(102.57, 2), (102.01, 6)],
            peaks(&highs, 4_usize, 1usize)
        );
    }

    #[test]
    #[should_panic]
    fn peaks_panic() {
        let highs = vec![101.26, 102.57, 102.57, 100.69, 100.83, 101.73, 102.01];
        peaks(&highs, 40_usize, 1usize);
    }

    #[test]
    fn valleys_single_valley() {
        let lows = vec![100.08, 98.75, 100.14, 98.98, 99.07, 100.1, 99.96];
        assert_eq!(vec![(98.75, 1)], valleys(&lows, 7_usize, 1usize));
    }

    #[test]
    fn valleys_multiple_valleys() {
        let lows = vec![100.08, 98.75, 100.14, 98.98, 99.07, 100.1, 99.96];
        assert_eq!(
            vec![(98.75, 1), (98.98, 3)],
            valleys(&lows, 4_usize, 1usize)
        );
    }

    #[test]
    fn valleys_multiple_valleys_same_period() {
        let lows = vec![98.75, 98.75, 100.14, 98.98, 99.07, 100.1, 99.96];
        assert_eq!(
            vec![(98.75, 1), (98.98, 3)],
            valleys(&lows, 4_usize, 1usize)
        );
    }

    #[test]
    #[should_panic]
    fn valleys_panic() {
        let lows = vec![98.75, 98.75, 100.14, 98.98, 99.07, 100.1, 99.96];
        valleys(&lows, 40_usize, 1usize);
    }

    #[test]
    fn peaks_trend() {
        let highs = vec![101.26, 102.57, 102.32, 100.69, 100.83, 101.73, 102.01];
        assert_eq!(
            (-0.11199999999999762, 102.68199999999999),
            peak_trend(&highs, 4_usize)
        );
    }

    #[test]
    fn valleys_trend() {
        let lows = vec![100.08, 98.75, 100.14, 98.98, 99.07, 100.1, 99.96];
        assert_eq!((0.11500000000000199, 98.635), valley_trend(&lows, 4_usize));
    }

    #[test]
    fn overall_trends() {
        let prices = vec![100.2, 100.46, 100.53, 100.38, 100.19];
        assert_eq!((-0.010000000000000852, 100.372), overall_trend(&prices));
    }

    #[test]
    fn break_down_trends_std_dev() {
        let prices = vec![100.2, 100.46, 100.53, 100.38, 100.19];
        let max_outliers = 1;
        let soft_r_squared_minimum = 0.75;
        let soft_r_squared_maximum = 1.0;
        let hard_r_squared_minimum = 0.5;
        let hard_r_squared_maximum = 1.5;
        let soft_standard_error_multiplier = 2.0;
        let hard_standard_error_multiplier = 3.0;
        let soft_reduced_chi_squared_multiplier = 2.0;
        let hard_reduced_chi_squared_multiplier = 3.0;
        let trend_break_down = break_down_trends(
            &prices,
            max_outliers,
            soft_r_squared_minimum,
            soft_r_squared_maximum,
            hard_r_squared_minimum,
            hard_r_squared_maximum,
            soft_standard_error_multiplier,
            hard_standard_error_multiplier,
            soft_reduced_chi_squared_multiplier,
            hard_reduced_chi_squared_multiplier,
        );
        assert_eq!(
            vec![
                (0, 2, 0.16499999999999915, 100.23166666666665),
                (2, 4, -0.1700000000000017, 100.87666666666668)
            ],
            trend_break_down
        );
    }
}
