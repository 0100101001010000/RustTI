//! # Moving Averages
//!
//! The `moving_average` module has functions used to calculate the moving average
//!
//! The moving average has three different types available (simple, smoothed, exponential) that the
//! caller can choose from.
//!
//! ## Bulk
//!
//! * [`mcginley_dynamic`](bulk::mcginley_dynamic) - Calculates the McGinley Dynamic
//! * [`moving_average`](bulk::moving_average) - Calculates different types Moving Averages
//!
//! ## Single
//!
//! * [`mcginley_dynamic`](single::mcginley_dynamic) - Calculates the McGinley Dynamic
//! * [`moving_average`](single::moving_average) - Calculates different types Moving Averages

/// `single` module holds functions that return a singular values
pub mod single {
    use crate::basic_indicators::single::mean;
    use crate::MovingAverageType;

    /// Calculates the Moving Average
    ///
    /// # Arguments
    ///
    /// * `prices` - Slice of prices
    /// * `moving_average_type` - Variant of [`MovingAverageType`]
    ///
    /// # Panics
    ///
    /// Panic if `prices.is_empty()`
    ///
    /// # Examples
    ///
    /// ```rust
    /// let prices = vec![100.0, 102.0, 103.0, 101.0, 100.0];
    ///
    /// let simple_moving_average =
    ///     rust_ti::moving_average::single::moving_average(
    ///         &prices,
    ///         rust_ti::MovingAverageType::Simple
    ///     );
    /// assert_eq!(101.2, simple_moving_average);
    ///
    /// let exponential_moving_average =
    ///     rust_ti::moving_average::single::moving_average(
    ///         &prices,
    ///         rust_ti::MovingAverageType::Exponential
    ///     );
    /// assert_eq!(100.99526066350714, exponential_moving_average);
    ///
    /// let smoothed_moving_average =
    ///     rust_ti::moving_average::single::moving_average(
    ///         &prices,
    ///         rust_ti::MovingAverageType::Smoothed
    ///     );
    /// assert_eq!(101.11375535459305, smoothed_moving_average);
    /// ```
    #[inline]
    pub fn moving_average(prices: &[f64], moving_average_type: MovingAverageType) -> f64 {
        if prices.is_empty() {
            panic!("Prices is empty")
        };
        match moving_average_type {
            MovingAverageType::Simple => mean(prices),
            MovingAverageType::Smoothed => personalised_moving_average(prices, 1.0, 0.0),
            MovingAverageType::Exponential => personalised_moving_average(prices, 2.0, 1.0),
            MovingAverageType::Personalised {
                alpha_num,
                alpha_den,
            } => personalised_moving_average(prices, alpha_num, alpha_den),
            _ => panic!("Unsupported MovingiAverageType"),
        }
    }

    /// Internal: Generic moving average with custom alpha parameters.
    ///
    /// # Panics
    /// Panics if prices is empty or denominator would be zero.
    #[inline]
    fn personalised_moving_average(
        prices: &[f64],
        alpha_nominator: f64,
        alpha_denominator: f64,
    ) -> f64 {
        let length = prices.len() as f64;
        if length == 1.0 {
            return prices[0];
        };

        if length + alpha_denominator == 0.0 {
            panic!(
                "The length of prices ({}) and the alpha_denominator ({}) add up to 0",
                length, alpha_denominator
            );
        };

        let alpha: f64 = alpha_nominator / (length + alpha_denominator);
        let multiplicator = 1.0 - alpha;
        let mut price_sum: f64 = 0.0;
        let mut denominator_sum: f64 = 0.0;
        for (index, price) in prices.iter().rev().enumerate() {
            let multiplactor_powd = multiplicator.powi(index as i32);
            denominator_sum += multiplactor_powd;
            price_sum += price * multiplactor_powd;
        }
        price_sum / denominator_sum
    }

    /// Calculates the McGinley dynamic
    ///
    /// # Arguments
    ///
    /// * `latest_price` - Most recent price
    /// * `previous_mcginley_dynamic` - Previous McGinley dynamic (if none 0.0)
    /// * `period` - Length of the observed period
    ///
    /// # Panics
    ///
    /// Panics if `period` <= 0
    ///
    /// # Examples
    ///
    /// ```rust
    /// let prices = vec![100.0, 102.0, 103.0, 101.0, 100.0];
    /// let period:usize = 5;
    /// let mcginley_dynamic =
    ///     rust_ti::moving_average::single::mcginley_dynamic(100.0, 0.0, period);
    /// assert_eq!(100.0, mcginley_dynamic);
    ///
    /// let next_prices = vec![102.0, 103.0, 101.0, 100.0, 99.0];
    /// let next_mcginley_dynamic =
    ///     rust_ti::moving_average::single::mcginley_dynamic(99.0, mcginley_dynamic, period);
    /// assert_eq!(99.79179592886295, next_mcginley_dynamic);
    /// ```
    #[inline]
    pub fn mcginley_dynamic(
        latest_price: f64,
        previous_mcginley_dynamic: f64,
        period: usize,
    ) -> f64 {
        if period <= 0 {
            panic!("Cannot have a 0 period");
        };
        if previous_mcginley_dynamic == 0.0 {
            return latest_price;
        };
        let base = latest_price / previous_mcginley_dynamic;
        previous_mcginley_dynamic
            + ((latest_price - previous_mcginley_dynamic) / (period as f64 * base.powi(4)))
    }
}

/// `bulk` module holds functions that return multiple valus for `moving_average`
pub mod bulk {
    use crate::moving_average::single;
    use crate::MovingAverageType;

    /// Calculates the moving average
    ///
    /// # Arguments
    ///
    /// * `prices` - Slice of prices
    /// * `moving_average_type` - Variant of [`MovingAverageType`]
    /// * `period` - Period over which to calculate the moving average
    ///
    /// # Panics
    ///
    /// Panics if `period` > `prices.len()`
    ///
    /// # Examples
    ///
    /// ```rust
    /// let prices = vec![100.0, 102.0, 103.0, 101.0, 100.0];
    /// let period: usize = 3;
    ///
    /// let simple_moving_average =
    ///     rust_ti::moving_average::bulk::moving_average(
    ///         &prices,
    ///         rust_ti::MovingAverageType::Simple,
    ///         period
    ///     );
    /// assert_eq!(
    ///     vec![101.66666666666667, 102.0, 101.33333333333333],
    ///     simple_moving_average
    /// );
    ///
    /// let exponential_moving_average =
    ///     rust_ti::moving_average::bulk::moving_average(
    ///         &prices,
    ///         rust_ti::MovingAverageType::Exponential,
    ///         period
    ///     );
    /// assert_eq!(
    ///     vec![102.28571428571429, 101.71428571428571, 100.71428571428571],
    ///     exponential_moving_average
    /// );
    ///
    /// let smoothed_moving_average =
    ///     rust_ti::moving_average::bulk::moving_average(
    ///         &prices,
    ///         rust_ti::MovingAverageType::Smoothed,
    ///         period
    ///     );
    /// assert_eq!(
    ///     vec![102.05263157894737, 101.8421052631579, 100.94736842105264],
    ///     smoothed_moving_average
    /// );
    /// ```
    #[inline]
    pub fn moving_average(
        prices: &[f64],
        moving_average_type: MovingAverageType,
        period: usize,
    ) -> Vec<f64> {
        let length = prices.len();
        if period > length {
            panic!(
                "Period ({}) cannot be longer than the length of provided prices ({})",
                period, length
            );
        };

        let mut moving_averages = Vec::with_capacity(length - period + 1);
        for window in prices.windows(period) {
            moving_averages.push(single::moving_average(window, moving_average_type));
        }
        moving_averages
    }

    /// Calculates the McGinley dynamic
    ///
    /// # Arguments
    ///
    /// * `prices` - Slice of prices
    /// * `previous_mcginley_dynamic` - Previous McGinley dynamic (if none 0.0)
    /// * `period` - Period over which to calculate the McGinley dynamic
    ///
    /// # Panics
    ///
    /// Panics if `period` > `prices.len()`
    ///
    /// # Examples
    ///
    /// ```rust
    /// let prices = vec![100.0, 102.0, 103.0, 101.0, 100.0];
    /// let period: usize = 3;
    ///
    /// let mcginley_dynamic =
    ///     rust_ti::moving_average::bulk::mcginley_dynamic(
    ///         &prices,
    ///         0.0,
    ///         period
    ///     );
    /// assert_eq!(vec![103.0, 102.2789387706985, 101.44764169058672], mcginley_dynamic);
    ///
    /// let mcginley_dynamic =
    ///     rust_ti::moving_average::bulk::mcginley_dynamic(
    ///         &prices[1..],
    ///         103.0,
    ///         period
    ///     );
    /// assert_eq!(vec![102.2789387706985, 101.44764169058672], mcginley_dynamic);
    /// ```
    #[inline]
    pub fn mcginley_dynamic(
        prices: &[f64],
        previous_mcginley_dynamic: f64,
        period: usize,
    ) -> Vec<f64> {
        let length = prices.len();
        if period > length {
            panic!(
                "Period ({}) cannot be longer than the length of provided prices ({})",
                period, length
            );
        };

        let mut mcginley_dynamics = Vec::with_capacity(length - period + 1);
        let mut mcginley_dynamic =
            single::mcginley_dynamic(prices[period - 1], previous_mcginley_dynamic, period);
        mcginley_dynamics.push(mcginley_dynamic);
        for &price in prices.iter().skip(period) {
            mcginley_dynamic = single::mcginley_dynamic(price, mcginley_dynamic, period);
            mcginley_dynamics.push(mcginley_dynamic);
        }
        mcginley_dynamics
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_simple_moving_average() {
        let prices = vec![100.2, 100.46, 100.53, 100.38, 100.19];
        let simple_ma = single::moving_average(&prices, crate::MovingAverageType::Simple);
        assert_eq!(100.352, simple_ma);
    }

    #[test]
    fn single_exponential_moving_average() {
        let prices = vec![100.2, 100.46, 100.53, 100.38, 100.19];
        let exponential_ma = single::moving_average(&prices, crate::MovingAverageType::Exponential);
        assert_eq!(100.32810426540287, exponential_ma);
    }

    #[test]
    fn single_smoothed_moving_average() {
        let prices = vec![100.2, 100.46, 100.53, 100.38, 100.19];
        let smoothed_ma = single::moving_average(&prices, crate::MovingAverageType::Smoothed);
        assert_eq!(100.34228938600666, smoothed_ma);
    }

    #[test]
    fn single_personalised_moving_average() {
        let prices = vec![100.2, 100.46, 100.53, 100.38, 100.19];
        let personalised_ma = single::moving_average(
            &prices,
            crate::MovingAverageType::Personalised {
                alpha_num: 5.0,
                alpha_den: 3.0,
            },
        );
        assert_eq!(100.27405995388162, personalised_ma)
    }

    #[test]
    #[should_panic]
    fn single_moving_average_panic() {
        let prices = Vec::new();
        single::moving_average(&prices, crate::MovingAverageType::Simple);
    }

    #[test]
    #[should_panic]
    fn single_moving_average_personalised_ma_panic() {
        let prices = vec![100.2, 100.46, 100.53, 100.38, 100.19];
        single::moving_average(
            &prices,
            crate::MovingAverageType::Personalised {
                alpha_num: 5.0,
                alpha_den: -5.0,
            },
        );
    }

    #[test]
    fn bulk_simlpe_moving_average() {
        let prices = vec![100.2, 100.46, 100.53, 100.38, 100.19];
        let period: usize = 3;
        let simple_ma = bulk::moving_average(&prices, crate::MovingAverageType::Simple, period);
        assert_eq!(
            vec![100.39666666666666, 100.456666666666666, 100.36666666666667],
            simple_ma
        );
    }

    #[test]
    fn bulk_exponential_moving_average() {
        let prices = vec![100.2, 100.46, 100.53, 100.38, 100.19];
        let period: usize = 3;
        let exponential_ma =
            bulk::moving_average(&prices, crate::MovingAverageType::Exponential, period);
        assert_eq!(
            vec![100.46285714285715, 100.4342857142857, 100.29285714285713],
            exponential_ma
        );
    }

    #[test]
    fn bulk_smoothed_moving_average() {
        let prices = vec![100.2, 100.46, 100.53, 100.38, 100.19];
        let period: usize = 3;
        let smoothed_ma = bulk::moving_average(&prices, crate::MovingAverageType::Smoothed, period);
        assert_eq!(
            vec![100.43842105263158, 100.4442105263158, 100.32157894736842],
            smoothed_ma
        );
    }

    #[test]
    fn bulk_personalised_moving_average() {
        let prices = vec![100.2, 100.46, 100.53, 100.38, 100.19];
        let period: usize = 3;
        let personalised_ma = bulk::moving_average(
            &prices,
            crate::MovingAverageType::Personalised {
                alpha_num: 5.0,
                alpha_den: 3.0,
            },
            period,
        );
        assert_eq!(
            vec![100.5125581395349, 100.40279069767443, 100.22441860465118],
            personalised_ma
        );
    }

    #[test]
    #[should_panic]
    fn bulk_moving_average_panic() {
        let prices = vec![100.2, 100.46, 100.53, 100.38, 100.19];
        let period: usize = 30;
        bulk::moving_average(&prices, crate::MovingAverageType::Simple, period);
    }

    #[test]
    fn single_mcginley_dynamic_no_previous() {
        let mcginley_dynamic = single::mcginley_dynamic(100.19, 0.0_f64, 5_usize);
        assert_eq!(100.19, mcginley_dynamic);
    }

    #[test]
    fn single_mcginley_dynamic_previous() {
        let mcginley_dynamic = single::mcginley_dynamic(100.21, 100.19, 5_usize);
        assert_eq!(100.19399680766176, mcginley_dynamic);
    }

    #[test]
    #[should_panic]
    fn single_mcginley_dynamic_panic() {
        single::mcginley_dynamic(100.0, 0.0_f64, 0_usize);
    }

    #[test]
    fn bulk_mcginley_dynamic() {
        let prices = vec![100.2, 100.46, 100.53, 100.38, 100.19];
        let period: usize = 3;
        assert_eq!(
            vec![100.53, 100.47970046511769, 100.38201189376744],
            bulk::mcginley_dynamic(&prices, 0.0_f64, period)
        );
        assert_eq!(
            vec![100.47970046511769, 100.38201189376744],
            bulk::mcginley_dynamic(&prices[1..], 100.53, period)
        );
    }

    #[test]
    #[should_panic]
    fn bulk_mcginley_dynamic_panic() {
        let prices = vec![100.2, 100.46, 100.53, 100.38, 100.19];
        let period: usize = 30;
        bulk::mcginley_dynamic(&prices, 0.0_f64, period);
    }
}
