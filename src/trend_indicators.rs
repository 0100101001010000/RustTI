//! Trend Indicators
//!
//! Trend indicators show the trend direction of an asset

/// `single` module holds functions that return a singular values
pub mod single {
    use crate::basic_indicators::single::{max, min};
    /// The `aroon_up` indicator tracks the uptrends in the `aroon_indicator` and is used to
    /// calculate the `aroon_oscillator`.
    ///
    /// The Aroon up is included in the return of the `aroon_indicator`. If the caller wants the
    /// Aroon up, Aroon down, and Aroon oscillator, then it is easier to call the `aroon_indicator`
    /// rather that the three seperately.
    ///
    /// Standard period to use is 25 but the caller determines how many prices they want to provide
    /// and the period will be determined from the length of `highs`.
    ///
    /// # Arguments
    ///
    /// * `highs` - Slice of highs
    ///
    /// # Examples
    ///
    /// ```
    /// let highs = vec![103.0, 102.0, 107.0, 104.0, 100.0];
    /// let aroon_up = rust_ti::trend_indicators::single::aroon_up(&highs);
    /// assert_eq!(50.0, aroon_up);
    /// ```
    pub fn aroon_up(highs: &[f64]) -> f64 {
        if highs.is_empty() {
            panic!("Highs cannot be empty")
        };

        let period = highs.len() - 1; // current period should be excluded from length
        let period_max = max(highs);
        let periods_since_max = period - highs.iter().rposition(|&x| x == period_max).unwrap();
        return 100.0 * ((period as f64 - periods_since_max as f64) / period as f64);
    }

    /// The `aroon_down` indicator tracks the downtrends in the `aroon_indicator` and is used to
    /// calculate the `aroon_oscillator`.
    ///
    /// The Aroon down is included in the return of the `aroon_indicator`. If the caller wants the
    /// Aroon up, Aroon down, and Aroon oscillator, then it is easier to call the `aroon_indicator`
    /// rather that the three seperately.
    ///
    /// Standard period to use is 25 but the caller determines how many prices they want to provide
    /// and the period will be determined from the length of `lows`.
    ///
    /// # Arguments
    ///
    /// * `low` - Slice of lows
    ///
    /// # Examples
    ///
    /// ```
    /// let lows = vec![98.0, 95.0, 101.0, 100.0, 97.0];
    /// let aroon_down = rust_ti::trend_indicators::single::aroon_down(&lows);
    /// assert_eq!(25.0, aroon_down);
    /// ```
    pub fn aroon_down(lows: &[f64]) -> f64 {
        if lows.is_empty() {
            panic!("Lows cannot be empty")
        };

        let period = lows.len() - 1; // current period should be excluded from length
        let period_min = min(lows);
        let periods_since_min = period - lows.iter().rposition(|&x| x == period_min).unwrap();
        return 100.0 * ((period as f64 - periods_since_min as f64) / period as f64);
    }

    /// The `aroon_oscillators` takes the difference between the Aroon up and the Aroon down to
    /// give a general sense of trend.
    ///
    /// The Aroon oscillator is returned in the `aroon_indicator` so it is easier to call that one
    /// function rather than 3.
    ///
    /// # Arguments
    ///
    /// * `aroon_up` - Aroon up for the period
    /// * `aroon_down` - Aroon down for the period
    ///
    /// # Examples
    ///
    /// ```
    /// let aroon_up = 50.0;
    /// let aroon_down = 25.0;
    /// let aroon_oscillator = rust_ti::trend_indicators::single::aroon_oscillator(&aroon_up,
    /// &aroon_down);
    /// assert_eq!(25.0, aroon_oscillator);
    /// ```
    pub fn aroon_oscillator(aroon_up: &f64, aroon_down: &f64) -> f64 {
        return aroon_up - aroon_down;
    }

    /// The `aroon_indicator` returns the Aroon up, Aroon down, and Aroon oscillator in that order
    ///
    /// # Arguments
    ///
    /// * `high` - Slice of highs
    /// * `low` - Slice of lows
    ///
    /// # Examples
    ///
    /// ```
    /// let highs = vec![103.0, 102.0, 107.0, 104.0, 100.0];
    /// let lows = vec![98.0, 95.0, 101.0, 100.0, 97.0];
    /// let aroon_indicator = rust_ti::trend_indicators::single::aroon_indicator(&highs, &lows);
    /// assert_eq!((50.0, 25.0, 25.0), aroon_indicator);
    /// ```
    pub fn aroon_indicator(highs: &[f64], lows: &[f64]) -> (f64, f64, f64) {
        if highs.len() != lows.len() {
            panic!(
                "Length of highs ({}) must match length of lows ({})",
                highs.len(),
                lows.len()
            )
        };

        let aroon_up = aroon_up(&highs);
        let aroon_down = aroon_down(&lows);
        let aroon_oscillaor = aroon_oscillator(&aroon_up, &aroon_down);
        return (aroon_up, aroon_down, aroon_oscillaor);
    }

    /// The `long_parabolic_time_price_system` returns Stop and Reverse (SaR) points based on a significant point
    /// (period low) or previous SaR, time, and price. The function is to be used when considering
    /// a long position.
    ///
    /// Welles uses the significant point as the SaR when no SaR points have been calculated
    /// afterwhich he uses the previous SaR. He also uses an acceleration factor that increases by 0.02 every
    /// day a new high is hit and never goes over 0.2, but in this function the value is determined
    /// by the function caller.
    ///
    /// Welles has a rule that the SaR at t+1 cannot be above the low for t and t-1.
    ///
    /// # Arguments
    ///
    /// * `previous_sar` - Previous Stop and Reverse, if none use the significant point, which is
    /// the period low.
    /// * `extreme_point` - Highest high for the period being observed
    /// * `acceleration_factor` - Factor used to multiply the difference between extreme point and
    /// previous SaR.
    /// * `low` - Lowest low for t or t-1
    ///
    /// # Examples
    ///
    /// ```
    /// let previous_sar = 50.09306;
    /// let extreme_point = 52.35;
    /// let acceleration_factor = 0.02;
    /// let low = 50.6;
    /// let parabolic_time_price_system =
    /// rust_ti::trend_indicators::single::long_parabolic_time_price_system(&previous_sar,
    /// &extreme_point, &acceleration_factor, &low);
    /// assert_eq!(50.1381988, parabolic_time_price_system);
    ///
    /// let previous_sar = 51.96;
    /// let extreme_point = 54.2;
    /// let acceleration_factor = 0.12;
    /// let low = 52.1;
    /// let parabolic_time_price_system =
    /// rust_ti::trend_indicators::single::long_parabolic_time_price_system(&previous_sar,
    /// &extreme_point, &acceleration_factor, &low);
    /// assert_eq!(52.1, parabolic_time_price_system);
    /// ```
    pub fn long_parabolic_time_price_system(
        previous_sar: &f64,
        extreme_point: &f64,
        acceleration_factor: &f64,
        low: &f64,
    ) -> f64 {
        let sar = previous_sar + acceleration_factor * (extreme_point - previous_sar);
        if &sar > low {
            return *low;
        };
        return sar;
    }

    /// The `short_parabolic_time_price_system` returns Stop and Reverse (SaR) points based on a significant point
    /// (period low) or previous SaR, time, and price. The function is to be used when considering
    /// a short position.
    ///
    /// Welles uses the significant point as the SaR when no SaR points have been calculated
    /// afterwhich he uses the previous SaR. He also uses an acceleration factor that increases by 0.02 every
    /// day a new low is hit and never goes over 0.2, but in this function the value is determined
    /// by the function caller.
    ///
    /// Welles has a rule that the SaR at t+1 cannot be above the high for t and t-1.
    ///
    /// # Arguments
    ///
    /// * `previous_sar` - Previous Stop and Reverse, if none use the significant point, which is
    ///  the period high.
    /// * `extreme_point` - Lowest low for the period being observed
    /// * `acceleration_factor` - Factor used to multiply the difference between extreme point and
    /// previous SaR.
    /// * `high` - Highest high for t or t-1
    ///
    /// # Examples
    ///
    /// ```
    /// let previous_sar = 58.0;
    /// let extreme_point = 56.3;
    /// let acceleration_factor = 0.02;
    /// let high = 50.6;
    /// let parabolic_time_price_system =
    /// rust_ti::trend_indicators::single::short_parabolic_time_price_system(&previous_sar,
    /// &extreme_point, &acceleration_factor, &high);
    /// assert_eq!(57.966, parabolic_time_price_system);
    ///
    /// let previous_sar = 57.7816384;
    /// let extreme_point = 55.5;
    /// let acceleration_factor = 0.08;
    /// let low = 58.1;
    /// let parabolic_time_price_system =
    /// rust_ti::trend_indicators::single::short_parabolic_time_price_system(&previous_sar,
    /// &extreme_point, &acceleration_factor, &low);
    /// assert_eq!(58.1, parabolic_time_price_system);
    /// ```
    pub fn short_parabolic_time_price_system(
        previous_sar: &f64,
        extreme_point: &f64,
        acceleration_factor: &f64,
        high: &f64,
    ) -> f64 {
        let sar = previous_sar - acceleration_factor * (previous_sar - extreme_point);
        if &sar < high {
            return *high;
        };
        return sar;
    }
}

/// `bulk` module holds functions that return multiple vaues
pub mod bulk {
    use crate::basic_indicators::single::{max, min};
    use crate::trend_indicators::single;
    /// The `aroon_up` indicator tracks the uptrends in the `aroon_indicator` and is used to
    /// calculate the `aroon_oscillator`.
    ///
    /// The Aroon up is included in the return of the `aroon_indicator`. If the caller wants the
    /// Aroon up, Aroon down, and Aroon oscillator, then it is easier to call the `aroon_indicator`
    /// rather that the three seperately.
    ///
    /// # Arguments
    ///
    /// * `highs` - Slice of highs
    /// * `period` - Period over which to calculate the Aroon up
    ///
    /// # Examples
    ///
    /// ```
    /// let highs = vec![103.0, 102.0, 107.0, 104.0, 100.0, 102.0, 99.0];
    /// let period: usize = 5;
    /// let aroon_up = rust_ti::trend_indicators::bulk::aroon_up(&highs, &period);
    /// assert_eq!(vec![50.0, 25.0, 0.0], aroon_up);
    /// ```
    pub fn aroon_up(highs: &[f64], period: &usize) -> Vec<f64> {
        let length = highs.len();
        if &length < period {
            panic!(
                "Period ({}) cannot be longer than length of highs ({})",
                period, length
            )
        };

        let mut aroon_ups = Vec::new();
        let loop_max = length - period + 1;
        for i in 0..loop_max {
            aroon_ups.push(single::aroon_up(&highs[i..i + period]));
        }
        return aroon_ups;
    }

    /// The `aroon_down` indicator tracks the downtrends in the `aroon_indicator` and is used to
    /// calculate the `aroon_oscillator`.
    ///
    /// The Aroon down is included in the return of the `aroon_indicator`. If the caller wants the
    /// Aroon up, Aroon down, and Aroon oscillator, then it is easier to call the `aroon_indicator`
    /// rather that the three seperately.
    ///
    /// # Arguments
    ///
    /// * `low` - Slice of lows
    /// * `period` - Period over which to calculate the Aroon down
    ///
    /// # Examples
    ///
    /// ```
    /// let lows = vec![98.0, 95.0, 101.0, 100.0, 97.0, 98.0, 97.0];
    /// let period: usize = 5;
    /// let aroon_down = rust_ti::trend_indicators::bulk::aroon_down(&lows, &period);
    /// assert_eq!(vec![25.0, 0.0, 100.0], aroon_down);
    /// ```
    pub fn aroon_down(lows: &[f64], period: &usize) -> Vec<f64> {
        let length = lows.len();
        if &length < period {
            panic!(
                "Period ({}) cannot be longer than length of lows ({})",
                period, length
            )
        };

        let mut aroon_downs = Vec::new();
        let loop_max = length - period + 1;
        for i in 0..loop_max {
            aroon_downs.push(single::aroon_down(&lows[i..i + period]));
        }
        return aroon_downs;
    }

    /// The `aroon_oscillators` takes the difference between the Aroon up and the Aroon down to
    /// give a general sense of trend.
    ///
    /// # Arguments
    ///
    /// * `aroon_up` - Slice of Aroon ups
    /// * `aroon_down` - Slice Aroon downs
    ///
    /// # Examples
    ///
    /// ```
    /// let aroon_up = vec![50.0, 25.0, 0.0];
    /// let aroon_down = vec![25.0, 0.0, 100.0];
    /// let aroon_oscillator = rust_ti::trend_indicators::bulk::aroon_oscillator(&aroon_up,
    /// &aroon_down);
    /// assert_eq!(vec![25.0, 25.0, -100.0], aroon_oscillator);
    /// ```
    pub fn aroon_oscillator(aroon_up: &[f64], aroon_down: &[f64]) -> Vec<f64> {
        let length = aroon_up.len();
        if length != aroon_down.len() {
            panic!(
                "Length of Aroon up ({}) and Aroon down ({}) must match",
                length,
                aroon_down.len()
            )
        };

        let mut aroon_oscillators = Vec::new();
        for i in 0..length {
            aroon_oscillators.push(single::aroon_oscillator(&aroon_up[i], &aroon_down[i]));
        }
        return aroon_oscillators;
    }

    /// The `aroon_indicator` returns the Aroon up, Aroon down, and Aroon oscillator in that order
    ///
    /// # Arguments
    ///
    /// * `high` - Slice of highs
    /// * `low` - Slice of lows
    /// * `period` - Period over which to calculate the Aroon indicator
    ///
    /// # Examples
    ///
    /// ```
    /// let highs = vec![103.0, 102.0, 107.0, 104.0, 100.0, 102.0, 99.0];
    /// let lows = vec![98.0, 95.0, 101.0, 100.0, 97.0, 98.0, 97.0];
    /// let period: usize = 5;
    /// let aroon_indicator = rust_ti::trend_indicators::bulk::aroon_indicator(&highs, &lows,
    /// &period);
    /// assert_eq!(vec![(50.0, 25.0, 25.0), (25.0, 0.0, 25.0), (0.0, 100.0, -100.0)], aroon_indicator);
    /// ```
    pub fn aroon_indicator(highs: &[f64], lows: &[f64], period: &usize) -> Vec<(f64, f64, f64)> {
        let length = highs.len();
        if length != lows.len() {
            panic!(
                "Length of highs ({}) must match length of lows ({})",
                highs.len(),
                lows.len()
            )
        };
        if &length < period {
            panic!(
                "Period ({}) cannot be longer than lengths of highs and lows ({})",
                period, length
            )
        };

        let mut aroon_indicators = Vec::new();
        let loop_max = length - period + 1;
        for i in 0..loop_max {
            aroon_indicators.push(single::aroon_indicator(
                &highs[i..i + period],
                &lows[i..i + period],
            ));
        }
        return aroon_indicators;
    }

    /// The `parabolic_time_price_system` returns Stop and Reverse (SaR) points given a slices of
    ///  high and low prices.
    ///
    /// Welles uses an acceleration factor that starts at 0.02, increases by 0.02 every
    /// day a new high/low is hit and never goes over 0.2, but in these values are determined
    /// by the function caller.
    ///
    /// Welles always makes the assumption that a short trade existed before starting to calculate
    /// the parabolic time price system. To reflect this the `parabolic_time_price_system` starts
    /// by calculating SaR from a short position.
    ///
    /// # Arguments
    ///
    /// * `highs` - Slice of highs.
    /// * `lows` - Slice of lows.
    /// * `acceleration_factor_start` - The value to start at for the acceleration factor. Standard
    /// is 0.02.
    /// * `acceleration_factor_max` - The maximum value the acceleration factor can be. Standard is
    /// 0.2.
    /// * `acceleration_factor_step` - The step by which to increase the acceleration factor every
    /// time it hits a new high/low. Standard is 0.02.
    ///
    /// # Examples
    ///
    /// ```
    /// let highs = vec![52.35, 52.1, 51.8, 52.1, 52.5, 52.8, 52.5, 53.5, 53.5, 53.8, 54.2, 53.4, 53.5, 54.4, 55.2, 55.7, 57.0, 57.5, 58.0, 57.7, 58.0, 57.5, 57.0, 56.7, 57.5, 56.7, 56.0, 56.2, 54.8, 55.5, 54.7, 54.0, 52.5, 51.0, 51.5, 51.7, 53.0];
    /// let lows = vec![51.5, 51.0, 50.5, 51.25, 51.7, 51.85, 51.5, 52.5, 52.5, 53.0, 52.5, 52.5, 52.1, 53.0, 54.0, 55.0, 56.0, 56.5, 57.0, 56.5, 57.3, 56.7, 56.3, 56.2, 56.0, 55.5, 55.0, 54.9, 54.0, 54.5, 53.8, 53.0, 51.5, 50.0, 50.5, 50.2, 51.5]; 
    /// let acceleration_factor_start = 0.02;
    /// let acceleration_factor_max = 0.2;
    /// let acceleration_factor_step = 0.02;
    /// let parabolic_time_price_system = rust_ti::trend_indicators::bulk::parabolic_time_price_system(&highs, &lows,
    /// &acceleration_factor_start, &acceleration_factor_max, &acceleration_factor_step);
    /// assert_eq!(
    ///     vec![50.0, 50.047, 50.09306, 50.1381988, 50.182434824, 50.27513743104, 50.4266291851776, 50.569031434067, 50.8035089193416, 51.0192282057943, 51.2973053852148, 51.6456287389891, 51.9521532903104, 52.1, 52.1, 52.596, 53.15472, 53.923776, 54.6390208, 55.31121664, 55.848973312, 56.2791786496, 57.966, 57.89536, 57.7816384, 57.599107328, 57.3391965952, 57.046493003776, 56.6199839832474, 56.2531862255927, 55.8606764294979, 55.3457546721883, 54.5766037377506, 53.6612829902005, 52.9290263921604, 52.3432211137283],
    ///     parabolic_time_price_system);
    /// ```
    pub fn parabolic_time_price_system(
        highs: &[f64],
        lows: &[f64],
        acceleration_factor_start: &f64,
        acceleration_factor_max: &f64,
        acceleration_factor_step: &f64,
    ) -> Vec<f64> {
        // TODO: Position is being produced here, if it were passed to function could have one functio
        // instead of two
        if highs.is_empty() || lows.is_empty() {
            panic!("Highs or lows cannot be empty")
        };
        let length = highs.len();
        if length != lows.len() {
            panic!(
                "Highs ({}) and lows ({}) must be the same length",
                length,
                lows.len()
            )
        };
        let mut acceleration_factor = *acceleration_factor_start;
        // TODO: Reverse this
        //let mut sars = vec![single::short_parabolic_time_price_system(
        let mut sars = vec![single::long_parabolic_time_price_system(
            &50.0,
            //&highs[0],
            &lows[0],
            &acceleration_factor,
            &highs[0],
        )];
        let mut position = 'l';
        let mut position_start = 0;
        for i in 1..length {
            let previous_sar = sars[i - 1];
            println!("SAR ({}), highs ({}), lows ({}), i ({})", previous_sar, highs[i], lows[i], i);
            if position == 's' && highs[i] > previous_sar {
                position = 'l';
                let period_max = max(&highs[..i]);
                let previous_min = min(&lows[i - 2..i - 1]);
                acceleration_factor = *acceleration_factor_start;
                let pivoted_sar = max(&highs[position_start..i-1]);
                position_start = i;
                sars.push(single::long_parabolic_time_price_system(
                    &pivoted_sar,
                    &period_max,
                    &acceleration_factor,
                    &previous_min,
                ));
            } else if position == 's' {
                let period_min = min(&lows[..i]);
                if period_min < lows[i] && &acceleration_factor < acceleration_factor_max {
                    acceleration_factor = acceleration_factor + acceleration_factor_step;
                };
                let previous_max = max(&highs[i - 2..i - 1]);
                sars.push(single::short_parabolic_time_price_system(
                    &previous_sar,
                    &period_min,
                    &acceleration_factor,
                    &previous_max,
                ));
            } else if position == 'l' && lows[i] < previous_sar {
                position = 's';
                let period_min = min(&lows[..i]);
                acceleration_factor = *acceleration_factor_start;
                let previous_max = max(&highs[i - 2..i - 1]);
                let pivoted_sar = min(&lows[position_start..i-1]);
                position_start = i;
                sars.push(single::short_parabolic_time_price_system(
                    &pivoted_sar,
                    &period_min,
                    &acceleration_factor,
                    &previous_max,
                ));
            } else if position == 'l' {
                let period_max = max(&highs[..i]);
                if period_max < highs[i] && &acceleration_factor < acceleration_factor_max {
                    acceleration_factor = acceleration_factor + acceleration_factor_step;
                };
                // TODO: when index 1 there is no day before yesterday...
                let previous_min = min(&lows[i - 2..i - 1]);
                println!("max: ({}), acc fac: ({}), prev min ({})", period_max, acceleration_factor, previous_min);
                sars.push(single::long_parabolic_time_price_system(
                    &previous_sar,
                    &period_max,
                    &acceleration_factor,
                    &previous_min,
                ));
            } else {
                panic!("Position ({}) is not a supported position!", position)
            }
        }
        return sars;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_aroon_up() {
        let highs = vec![101.26, 102.57, 102.32, 100.69];
        assert_eq!(33.33333333333333, single::aroon_up(&highs));
    }

    #[test]
    #[should_panic]
    fn singe_aroon_up_panic() {
        let highs = Vec::new();
        single::aroon_up(&highs);
    }

    #[test]
    fn bulk_aroon_up() {
        let highs = vec![101.26, 102.57, 102.32, 100.69, 100.83, 101.73, 102.01];
        assert_eq!(
            vec![33.33333333333333, 0.0, 0.0, 100.0],
            bulk::aroon_up(&highs, &4)
        );
    }

    #[test]
    #[should_panic]
    fn bulk_aroon_up_panic() {
        let highs = vec![101.26, 102.57, 102.32, 100.69, 100.83, 101.73, 102.01];
        bulk::aroon_up(&highs, &40);
    }

    #[test]
    fn single_aroon_down() {
        let lows = vec![100.08, 98.75, 100.14, 98.98];
        assert_eq!(33.33333333333333, single::aroon_down(&lows));
    }

    #[test]
    #[should_panic]
    fn single_aroon_down_panic() {
        let lows = Vec::new();
        single::aroon_down(&lows);
    }

    #[test]
    fn bulk_aroon_down() {
        let lows = vec![100.08, 98.75, 100.14, 98.98, 99.07, 100.1, 99.96];
        assert_eq!(
            vec![33.33333333333333, 0.0, 33.33333333333333, 0.0],
            bulk::aroon_down(&lows, &4)
        );
    }

    #[test]
    #[should_panic]
    fn bulk_aroon_down_panic() {
        let lows = vec![100.08, 98.75, 100.14, 98.98, 99.07, 100.1, 99.96];
        bulk::aroon_down(&lows, &40);
    }

    #[test]
    fn single_aroon_oscillator() {
        assert_eq!(
            0.0,
            single::aroon_oscillator(&33.33333333333333, &33.33333333333333)
        );
    }

    #[test]
    fn bulk_aroon_oscillator() {
        let aroon_up = vec![33.33333333333333, 0.0, 0.0, 100.0];
        let aroon_down = vec![33.33333333333333, 0.0, 33.33333333333333, 0.0];
        assert_eq!(
            vec![0.0, 0.0, -33.33333333333333, 100.0],
            bulk::aroon_oscillator(&aroon_up, &aroon_down)
        );
    }

    #[test]
    #[should_panic]
    fn bulk_aroon_oscillator_up_panic() {
        let aroon_up = vec![33.33333333333333, 0.0, 0.0];
        let aroon_down = vec![33.33333333333333, 0.0, 33.33333333333333, 0.0];
        bulk::aroon_oscillator(&aroon_up, &aroon_down);
    }

    #[test]
    #[should_panic]
    fn bulk_aroon_oscillator_down_panic() {
        let aroon_up = vec![33.33333333333333, 0.0, 0.0, 100.0];
        let aroon_down = vec![33.33333333333333, 0.0, 33.33333333333333];
        bulk::aroon_oscillator(&aroon_up, &aroon_down);
    }

    #[test]
    fn single_aroon_indicator() {
        let lows = vec![100.08, 98.75, 100.14, 98.98];
        let highs = vec![101.26, 102.57, 102.32, 100.69];
        assert_eq!(
            (33.33333333333333, 33.33333333333333, 0.0),
            single::aroon_indicator(&highs, &lows)
        );
    }

    #[test]
    #[should_panic]
    fn single_aroon_indicator_high_panic() {
        let lows = vec![100.08, 98.75, 100.14, 98.98];
        let highs = vec![101.26, 102.57, 102.32];
        single::aroon_indicator(&highs, &lows);
    }

    #[test]
    #[should_panic]
    fn single_aroon_indicator_low_panic() {
        let lows = vec![100.08, 98.75, 100.14];
        let highs = vec![101.26, 102.57, 102.32, 100.69];
        single::aroon_indicator(&highs, &lows);
    }

    #[test]
    fn bulk_aroon_indicator() {
        let highs = vec![101.26, 102.57, 102.32, 100.69, 100.83, 101.73, 102.01];
        let lows = vec![100.08, 98.75, 100.14, 98.98, 99.07, 100.1, 99.96];
        assert_eq!(
            vec![
                (33.33333333333333, 33.33333333333333, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 33.33333333333333, -33.33333333333333),
                (100.0, 0.0, 100.0)
            ],
            bulk::aroon_indicator(&highs, &lows, &4)
        );
    }

    #[test]
    #[should_panic]
    fn bulk_aroon_indicator_high_panic() {
        let highs = vec![102.57, 102.32, 100.69, 100.83, 101.73, 102.01];
        let lows = vec![100.08, 98.75, 100.14, 98.98, 99.07, 100.1, 99.96];
        bulk::aroon_indicator(&highs, &lows, &4);
    }

    #[test]
    #[should_panic]
    fn bulk_aroon_indicator_low_panic() {
        let highs = vec![101.26, 102.57, 102.32, 100.69, 100.83, 101.73, 102.01];
        let lows = vec![98.75, 100.14, 98.98, 99.07, 100.1, 99.96];
        bulk::aroon_indicator(&highs, &lows, &4);
    }

    #[test]
    #[should_panic]
    fn bulk_aroon_indicator_period_panic() {
        let highs = vec![101.26, 102.57, 102.32, 100.69, 100.83, 101.73, 102.01];
        let lows = vec![100.08, 98.75, 100.14, 98.98, 99.07, 100.1, 99.96];
        bulk::aroon_indicator(&highs, &lows, &40);
    }
}