# Contributing to RustTI

Thank you for considering contributing—every improvement is appreciated, big or small!

---

## 🙌 How to Contribute

1. **Fork and clone the repository**  
2. **Make your changes**  
3. **Open a Pull Request (PR)** on [GitHub](https://github.com/0100101001010000/RustTI/pulls)
4. **Tag @0100101001010000** in your PR for a review

See [open issues](https://github.com/0100101001010000/RustTI/issues) if you want to start with something small.

---

## 🛠️ What to Work On?

- Remove `test_` prefix from all test functions
- Refactor parabolic T/P system: remove single functions; let bulk function determine trend (like `volatility` system)
- Allow indicators to use non-normal distributions (for variance, standard deviation, etc.)

### New Indicator Ideas

- Indicator showing price distribution over a period (count unique prices, order ascending)
- Median/mode absolute deviation based on deviation from mode/median (not mean)
- McGinley dynamic versions of indicators using MA (Donchian, Keltner, Supertrend, PPO, RVI, etc.)

---

## ➕ Adding a New Indicator

1. **Open an Issue**  
   - Describe the indicator, what it’s used for, and how it’s calculated (with source/reference)
2. **Implement the indicator**  
   - Add documentation and unit tests
3. **Verify results**  
   - Add a tab to `assets/rust_ti_hand_calcs.ods` with hand calculations to ensure test accuracy

---

## 🧪 Code Style & Testing

- Format code with `rustfmt` or `cargo fmt`
- Run tests with `cargo test` before submitting your PR

---

Thanks again for your interest and contributions!

