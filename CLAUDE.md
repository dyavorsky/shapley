# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

`shapley` is an R package (R >= 4.1) that computes Shapley/LMG variable importance values for regression models. It supports continuous, binary, and ordinal outcomes. Dependencies: `MASS`, `stats`.

## Development Commands

```r
# Load package in development
devtools::load_all()

# Run R CMD check
devtools::check()

# Install locally
devtools::install()

# Run tests (tests/testthat/ when they exist)
devtools::test()
testthat::test_file("tests/testthat/test-shapley.R")  # single test file
```

## Architecture

All source lives in a single file: `R/shapley.R` (133 lines).

**Call chain:**

```
shapley()           # public API — validates inputs, dispatches
  ├── .make_r2_fn() # returns a closure that computes R² for a given predictor subset
  │     branches on y_type: lm (continuous) / glm+binomial (binary) / polr (ordinal)
  │     binary & ordinal use McFadden pseudo-R²; continuous uses standard R²
  └── .shapley_engine()  # enumerates all 2^p predictor subsets, caches R² values,
                          # applies Shapley weighting: |S|! × (p−|S|−1)! / p!
```

**Return value:** `shapley_result` S3 object with fields `values`, `r2_full`, `y_type`, `n`. `print.shapley_result()` is the only S3 method.

**Key behavioral notes:**
- `z_vars` (control variables) are included in every subset model so Shapley values measure incremental fit above controls.
- Complete-case analysis: rows with any missing values across `y_var`, `x_vars`, `z_vars` are dropped before fitting.
- Warns if `length(x_vars) > 15` due to O(2^p) complexity.
