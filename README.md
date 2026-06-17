# shapley

An R package for Shapley/LMG variable importance in regression models.

`shapley()` computes each predictor's average marginal R² contribution over all possible orderings — the same LMG method used in [`relaimpo`](https://cran.r-project.org/package=relaimpo) for linear regression — and extends it to **binary outcomes** (logistic regression, McFadden pseudo-R²) and **ordinal outcomes** (proportional odds regression, McFadden pseudo-R²). For continuous outcomes, R² for each predictor subset is computed directly from the correlation matrix rather than by refitting a model, making it fast even at 20+ predictors.

## Installation

```r
# install.packages("pak")
pak::pak("dyavorsky/shapley")
```

## Usage

```r
library(shapley)
data(swiss)

x_vars <- c("Agriculture", "Examination", "Education", "Catholic", "Infant.Mortality")

# Continuous outcome (linear regression)
shapley("Fertility", x_vars, data = swiss, y_type = "continuous")

# Binary outcome (logistic regression)
swiss$high_fert <- as.integer(swiss$Fertility > median(swiss$Fertility))
shapley("high_fert", x_vars, data = swiss, y_type = "binary")

# Ordinal outcome (proportional odds regression)
swiss$fert_cat <- cut(swiss$Fertility,
                      breaks = quantile(swiss$Fertility, c(0, 1/3, 2/3, 1)),
                      include.lowest = TRUE, labels = c("low", "mid", "high"))
swiss$fert_cat <- factor(swiss$fert_cat, ordered = TRUE)
shapley("fert_cat", x_vars, data = swiss, y_type = "ordinal")
```

Control variables (`z_vars`) can be specified and are held fixed across all subset models, so Shapley values reflect only the incremental fit above those controls.

See `vignette("getting-started", package = "shapley")` for a comparison against `relaimpo` and side-by-side results across all three outcome types on the `swiss` dataset.
