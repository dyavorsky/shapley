# ## shapley: Shapley/LMG variable importance for regression models
#
# Computes Shapley values by averaging each predictor's marginal R² contribution
# over all possible orderings (LMG method). Supports lm, glm(binomial), and polr.
# Control variables (z_vars) are included in every subset model, so values
# attribute only the incremental fit above the controls.
# Computational cost is O(2^p) model fits; warns if p > 15.


shapley <- function(y_var, x_vars, z_vars = character(0), data,
                    y_type = c("continuous", "binary", "ordinal")) {

  y_type <- match.arg(y_type)

  all_vars      <- c(y_var, x_vars, z_vars)
  complete_data <- data[complete.cases(data[, all_vars, drop=FALSE]), ]
  n             <- nrow(complete_data)
  n_missing     <- nrow(data) - n

  if (n < length(x_vars) + length(z_vars) + 1)
    stop("Insufficient complete cases for Shapley analysis")

  if (y_type == "continuous") {
    non_numeric <- Filter(function(v) !is.numeric(complete_data[[v]]), c(x_vars, z_vars))
    if (length(non_numeric) > 0)
      stop("y_type = 'continuous' requires numeric predictors; convert these first: ",
           paste(non_numeric, collapse = ", "))
  }

  get_r2_fn  <- .make_r2_fn(y_var, complete_data, y_type)
  model_type <- switch(y_type,
    continuous = "lm",
    binary     = "glm(binomial)",
    ordinal    = "polr"
  )

  values  <- .shapley_engine(x_vars, z_vars, get_r2_fn)
  r2_full <- get_r2_fn(c(z_vars, x_vars))

  structure(
    list(
      values     = values,
      r2         = r2_full,
      model_type = model_type,
      n          = n,
      n_missing  = n_missing
    ),
    class = "shapley_result"
  )
}


print.shapley_result <- function(x, digits = 3, ...) {
  cat("Shapley/LMG Importance\n")
  cat(strrep("-", 40), "\n")
  cat("Model: ", x$model_type, ", n = ", x$n, sep="")
  if (x$n_missing > 0) cat(", ", x$n_missing, " missing", sep="")
  cat(", R² = ", round(x$r2, digits), "\n\n", sep="")
  print(round(x$values, digits))
  invisible(x)
}


# Build the R²-computing closure for a given y_type.
# The closure captures complete_data and y_var; called with a vector of predictor names.
.make_r2_fn <- function(y_var, complete_data, y_type) {

  if (y_type == "continuous") {
    # Pre-compute correlation matrix once; each subset's R² = r_yS' solve(R_SS) r_yS.
    # Replaces 2^p lm() calls with 2^p small matrix solves — ~1000x faster at p=15.
    num_vars <- names(complete_data)[vapply(complete_data, is.numeric, logical(1L))]
    cor_mat  <- cor(complete_data[, num_vars, drop = FALSE])
    function(preds) {
      if (length(preds) == 0) return(0)
      idx  <- c(y_var, preds)
      cm   <- cor_mat[idx, idx, drop = FALSE]
      r_yS <- cm[1L, -1L, drop = FALSE]
      R_SS <- cm[-1L, -1L, drop = FALSE]
      as.numeric(r_yS %*% solve(R_SS) %*% t(r_yS))
    }

  } else if (y_type == "binary") {
    null_dev <- glm(as.formula(paste(y_var, "~ 1")), family=binomial,
                    data=complete_data)$deviance
    function(preds) {
      if (length(preds) == 0) return(0)
      m <- glm(as.formula(paste(y_var, "~", paste(preds, collapse=" + "))),
               family=binomial, data=complete_data)
      1 - m$deviance / null_dev
    }

  } else if (y_type == "ordinal") {
    null_dev <- polr(as.formula(paste(y_var, "~ 1")), data=complete_data,
                     Hess=TRUE)$deviance
    function(preds) {
      if (length(preds) == 0) return(0)
      m <- polr(as.formula(paste(y_var, "~", paste(preds, collapse=" + "))),
                data=complete_data, Hess=TRUE)
      1 - m$deviance / null_dev
    }
  }
}


# Core Shapley/LMG engine. Enumerates all 2^p subsets of x_vars, caches R² for
# each, then applies the Shapley weighting formula. z_vars are passed into
# get_r2_fn and held constant across all subsets.
.shapley_engine <- function(x_vars, z_vars, get_r2_fn) {
  p <- length(x_vars)

  if (p > 15)
    warning("Shapley computation with ", p, " x variables requires up to ", 2^p,
            " model fits and may be slow.")

  # r2_cache[k+1] = R² for the subset of x_vars encoded by binary integer k.
  # Bit j (0-indexed) set means x_vars[j+1] is included; z_vars always included.
  r2_cache <- numeric(2^p)
  for (k in 0:(2^p - 1)) {
    bits            <- as.logical(intToBits(k)[seq_len(p)])
    r2_cache[k + 1] <- get_r2_fn(c(z_vars, x_vars[bits]))
  }

  # For each variable i, sum weighted marginal contributions R²(S∪{i}) - R²(S)
  # over all subsets S of x_vars not containing i.
  # Weight = |S|! * (p-|S|-1)! / p!  =  fraction of orderings where i follows exactly S.
  shapley_vals <- setNames(numeric(p), x_vars)

  for (i in seq_len(p)) {
    sv        <- 0
    other_idx <- setdiff(seq_len(p), i)
    q         <- p - 1

    for (k in 0:(2^q - 1)) {
      S_idx  <- other_idx[as.logical(intToBits(k)[seq_len(q)])]
      weight <- factorial(length(S_idx)) * factorial(p - length(S_idx) - 1) / factorial(p)
      key_S  <- sum(2^(S_idx - 1))
      key_Si <- sum(2^(c(S_idx, i) - 1))
      sv     <- sv + weight * (r2_cache[key_Si + 1] - r2_cache[key_S + 1])
    }

    shapley_vals[[i]] <- sv
  }

  shapley_vals
}
