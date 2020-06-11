#' @import purrr
#' @import stats
#' @import furrr
#' @import future
#' @import vroom
#' @importFrom magrittr %>%
#' @importFrom utils capture.output
#' @aliases NULL
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))

#' Build Linear Regression Model
#'
#' Build the linear Regression model,
#' using the bag of little bootstrap method to find the coefficients
#'
#' @param formula linear regression model
#'
#' @param data data frame or vector of file names
#' @param m number of subsamples
#' @param B number of bootstraps
#'
#' @export
blblm <- function(formula, data, m = 10, B = 5000) {
  if (is.data.frame(data)){
    data_list <- split_data(data, m)
    n_tol <- nrow(data)
  } else{
    data_list <- data %>%
    map(vroom, col_types = cols())
    n_tol <- data_list %>%
      map_dbl(nrow) %>%
      reduce(`+`)
  }

  parallel <- readline(prompt = "Do you want to use parallelization to generate the model? [y/n]\n")
  if (parallel == "y"){
    n_cores <- readline(prompt = "How many cores you want to use for parallellization:") %>%
      as.integer()
    stopifnot(n_cores <= availableCores())
    suppressWarnings(plan(multiprocess, workers = n_cores))

    estimates <- future_map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = n_tol, B = B)
    )

  } else{
    estimates <- map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = n_tol, B = B))
  }

  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}


#' split data
#'
#' split data into m parts of approximated equal sizes using random sampling
#'
#' @param data input data
#' @param m number of subsamples
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' compute the estimates
#'
#' Use bootstrap on each subsample to compute its coefficients
#'
#' @param formula linear regression model
#' @param data input data
#' @param n size of each blb dataset
#' @param B number of bootstraps
lm_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}


#' compute the regression estimates for a blb dataset
#'
#' @param formula linear regression model
#' @param data input data
#' @param n number of rows in input data
lm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs)
}


#' estimate the regression estimates based on given the number of repetitions
#'
#' @param formula linear regression model
#' @param data input data
#' @param freqs number of times of each row of a subsample appears in a blb dataset
lm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wrong variable from the global scope.
  environment(formula) <- environment()
  fit <- lm(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


#' Obtain model's coefficient
#'
#' compute the coefficients from fit
#'
#' @param fit fitted model
blbcoef <- function(fit) {
  coef(fit)
}


#' Obtain model's sigma
#'
#' compute sigma from fit
#'
#' @param fit fitted model
blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' Print the linear model
#'
#' @param x the fitted model
#'
#' @param ... further arguments passed to or from other methods
#'
#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}


#' Print the sigma
#'
#' Compute the sigma and its confidence interval if confidence = TRUE
#'
#' @param object the fitted model
#'
#' @param confidence boolean
#' @param level level of confidence
#' @param ... further arguments passed to or from other methods
#'
#' @export
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' Get the coefficient of fitted model
#'
#' @param object the fitted model
#'
#' @param ... further arguments passed to or from other methods
#'
#' @export
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}


#' Confidence interval of parameters
#'
#' compute all parameters' CIs if parm is not specified,
#' otherwise compute specific parm's CI
#'
#' @param object fitted model
#'
#' @param parm parameters
#' @param level level of confidence
#' @param ... further arguments passed to or from other methods
#'
#' @export
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' Predict outcome with new data
#'
#' Compute response variable and its predicted confidence interval given a set of new data
#'
#' @param object fitted model
#'
#' @param new_data data matrix
#' @param confidence boolean
#' @param level level of confidence
#' @param ... further arguments passed to or from other methods
#'
#' @export
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms.formula(object$formula,data = new_data), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}

##################################################################################################
##################################################################################################

#' Build generalized linear model
#'
#' Set up glm using blb method, default family is gaussian
#'
#' @param formula linear model
#' @param data data.frame/string vector
#' @param family string type of glm
#' @param m integer
#' @param B integer
#'
#' @export
blbglm <- function(formula, data, family = gaussian(), m = 10, B = 5000){
  if (is.data.frame(data)){
    data_list <- split_data(data, m)
    n_tol <- nrow(data)

  } else{
    data_list <- data %>%
      map(vroom, col_types = cols())
      n_tol <- data_list %>%
      map_dbl(nrow) %>%
      reduce(`+`)
  }

  parallel <- readline(prompt = "Do you want to use parallelization to generate the model? [y/n]\n")
  if (parallel == "y"){
    n_cores <- readline(prompt = "How many cores you want to use for parallellization:") %>%
      as.integer()
    stopifnot(n_cores <= availableCores())
    suppressWarnings(plan(multiprocess, workers = n_cores))

    estimates <- future_map(
      data_list,
      ~ glm_each_subsample(formula = formula, data = ., family = family, n = n_tol, B = B))

  } else{
    estimates <- map(
      data_list,
      ~ glm_each_subsample(formula = formula, data = ., family = family, n = n_tol, B = B))
  }

  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}



#' Set up bootstraps for each sub-samples
#'
#' @param formula linear model
#' @param data data.frame
#' @param family string type of glm
#' @param n integer each boot's size
#' @param B integer number of bootstraps
glm_each_subsample <- function(formula, data, family, n, B) {
  replicate(B, glm_each_boot(formula, data, family, n), simplify = FALSE)
}


#' Build glm for each bootstrap dataset
#'
#' @param formula linear model
#' @param data data.frame
#' @param family string type of glm
#' @param n integer each boot's size
glm_each_boot <- function(formula, data, family, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  glm1(formula, data, family, freqs)
}


#' estimate the regression estimates based on given the number of repetitions
#'
#' @param formula linear model
#' @param data data.frame
#' @param family string type of glm
#' @param freqs weights of each row in sub-sample
glm1 <- function(formula, data, family, freqs){
  # drop the original closure of formula,
  # otherwise the formula will pick a wrong variable from the global scope.
  environment(formula) <- environment()
  fit <- glm(formula, data, family = family, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}
