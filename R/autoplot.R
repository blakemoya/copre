#' Create a CopRe Result ggplot
#'
#' @param x A `copre_result` object.
#' @param ... Additional arguments discarded from `plot`.
#' @param func Either 'distribution', 'density', or 'gradient'.
#' @param confint A decimal value indicating the confidence interval width (e.g.
#'   0.95 for a 95% confidence interval). Defaults to `NULL`, in which case no
#'   confidence intervals will be drawn.
#'
#' @return A `ggplot` object.
autoplot.copre_result <- function(x, ..., func = 'density', confint = NULL) {
  autoplot.grideval_result(grideval(x, func = func), confint = confint)
}

#' Create a SeqRe Result ggplot
#'
#' @param x A `seqre_result` object.
#' @param ... Additional arguments discarded from `plot`.
#' @param func Either 'distribution', 'density', or 'gradient'.
#' @param confint A decimal value indicating the confidence interval width (e.g.
#'   0.95 for a 95% confidence interval). Defaults to `NULL`, in which case no
#'   confidence intervals will be drawn.
#'
#' @return A `ggplot` object.
autoplot.seqre_result <- function(x, ..., func = 'density', confint = NULL) {
  autoplot.grideval_result(grideval(x, func = func), confint = confint)
}

#' Create a ggplot of a `grideval_result` Object
#'
#' @param x A `grideval_result` object.
#' @param ... Additional arguments discarded from `plot`.
#' @param confint A decimal value indicating the confidence interval width (e.g.
#'   0.95 for a 95 percent confidence interval). Defaults to `NULL`, in which
#'   case no confidence intervals will be drawn.
#'
#' @return A `ggplot` object.
autoplot.grideval_result <- function(x, ..., confint = NULL) {
  grd <- x$grid
  df <- data.frame(Value = rep(grd, each = nrow(x)),
                   K = rep(1:nrow(x), length(grd)),
                   X = as.numeric(x))
  K <- Value <- X <- NULL # Avoid unbound global note
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Value, y = X, group = K)) +
    ggplot2::ylab(paste0(toupper(substring(x$func, 1, 1)),
                         substring(x$func, 2))) +
    ggplot2::geom_line(alpha = 0.25, color = 'grey') +
    ggplot2::geom_line(data = data.frame(Value = grd, X = apply(x, 2, mean)),
                       ggplot2::aes(group = 0),
                       color = 'red') +
    ggplot2::theme_bw()
  data <- sort(x$args$data)
  n <- length(data)
  if(x$func == 'density') {
    p <- p +
      ggplot2::geom_point(
        data = data.frame(Value = data + stats::runif(n, -0.001, 0.001),
                          X = stats::runif(n, -max(df$X) / 50, 0),
                          K = 0),
        shape = 16, size = 0.5, alpha = 0.5
        )
  } else if (x$func == 'distribution') {
    e_cdf <- stats::ecdf(data)
    p <- p + ggplot2::stat_function(fun = e_cdf, ggplot2::aes(group = 0),
                                    geom = 'step', n = 1001)
    if (!is.null(confint)) {
      err_int <- 1 - confint
      eps_dkw <- sqrt(log(2 / err_int) / (2 * n))
      upper_dkw <- stats::stepfun(data, c(0, pmin(e_cdf(data) + eps_dkw, 1)))
      lower_dkw <- stats::stepfun(data, c(0, pmax(e_cdf(data) - eps_dkw, 0)))
      eps_clt <- stats::qnorm(1 - (err_int / 2)) *
        sqrt(e_cdf(data) * (1 - e_cdf(data)) / n)
      upper_clt <- stats::stepfun(data, c(0, e_cdf(data) + eps_clt))
      lower_clt <- stats::stepfun(data, c(0, e_cdf(data) - eps_clt))
      p <- p +
        ggplot2::stat_function(fun = upper_clt, ggplot2::aes(group = 0),
                               geom = 'step', n = 1001,
                               size = 0.25, alpha = 0.5) +
        ggplot2::stat_function(fun = lower_clt, ggplot2::aes(group = 0),
                               geom = 'step', n = 1001,
                               size = 0.25, alpha = 0.5) +
        ggplot2::stat_function(fun = lower_dkw, ggplot2::aes(group = 0),
                               geom = 'step', n = 1001,
                               size = 0.25, color = 'grey50',
                               linetype = 'longdash', alpha = 0.5) +
        ggplot2::stat_function(fun = upper_dkw, ggplot2::aes(group = 0),
                               geom = 'step', n = 1001,
                               size = 0.25, color = 'grey50',
                               linetype = 'longdash', alpha = 0.5)
    }
  }
  return(p)
}

#' Register `autoplot` methods to `ggplot2`
#'
#' @return None
register_autoplot_s3_methods = function() {
  if (requireNamespace('ggplot2', quietly = TRUE)) {
    register_s3_method('ggplot2', 'autoplot', 'copre_result')
    register_s3_method('ggplot2', 'autoplot', 'seqre_result')
    register_s3_method('ggplot2', 'autoplot', 'grideval_result')
  }
}
