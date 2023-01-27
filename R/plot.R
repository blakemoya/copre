#' Create a CopRe Result Plot
#'
#' @param x A `copre_result` object.
#' @param ... Additional arguments discarded from `plot`.
#' @param func Either 'distribution', 'density', or 'gradient'.
#' @param confint A decimal value indicating the confidence interval width (e.g.
#'   0.95 for a 95% confidence interval). Defaults to `NULL`, in which case no
#'   confidence intervals will be drawn.
#' @param use_ggplot A logical value indicating whether to use `ggplot2` instead
#'   of the base `plot` function.
#'
#' @return None.
#' @export
plot.copre_result <- function(x, ..., func = 'density', confint = NULL,
                              use_ggplot = TRUE) {
  plot.grideval_result(grideval(x, func = func), confint = confint,
                       use_ggplot = use_ggplot)
}

#' Create a SeqRe Result Plot
#'
#' @param x A `seqre_result` object.
#' @param ... Additional arguments discarded from `plot`.
#' @param func Either 'distribution', 'density', or 'gradient'.
#' @param confint A decimal value indicating the confidence interval width (e.g.
#'   0.95 for a 95% confidence interval). Defaults to `NULL`, in which case no
#'   confidence intervals will be drawn.
#' @param use_ggplot A logical value indicating whether to use `ggplot2` instead
#'   of the base `plot` function.
#'
#' @return None.
#' @export
plot.seqre_result <- function(x, ..., func = 'density', confint = NULL,
                              use_ggplot = TRUE) {
  plot.grideval_result(grideval(x, func = func), confint = confint,
                       use_ggplot = use_ggplot)
}

#' Create a Plot of a `grideval_result` Object
#'
#' @param x A `grideval_result` object.
#' @param ... Additional arguments discarded from `plot`.
#' @param confint A decimal value indicating the confidence interval width (e.g.
#'   0.95 for a 95 percent confidence interval). Defaults to `NULL`, in which
#'   case no confidence intervals will be drawn.
#' @param use_ggplot A logical value indicating whether to use `ggplot2` instead
#'   of the base `plot` function.
#'
#' @return None.
#' @export
plot.grideval_result <- function(x, ..., confint = NULL, use_ggplot = TRUE) {
  if (requireNamespace('ggplot2', quietly = TRUE) & use_ggplot) {
    p <- autoplot.grideval_result(x, ..., confint = confint)
    print(p)
  } else {
    grd <- x$grid
    df <- data.frame(Value = rep(grd, each = nrow(x)),
                     K = rep(1:nrow(x), length(grd)),
                     X = as.numeric(x))
    K <- Value <- X <- NULL # Avoid unbound global note
    with(df, {
      plot(Value, X, type = 'n', ylab = paste0(toupper(substring(x$func, 1, 1)),
                                               substring(x$func, 2)))
      for (k in 1:max(K)) {
        k_idx = which(K == k)
        lines(Value[k_idx], X[k_idx], col = rgb(0.75, 0.75, 0.75, 0.25))
      }
    })
    lines(grd, apply(x, 2, mean), col = 'red')
    data <- sort(x$args$data)
    n <- length(data)
    if(x$func == 'density') {
      points(data + stats::runif(n, -0.001, 0.001),
             stats::runif(n, -max(df$X) / 50, 0),
             pch = 16, cex = 0.5, col = rgb(0, 0, 0, 0.5))
    } else if (x$func == 'distribution') {
      e_cdf <- stats::ecdf(data)
      lines(grd, e_cdf(grd))
      if (!is.null(confint)) {
        err_int <- 1 - confint
        eps_dkw <- sqrt(log(2 / err_int) / (2 * n))
        upper_dkw <- stats::stepfun(data, c(0, pmin(e_cdf(data) + eps_dkw, 1)))
        lower_dkw <- stats::stepfun(data, c(0, pmax(e_cdf(data) - eps_dkw, 0)))
        eps_clt <- stats::qnorm(1 - (err_int / 2)) *
          sqrt(e_cdf(data) * (1 - e_cdf(data)) / n)
        upper_clt <- stats::stepfun(data, c(0, e_cdf(data) + eps_clt))
        lower_clt <- stats::stepfun(data, c(0, e_cdf(data) - eps_clt))
        lines(grd, upper_clt(grd), col = rgb(0, 0, 0, 0.5))
        lines(grd, lower_clt(grd), col = rgb(0, 0, 0, 0.5))
        lines(grd, lower_dkw(grd), lty = 'longdash',
              col = rgb(0.5, 0.5, 0.5, 0.5))
        lines(grd, upper_dkw(grd), lty = 'longdash',
              col = rgb(0.5, 0.5, 0.5, 0.5))
      }
    }
  }
}
