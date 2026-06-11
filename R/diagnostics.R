#' Kish effective sample size
#'
#' @param w Numeric vector of mean-normalised weights (sum = N).
#' @return Scalar Neff.
#' @export
compute_neff <- function(w) {
  w <- w[!is.na(w) & w > 0]
  if (length(w) == 0) return(NA_real_)
  (sum(w))^2 / sum(w^2)
}

#' Standardised mean difference per covariate
#'
#' @param df       Data frame containing the covariates and a binary grouping
#'   column.
#' @param covnames Character vector of covariate column names.
#' @param group_col Name of the binary grouping column (default
#'   \code{"marker_class"}).
#' @param treated_level Value that identifies the treated group (default
#'   \code{"up"}).
#' @return Named numeric vector of SMD values.
#' @export
compute_smd <- function(df, covnames,
                        group_col     = "marker_class",
                        treated_level = "up") {
  sapply(covnames, function(v) {
    x  <- df[[v]]
    g  <- df[[group_col]] == treated_level
    m1 <- mean(x[g],  na.rm = TRUE)
    m0 <- mean(x[!g], na.rm = TRUE)
    s1 <- stats::var(x[g],  na.rm = TRUE)
    s0 <- stats::var(x[!g], na.rm = TRUE)
    pooled_sd <- sqrt((s1 + s0) / 2)
    if (pooled_sd == 0) return(NA_real_)
    (m1 - m0) / pooled_sd
  })
}

#' Summarise pre/post matching SMD
#'
#' @param pre_smd  Named numeric vector of pre-matching SMD values.
#' @param post_smd Named numeric vector of post-matching SMD values.
#' @return Data frame with columns \code{covariate}, \code{smd_pre},
#'   \code{smd_post}.
#' @export
smd_summary <- function(pre_smd, post_smd) {
  data.frame(
    covariate = names(pre_smd),
    smd_pre   = unname(pre_smd),
    smd_post  = unname(post_smd),
    row.names = NULL
  )
}

#' Love plot: pre vs post matching SMD
#'
#' Requires \pkg{ggplot2}.
#'
#' @param smd_df    Output of \code{\link{smd_summary}}.
#' @param title     Plot title.
#' @param threshold Dashed reference line (default \code{0.1}).
#' @return A \code{ggplot} object.
#' @export
love_plot <- function(smd_df, title = "Love Plot", threshold = 0.1) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for love_plot(). ",
         "Install with: install.packages('ggplot2')")

  df_long <- rbind(
    data.frame(covariate = smd_df$covariate,
               smd       = abs(smd_df$smd_pre),
               timing    = "Before matching"),
    data.frame(covariate = smd_df$covariate,
               smd       = abs(smd_df$smd_post),
               timing    = "After matching")
  )
  df_long$timing <- factor(df_long$timing,
                            levels = c("Before matching", "After matching"))

  ggplot2::ggplot(df_long,
                  ggplot2::aes(x = .data$smd, y = .data$covariate,
                               colour = .data$timing,
                               shape  = .data$timing)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_vline(xintercept = threshold,
                        linetype = "dashed", colour = "grey40") +
    ggplot2::scale_colour_manual(
      values = c("Before matching" = "#E74C3C",
                 "After matching"  = "#2ECC71")) +
    ggplot2::labs(title  = title,
                  x      = "Absolute Standardised Mean Difference",
                  y      = NULL,
                  colour = NULL, shape = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}
