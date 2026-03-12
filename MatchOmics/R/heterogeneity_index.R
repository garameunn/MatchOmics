#' Prepare heterogeneity indices
#'
#' @param omics_type omics type
#' @param data feature matrix (samples x features)
#' @param index optional user-supplied heterogeneity matrix
#' @return data.frame of heterogeneity covariates
#' @export
prepare_heterogeneity <- function(omics_type, data = NULL, index = NULL) {
  if (!is.null(index)) {
    return(as.data.frame(index))
  }

  omics_type <- match.arg(
    omics_type,
    c("microbiome", "proteome", "metabolome", "transcriptome", "custom")
  )

  if (is.null(data)) {
    stop("Either `data` or `index` must be provided.")
  }

  data <- as.matrix(data)

  if (omics_type == "microbiome") {
    rel <- sweep(data, 1, rowSums(data, na.rm = TRUE), "/")
    shannon <- -rowSums(rel * log(rel), na.rm = TRUE)
    richness <- rowSums(data > 0, na.rm = TRUE)
    return(data.frame(shannon = shannon, richness = richness))
  }

  if (omics_type %in% c("proteome", "metabolome", "transcriptome")) {
    pcs <- stats::prcomp(data, center = TRUE, scale. = TRUE)$x[, 1:2, drop = FALSE]
    return(as.data.frame(pcs))
  }

  stop("For omics_type = 'custom', please provide `index` directly.")
}
