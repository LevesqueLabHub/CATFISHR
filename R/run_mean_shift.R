#' Run mean shift clustering on RNA and CNV Mahalanobis distances for one sample
#'
#' Performs mean shift clustering on log2-transformed RNA and CNV Mahalanobis
#' distances for a single sample. Bandwidths are tested by mean silhouette
#' score, and the best bandwidth is selected while optionally limiting the
#' maximum number of clusters.
#'
#' @param mahal_df Data frame containing Mahalanobis distance results for one sample.
#' @param bandwidths Numeric vector of bandwidths to test.
#' @param max_clusters Maximum number of mean shift clusters allowed when
#'    choosing the best bandwidth. (Default: 7)
#' @param iterations Number of mean shift iterations. (Default: 500)
#' @param sample_name Optional sample name used in tuning output. (Default: NULL)
#' @param sample_col Column containing cell/sample barcodes. (Default:
#'   "sample_barcode")
#' @param rna_dist_col Column containing closest RNA Mahalanobis distance.
#'   (Default: "Mahal_Dist_RNA")
#' @param cnv_dist_col Column containing closest CNV Mahalanobis distance.
#'   (Default: "Mahal_Dist_CNV")
#'
#' @return A list with:
#' \describe{
#'   \item{\code{data}}{input data frame with mean shift Assignment added.}
#'   \item{\code{parameter_tuning}}{bandwidth tuning results.}
#'   }
#'
#' @export
run_mean_shift <- function(mahal_df,
                                 bandwidths = c(0.4, 0.5, 0.6, 0.7, 0.8),
                                 max_clusters = 7,
                                 iterations = 500,
                                 sample_name = NULL,
                                 sample_col = "sample_barcode",
                                 rna_dist_col = "Mahal_Dist_RNA",
                                 cnv_dist_col = "Mahal_Dist_CNV") {

  # --- Input validation ---
  mahal_df <- as.data.frame(mahal_df)

  required_cols <- c(sample_col, rna_dist_col, cnv_dist_col)
  missing_cols <- setdiff(required_cols, colnames(mahal_df))

  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if (anyNA(mahal_df[[sample_col]])) {
    stop(sprintf("%s contains NA values", sample_col))
  }

  if (anyDuplicated(mahal_df[[sample_col]]) > 0) {
    stop(sprintf("Duplicated values found in %s", sample_col))
  }

  if (!is.numeric(bandwidths) || length(bandwidths) == 0 ||
      anyNA(bandwidths) || any(!is.finite(bandwidths)) ||
      any(bandwidths <= 0)) {
    stop("bandwidths must be a positive numeric vector")
  }

  if (!is.numeric(max_clusters) || length(max_clusters) != 1 ||
      !is.finite(max_clusters) || max_clusters < 1) {
    stop("max_clusters must be a single positive integer")
  }

  if (!is.numeric(iterations) || length(iterations) != 1 ||
      !is.finite(iterations) || iterations < 1) {
    stop("iterations must be a single positive integer")
  }

  max_clusters <- as.integer(max_clusters)
  iterations <- as.integer(iterations)

  if (is.null(sample_name)) {
    sample_name <- NA_character_
  }

  # --- Create matrix to cluster on ---
  log2_rna <- log2(mahal_df[[rna_dist_col]] + 1)
  log2_cnv <- log2(mahal_df[[cnv_dist_col]] + 1)
  cell_ids <- mahal_df[[sample_col]]

  temp_dat <- cbind(RNA_M_Dist = log2_rna, CNV_M_Dist = log2_cnv)
  rownames(temp_dat) <- cell_ids

  # Remove cells with non-finite values
  keep_cells <- is.finite(temp_dat[, "RNA_M_Dist"]) &
    is.finite(temp_dat[, "CNV_M_Dist"])

  n_removed <- sum(!keep_cells)
  if (n_removed > 0) {
    warning(sprintf("Removed %d cells with non-finite distances", n_removed))
  }

  temp_dat <- temp_dat[keep_cells, , drop = FALSE]

  # --- If too few usable cells, return NA assignments ---
  if (nrow(temp_dat) < 3) {
    warning("Fewer than 3 usable cells. Skipping mean shift.")

    mahal_df$Assignment <- NA

    parameter_tuning <- data.frame(
      Samples = sample_name,
      Num_Clusters = NA_integer_,
      Bandwidth = bandwidths,
      Mean_Silh = NA_real_,
      Best = "No"
    )

    return(list(
      data = mahal_df,
      parameter_tuning = parameter_tuning
    ))
  }

  # --- Distance matrix for silhouette calculation ---
  dist_mat <- stats::dist(temp_dat)

  mean_shifts <- list()
  mean_silh <- rep(NA_real_, length(bandwidths))
  num_clusts <- rep(NA_integer_, length(bandwidths))

  # --- Run mean shift across bandwidths ---
  for (j in seq_along(bandwidths)) {

    mean_shifts[[j]] <- meanShiftR::meanShift(
      queryData = temp_dat,
      trainData = temp_dat,
      nNeighbors = nrow(temp_dat),
      bandwidth = rep(bandwidths[j], ncol(temp_dat)),
      iterations = iterations,
      alpha = 0
    )

    num_clusts[j] <- length(unique(mean_shifts[[j]]$assignment))

    # Silhouette is only meaningful when 1 < clusters < cells
    if (num_clusts[j] > 1 && num_clusts[j] < nrow(temp_dat)) {
      silh <- cluster::silhouette(
        x = mean_shifts[[j]]$assignment,
        dist = dist_mat
      )

      mean_silh[j] <- mean(silh[, "sil_width"])
    }
  }

  names(mean_shifts) <- paste0("Bandwidth_", bandwidths)

  # --- Choose best bandwidth ---
  chosen_vec <- rep("No", length(bandwidths))

  valid_models <- which(!is.na(mean_silh) & num_clusts <= max_clusters)

  if (length(valid_models) > 0) {
    model <- valid_models[which.max(mean_silh[valid_models])]
  } else if (!all(is.na(mean_silh))) {
    model <- which.max(mean_silh)
    warning(sprintf(
      "No bandwidth produced <= %d clusters. Using best silhouette.",
      max_clusters
    ))
  } else {
    model <- 1
    warning("All bandwidths produced single clusters. Using first bandwidth.")
  }

  chosen_vec[model] <- "Yes"

  parameter_tuning <- data.frame(
    Samples = sample_name,
    Num_Clusters = num_clusts,
    Bandwidth = bandwidths,
    Mean_Silh = mean_silh,
    Best = chosen_vec
  )

  # --- Assign clusters back to original data ---
  assignments <- rep(NA_integer_, nrow(mahal_df))
  names(assignments) <- cell_ids
  assignments[rownames(temp_dat)] <- mean_shifts[[model]]$assignment

  mahal_df$Assignment <- factor(
    assignments[cell_ids],
    levels = sort(unique(mean_shifts[[model]]$assignment))
  )

  list(
    data = mahal_df,
    parameter_tuning = parameter_tuning
  )
}
