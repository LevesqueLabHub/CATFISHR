#' Calculate Mahalanobis distance to reference cluster centroids
#'
#' Computes squared Mahalanobis distance from every cell to each reference
#' cluster centroid in a shared PCA space. All reference clusters use the same
#' number of PCs. Covariance singularity is handled by ridge regularization.
#'
#' @param pca_matrix A cells x PCs matrix of PCA embeddings.
#' @param clusters A vector of cluster assignments, one per cell. Order must
#'   match rows of pca_matrix, or cell_ids if provided.
#' @param ref_clusters Character or numeric vector of reference cluster IDs.
#'   The unassigned cluster "-1" is automatically excluded.
#' @param cell_ids Character vector of cell barcodes. If NULL, falls back to
#'   rownames of pca_matrix. (Default: NULL)
#' @param n_pcs Number of PCs to use. If NULL, uses select_pcs(). (Default: NULL)
#' @param min_cells Minimum cells required in a reference cluster. If NULL,
#'   automatically set to max(10, n_pcs + 2). (Default: NULL)
#' @param ridge_lambda Ridge multiplier added to covariance diagonal.
#'   (Default: 1e-6)
#' @param cov_tol Minimum reciprocal condition number accepted after ridge
#'   regularization. (Default: sqrt(.Machine$double.eps))
#'
#' @return A list with:
#' \describe{
#'   \item{\code{min_ref_distances}}{Named numeric vector of shortest squared Mahalanobis distance per cell.}
#'   \item{\code{all_distances}}{Cells by reference clusters squared distance matrix.}
#'   \item{\code{clusters}}{Named character vector of input cluster assignments.}
#'   \item{\code{is_reference}}{Named logical vector indicating whether each cell belongs to a reference cluster.}
#'   \item{\code{cluster_status}}{Named character vector: \code{"reference"} or \code{"query"}.}
#'   \item{\code{ref_clusters}}{Reference clusters actually used.}
#'   \item{\code{n_pcs}}{Number of PCs used.}
#' }
#'
#' @export
calc_mahalanobis <- function(pca_matrix,
                             clusters,
                             ref_clusters,
                             cell_ids = NULL,
                             n_pcs = NULL,
                             min_cells = NULL,
                             ridge_lambda = 1e-6,
                             cov_tol = sqrt(.Machine$double.eps)) {

  # --- Input validation: pca_matrix ---
  if (inherits(pca_matrix, "data.table")) {
    stop("pca_matrix is a data.table which does not support rownames. ",
         "Convert to matrix and provide cell_ids separately.")
  }

  if (is.data.frame(pca_matrix)) {
    pca_matrix <- as.matrix(pca_matrix)
  }

  if (!is.matrix(pca_matrix)) {
    stop("pca_matrix must be a matrix or data.frame")
  }

  if (!is.numeric(pca_matrix)) {
    stop("pca_matrix must contain numeric values")
  }

  if (anyNA(pca_matrix)) {
    stop("pca_matrix contains NA values")
  }

  if (any(!is.finite(pca_matrix))) {
    stop("pca_matrix contains non-finite values")
  }

  if (nrow(pca_matrix) == 0 || ncol(pca_matrix) == 0) {
    stop("pca_matrix has zero rows or columns")
  }

  # --- Resolve cell IDs ---
  if (is.null(cell_ids)) {
    cell_ids <- rownames(pca_matrix)

    if (is.null(cell_ids)) {
      stop("No cell_ids provided and pca_matrix has no rownames. ",
           "Provide cell_ids explicitly.")
    }
  }

  if (length(cell_ids) != nrow(pca_matrix)) {
    stop(sprintf(
      "cell_ids length (%d) does not match pca_matrix rows (%d)",
      length(cell_ids), nrow(pca_matrix)
    ))
  }

  if (anyNA(cell_ids)) {
    stop("cell_ids contains NA values")
  }

  if (any(duplicated(cell_ids))) {
    stop("cell_ids contains duplicates")
  }

  rownames(pca_matrix) <- cell_ids

  # --- Input validation: clusters ---
  if (length(clusters) != length(cell_ids)) {
    stop(sprintf(
      "clusters length (%d) does not match number of cells (%d)",
      length(clusters), length(cell_ids)
    ))
  }

  if (anyNA(clusters)) {
    stop("clusters contains NA values")
  }

  clusters <- as.character(clusters)
  names(clusters) <- cell_ids

  # --- Input validation: ref_clusters ---
  if (missing(ref_clusters) || is.null(ref_clusters) || length(ref_clusters) == 0) {
    stop("ref_clusters must contain at least one cluster ID")
  }

  if (anyNA(ref_clusters)) {
    warning("ref_clusters contains NA values; removing them")
    ref_clusters <- ref_clusters[!is.na(ref_clusters)]
  }

  ref_clusters <- setdiff(as.character(ref_clusters), "-1")

  if (length(ref_clusters) == 0) {
    stop("ref_clusters is empty or contains only the unassigned cluster (-1)")
  }

  available_refs <- intersect(ref_clusters, unique(clusters))
  missing_refs <- setdiff(ref_clusters, available_refs)

  if (length(available_refs) == 0) {
    stop(sprintf(
      "None of ref_clusters (%s) found in cluster assignments. Available: %s",
      paste(ref_clusters, collapse = ", "),
      paste(unique(clusters), collapse = ", ")
    ))
  }

  if (length(missing_refs) > 0) {
    warning(sprintf(
      "Reference clusters not found in data: %s",
      paste(missing_refs, collapse = ", ")
    ))
  }

  # --- Validate ridge parameters ---
  if (!is.numeric(ridge_lambda) || length(ridge_lambda) != 1 ||
      !is.finite(ridge_lambda) || ridge_lambda <= 0) {
    stop("ridge_lambda must be a single positive numeric value")
  }

  if (!is.numeric(cov_tol) || length(cov_tol) != 1 ||
      !is.finite(cov_tol) || cov_tol <= 0) {
    stop("cov_tol must be a single positive numeric value")
  }

  # --- Select PCs ---
  if (is.null(n_pcs)) {
    n_pcs <- select_pcs(pca_matrix)
    message(sprintf("Selected %d PCs by variance heuristic", n_pcs))
  }

  if (!is.numeric(n_pcs) || length(n_pcs) != 1 ||
      !is.finite(n_pcs) || n_pcs < 1) {
    stop("n_pcs must be a single positive integer")
  }

  if (n_pcs != as.integer(n_pcs)) {
    stop("n_pcs must be an integer")
  }

  if (n_pcs > ncol(pca_matrix)) {
    warning(sprintf(
      "n_pcs (%d) exceeds number of columns in pca_matrix (%d). Using %d PCs.",
      n_pcs, ncol(pca_matrix), ncol(pca_matrix)
    ))
    n_pcs <- ncol(pca_matrix)
  }

  n_pcs <- as.integer(n_pcs)

  # --- Set or validate min_cells after n_pcs is known ---
  if (is.null(min_cells)) {
    min_cells <- max(10L, n_pcs + 2L)
    message(sprintf(
      "min_cells was NULL; using min_cells = max(10, n_pcs + 2) = %d",
      min_cells
    ))
  } else {
    if (!is.numeric(min_cells) || length(min_cells) != 1 ||
        !is.finite(min_cells) || min_cells < 2) {
      stop("min_cells must be NULL or a single integer >= 2")
    }

    if (min_cells != as.integer(min_cells)) {
      stop("min_cells must be an integer")
    }

    min_cells <- as.integer(min_cells)
  }

  pca_sub <- pca_matrix[, seq_len(n_pcs), drop = FALSE]

  # --- Compute centroids and regularized covariance per reference cluster ---
  centroids <- list()
  covariances <- list()
  skipped <- character(0)

  for (ref in available_refs) {
    ref_cells <- cell_ids[clusters == ref]
    ref_data <- pca_sub[ref_cells, , drop = FALSE]

    if (nrow(ref_data) < min_cells) {
      warning(sprintf(
        "Reference cluster '%s' has %d cells, minimum required is %d. Skipping.",
        ref, nrow(ref_data), min_cells
      ))
      skipped <- c(skipped, ref)
      next
    }

    if (nrow(ref_data) <= n_pcs) {
      warning(sprintf(
        "Reference cluster '%s' has %d cells and %d PCs. Covariance is likely rank-deficient; ridge regularization will be used.",
        ref, nrow(ref_data), n_pcs
      ))
    }

    # Median centroid per PC
    centroids[[ref]] <- apply(ref_data, 2, stats::median)

    # Covariance in the same shared PCA space for every reference cluster
    cov_mat <- as.matrix(stats::cov(ref_data))

    if (anyNA(cov_mat) || any(!is.finite(cov_mat))) {
      warning(sprintf(
        "Reference cluster '%s': covariance contains NA or non-finite values. Skipping.",
        ref
      ))
      skipped <- c(skipped, ref)
      next
    }

    # Ridge regularization
    diag_mean <- mean(diag(cov_mat))

    if (!is.finite(diag_mean) || diag_mean <= 0) {
      diag_mean <- 1
    }

    lambda <- ridge_lambda * diag_mean
    cov_reg <- cov_mat + diag(lambda, nrow = n_pcs)

    # Increase ridge strength if still poorly conditioned
    ridge_multiplier <- 1

    while (base::rcond(cov_reg) < cov_tol && ridge_multiplier < 1e8) {
      ridge_multiplier <- ridge_multiplier * 10
      cov_reg <- cov_mat + diag(lambda * ridge_multiplier, nrow = n_pcs)
    }

    if (base::rcond(cov_reg) < cov_tol) {
      warning(sprintf(
        "Reference cluster '%s': covariance remains poorly conditioned after ridge regularization. Skipping.",
        ref
      ))
      skipped <- c(skipped, ref)
      next
    }

    covariances[[ref]] <- cov_reg
  }

  usable_refs <- setdiff(available_refs, skipped)

  if (length(usable_refs) == 0) {
    stop("No reference clusters usable after validation. ",
         "Check reference cluster cell counts, n_pcs, and covariance conditioning.")
  }

  # --- Calculate distances from all cells to each reference centroid ---
  all_distances <- matrix(
    NA_real_,
    nrow = nrow(pca_sub),
    ncol = length(usable_refs),
    dimnames = list(cell_ids, usable_refs)
  )

  for (ref in usable_refs) {
    all_distances[, ref] <- stats::mahalanobis(
      x = pca_sub,
      center = centroids[[ref]],
      cov = covariances[[ref]]
    )
  }

  # --- Extract minimum squared distance to any reference ---
  # --- Extract minimum squared distance to any reference ---
  min_ref_distances <- apply(all_distances, 1, min)

  names(min_ref_distances) <- cell_ids

  # --- Mark whether each cell belongs to a reference cluster ---
  is_reference <- clusters %in% usable_refs
  names(is_reference) <- cell_ids

  cluster_status <- ifelse(is_reference, "reference", "query")
  names(cluster_status) <- cell_ids

  list(
    min_ref_distances = min_ref_distances,
    all_distances = all_distances,
    clusters = clusters,
    is_reference = is_reference,
    cluster_status = cluster_status,
    ref_clusters = usable_refs,
    n_pcs = n_pcs
  )
}


#' Select number of PCs using variance drop-off heuristic
#'
#' @param pca_matrix A cells x PCs matrix of PCA embeddings.
#' @param drop_threshold Minimum percent variance drop between consecutive PCs.
#' @param default_pcs Number of PCs to use if no drop exceeds threshold.
#'
#' @return Integer number of PCs to use.
select_pcs <- function(pca_matrix,
                       drop_threshold = 0.1,
                       default_pcs = 30) {

  if (is.data.frame(pca_matrix)) {
    pca_matrix <- as.matrix(pca_matrix)
  }

  if (!is.matrix(pca_matrix) || !is.numeric(pca_matrix)) {
    stop("pca_matrix must be a numeric matrix or data.frame")
  }

  if (ncol(pca_matrix) < 2) {
    return(1L)
  }

  # Estimate variance explained from PCA embeddings
  pc_var <- apply(pca_matrix, 2, stats::var)

  if (anyNA(pc_var) || any(!is.finite(pc_var)) || sum(pc_var) <= 0) {
    warning("Could not estimate PC variance. Using default_pcs.")
    return(as.integer(min(default_pcs, ncol(pca_matrix))))
  }

  pct <- pc_var / sum(pc_var) * 100

  delta_pct <- pct[seq_len(length(pct) - 1)] - pct[2:length(pct)]

  hits <- which(delta_pct > drop_threshold)

  if (length(hits) == 0) {
    n_pcs <- min(default_pcs, ncol(pca_matrix))
  } else {
    n_pcs <- max(hits) + 1
  }

  as.integer(min(n_pcs, ncol(pca_matrix)))
}

