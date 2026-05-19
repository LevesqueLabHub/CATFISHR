#' Format RNA and CNV Mahalanobis distance outputs
#'
#' Combines the minimum Mahalanobis distances from RNA and CNV
#' `calc_mahalanobis()` outputs into a single data frame for downstream
#' mean shift clustering and plotting.
#'
#' @param mahal_RNA Output list from [calc_mahalanobis()] for RNA PCA space.
#' @param mahal_CNV Output list from [calc_mahalanobis()] for CNV PCA space.
#'   outputs are kept. If `TRUE`, all cells from either output are kept and
#'   missing values are filled with `NA`. Default `FALSE`.
#'
#' @return A data frame with one row per cell and columns:
#' \describe{
#'   \item{sample_barcode}{Cell barcode.}
#'   \item{Mahal_Dist_RNA}{Minimum squared Mahalanobis distance in RNA space.}
#'   \item{Mahal_Dist_CNV}{Minimum squared Mahalanobis distance in CNV space.}
#'   \item{RNA_ref_status}{If cell in RNA reference cluster.}
#'   \item{CNV_ref_status}{If cell in CNV reference cluster.}
#'   \item{RNA_ori_cluster}{Original RNA cluster assignments.}
#'   \item{CNV_ori_cluster}{Original CNV cluster assignments.}
#' }
#'
#' @export
format_mahal_output <- function(mahal_RNA, mahal_CNV) {

  common_cells <- intersect(
    names(mahal_RNA$min_ref_distances),
    names(mahal_CNV$min_ref_distances)
  )

  if (length(common_cells) == 0) {
    stop("No overlapping cells between RNA and CNV outputs")
  }

  data.frame(
    sample_barcode = common_cells,
    Mahal_Dist_RNA = unname(mahal_RNA$min_ref_distances[common_cells]),
    Mahal_Dist_CNV = unname(mahal_CNV$min_ref_distances[common_cells]),
    RNA_ref_status = unname(mahal_RNA$cluster_status[common_cells]),
    CNV_ref_status = unname(mahal_CNV$cluster_status[common_cells]),
    RNA_ori_cluster = unname(mahal_RNA$clusters[common_cells]),
    CNV_ori_cluster = unname(mahal_CNV$clusters[common_cells]),
    stringsAsFactors = FALSE
  )
}
