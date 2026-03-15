################################################################################
#' Convert a biom object to a SummarizedExperiment
#'
#' Converts a \code{\link{biom-class}} object into a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}}, placing the
#' count/value matrix in \code{assay(., "counts")}, sample metadata in
#' \code{colData()}, and observation/feature metadata in \code{rowData()}.
#'
#' Both \code{colData} and \code{rowData} are coerced to
#' \code{\link[S4Vectors]{DataFrame}}.  When a biom object carries no sample or
#' observation metadata (i.e. the corresponding accessor returns \code{NULL}),
#' an empty \code{DataFrame} with the correct row names is substituted so that
#' the resulting \code{SummarizedExperiment} is always fully valid.
#'
#' @param x A \code{\link{biom-class}} object.
#'
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} with:
#' \describe{
#'   \item{\code{assay("counts")}}{The feature-by-sample count (or value)
#'     matrix returned by \code{\link{biom_data}}.}
#'   \item{\code{colData()}}{Per-sample metadata from
#'     \code{\link{sample_metadata}}, or an empty \code{DataFrame}.}
#'   \item{\code{rowData()}}{Per-feature metadata from
#'     \code{\link{observation_metadata}}, or an empty \code{DataFrame}.}
#' }
#'
#' @seealso
#' \code{\link{read_biom}}, \code{\link{biom_data}},
#' \code{\link{sample_metadata}}, \code{\link{observation_metadata}}
#'
#' @importFrom methods setAs
#' @export
#'
#' @examples
#' biom_file <- system.file("extdata", "rich_sparse_otu_table.biom",
#'                          package = "biomformat")
#' x <- read_biom(biom_file)
#' if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'   se <- biom_to_SummarizedExperiment(x)
#'   se
#'   # S4 coercion syntax also works:
#'   se2 <- as(x, "SummarizedExperiment")
#'   identical(SummarizedExperiment::assay(se), SummarizedExperiment::assay(se2))
#' }
biom_to_SummarizedExperiment <- function(x) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
    stop("SummarizedExperiment is required. ",
         "Install with: BiocManager::install('SummarizedExperiment')")
  if (!requireNamespace("S4Vectors", quietly = TRUE))
    stop("S4Vectors is required. ",
         "Install with: BiocManager::install('S4Vectors')")

  mat <- biom_data(x)
  smd <- sample_metadata(x)
  omd <- observation_metadata(x)

  # colData: per-sample metadata - must be a DataFrame with sample names as rows
  col_data <- if (is.null(smd)) {
    S4Vectors::DataFrame(row.names = colnames(mat))
  } else {
    S4Vectors::DataFrame(smd)
  }

  # rowData: per-feature metadata - must be a DataFrame with feature names as rows
  row_data <- if (is.null(omd)) {
    S4Vectors::DataFrame(row.names = rownames(mat))
  } else {
    S4Vectors::DataFrame(omd)
  }

  SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = mat),
    colData = col_data,
    rowData = row_data
  )
}

# S4 coercion: as(x, "SummarizedExperiment")
# Registered only when SummarizedExperiment is available at load time.
# The .onLoad hook is the correct place for conditional S4 registration, but
# since the biom class is defined before SE is even loaded, using setAs() inside
# an if() block at namespace-load time is the accepted Bioconductor pattern.
if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  setAs("biom", "SummarizedExperiment", function(from) {
    biom_to_SummarizedExperiment(from)
  })
}
################################################################################
