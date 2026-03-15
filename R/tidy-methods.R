################################################################################
#' Coerce a biom object to a long-format tidy data frame
#'
#' Returns one row per (feature, sample) pair with columns
#' \code{feature_id}, \code{sample_id}, and \code{count} (or \code{value}
#' for character-type BIOM data), plus any sample and observation metadata
#' columns appended via a left-join.
#'
#' Column name conflicts between the base triplet and joined metadata are
#' resolved automatically by \code{\link{merge}} using the suffixes
#' \code{""} and \code{"_sample"} (for sample metadata) or
#' \code{"_feature"} (for observation metadata).
#'
#' The implementation is pure base R -- no tidyverse dependency.
#' For a tibble-wrapped version see \code{\link{as_tibble.biom}}.
#'
#' @param x A \code{\link{biom-class}} object.
#' @param ... Ignored; present for S3 generic compatibility.
#'
#' @return A \code{\link{data.frame}} in long (tidy) format with at least
#'   three columns: \code{feature_id}, \code{sample_id}, \code{count}.
#'   Additional columns from \code{\link{sample_metadata}} and
#'   \code{\link{observation_metadata}} are appended when present.
#'
#' @seealso \code{\link{biom_data}}, \code{\link{as_tibble.biom}}
#'
#' @export
#' @examples
#' biom_file <- system.file("extdata", "rich_dense_otu_table.biom",
#'                          package = "biomformat")
#' x <- read_biom(biom_file)
#' df <- as.data.frame(x)
#' head(df)
#' dim(df)   # 30 rows (5 features x 6 samples), 14 columns
as.data.frame.biom <- function(x, ...) {
  mat  <- as(biom_data(x), "matrix")
  smd  <- sample_metadata(x)
  omd  <- observation_metadata(x)

  # Pivot to long form using base R (no tidyr dependency).
  # as.data.frame(as.table(mat)) gives one row per cell with columns
  # Var1 (feature), Var2 (sample), Freq (value).
  long <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(long) <- c("feature_id", "sample_id", "count")

  # Left-join sample metadata keyed on sample_id
  if (!is.null(smd)) {
    smd_df <- data.frame(sample_id  = rownames(smd),
                         smd,
                         check.names    = FALSE,
                         stringsAsFactors = FALSE)
    long <- merge(long, smd_df,
                  by = "sample_id", all.x = TRUE,
                  suffixes = c("", "_sample"))
  }

  # Left-join observation/feature metadata keyed on feature_id
  if (!is.null(omd)) {
    omd_df <- data.frame(feature_id = rownames(omd),
                         omd,
                         check.names    = FALSE,
                         stringsAsFactors = FALSE)
    long <- merge(long, omd_df,
                  by = "feature_id", all.x = TRUE,
                  suffixes = c("", "_feature"))
  }

  long
}

################################################################################
#' Coerce a biom object to a tidy tibble
#'
#' A thin wrapper around \code{\link{as.data.frame.biom}} that returns a
#' \code{\link[tibble]{tibble}} instead of a plain \code{data.frame}.
#' Requires the \pkg{tibble} package (in \code{Suggests}).
#'
#' @param x A \code{\link{biom-class}} object.
#' @param ... Passed to \code{\link[tibble]{as_tibble}}.
#'
#' @return A \code{\link[tibble]{tibble}} in long (tidy) format.
#'   See \code{\link{as.data.frame.biom}} for column details.
#'
#' @seealso \code{\link{as.data.frame.biom}}
#'
#' @export
#' @examples
#' biom_file <- system.file("extdata", "rich_dense_otu_table.biom",
#'                          package = "biomformat")
#' x <- read_biom(biom_file)
#' if (requireNamespace("tibble", quietly = TRUE)) {
#'   # Call the method directly to ensure correct S3 dispatch
#'   tbl <- as_tibble.biom(x)
#'   tbl
#' }
as_tibble.biom <- function(x, ...) {
  if (!requireNamespace("tibble", quietly = TRUE))
    stop("tibble is required. Install with: install.packages('tibble')")
  tibble::as_tibble(as.data.frame(x), ...)
}
################################################################################
