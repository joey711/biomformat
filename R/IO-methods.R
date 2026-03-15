################################################################################
################################################################################
#' Read a biom-format file, returning a \code{biom-class}.
#'
#' Import the data from a biom-format file into R, represented as an instance
#' of the \code{\link{biom-class}}; essentially a \code{\link{list}} with 
#' special constraints that map to \href{http://biom-format.org/documentation/biom_format.html}{the biom-format definition}.
#' 
#' The BIOM file format (canonically pronounced biome) is designed to be a general-use format for representing biological sample by observation contingency tables. BIOM is a recognized standard for the \href{http://www.earthmicrobiome.org/}{Earth Microbiome Project} and is a \href{http://gensc.org/}{Genomics Standards Consortium} candidate project. Please see \href{http://biom-format.org/}{the biom-format home page} for more details.
#' 
#' It is tempting to include an argument identifying the 
#' biom-format version number of the data file being imported.
#' However, the biom-format version number is a required
#' field in the biom-format definition. 
#' Rather than duplicate this formal specification
#' and allow the possibility of a conflict, the version 
#' number of the biom format will be referred to only by
#' the "format" field in the biom formatted data,
#' or its representation in R.
#'
#' @usage read_biom(biom_file)
#'
#' @param biom_file (Required). A character string indicating the 
#'  file location of the biom formatted file. This is a HDF5 or JSON formatted file
#'  specific to biological datasets. 
#'  The format is formally defined at \href{http://biom-format.org/documentation/biom_format.html}{the biom-format definition}
#'  and depends on the versioning.
#'
#' @return An instance of the \code{biom-class}.
#'
#' @seealso 
#' 
#' Function to create a biom object from R data,
#' \code{\link{make_biom}}.
#' 
#' Definition of the
#' \code{\link{biom-class}}. 
#' 
#' Function to write a biom format file from a biom object,
#' \code{\link{write_biom}}
#'
#' Accessor functions like \code{\link{header}}.
#'
#' @references \url{http://biom-format.org/}
#'
#' @importFrom jsonlite fromJSON
#' @export
#' @examples
#' # # # import with default parameters, specify a file
#' biom_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "biomformat")
#' biom_file
#' read_biom(biom_file)
#' biom_file <- system.file("extdata", "min_sparse_otu_table.biom", package = "biomformat")
#' biom_file
#' read_biom(biom_file)
#' ## The previous examples use system.file() because of constraints in specifying a fixed
#' ##   path within a reproducible example in a package. 
#' ## In practice, however, you can simply provide "hard-link"
#' ## character string path to your file:
#' # mybiomfile <- "path/to/my/biomfile.biom"
#' # read_biom(mybiomfile)
read_biom <- function(biom_file){
  # -- Detect file format by magic bytes (first 4 bytes) 
  # HDF5 files always begin with the 8-byte signature:
  #   \x89 H D F \r \n \x1a \n
  # Checking the first 4 bytes is sufficient and avoids reading large files.
  magic <- tryCatch(
    readBin(biom_file, what = "raw", n = 4L),
    error = function(e) raw(0)
  )
  is_hdf5 <- length(magic) == 4L &&
              magic[1] == as.raw(0x89) &&
              magic[2] == as.raw(0x48) &&   # 'H'
              magic[3] == as.raw(0x44) &&   # 'D'
              magic[4] == as.raw(0x46)      # 'F'

  if (is_hdf5) {
    # -- Route exclusively to HDF5 parser; never fall back to jsonlite --
    # This eliminates the confusing "lexical error" that appeared when an HDF5
    # file was accidentally parsed by jsonlite first.
    trash <- tryCatch(
      { x <- read_hdf5_biom(biom_file) },
      error = function(e) structure(class = "try-error",
                                    conditionMessage(e))
    )
    if (inherits(trash, "try-error")) {
      stop("File detected as HDF5 (BIOM v2) by magic bytes, but parsing failed.\n",
           "File: ", biom_file, "\n",
           "Error: ", as.character(trash))
    }
  } else {
    # -- Route exclusively to JSON parser (BIOM v1) --
    trash <- tryCatch(
      { x <- fromJSON(biom_file, simplifyDataFrame = FALSE, simplifyMatrix = FALSE) },
      error = function(e) structure(class = "try-error",
                                    conditionMessage(e))
    )
    if (inherits(trash, "try-error")) {
      stop("File detected as JSON (BIOM v1) but parsing failed.\n",
           "File: ", biom_file, "\n",
           "If this is an HDF5/BIOM-v2 file, ensure its magic bytes are intact.\n",
           "Error: ", as.character(trash))
    }
  }
  # Use the biom() constructor function to
  # instantiate a biom-class, perform validity checks. Return.
  return( biom(x) )
}
################################################################################
#' Write a biom-format v1 file, returning a \code{biom-class}.
#'
#' @param x (Required). A biom object that is going to be written to file
#'  as a proper biom formatted file, adhering to 
#'  \href{http://biom-format.org/documentation/biom_format.html}{the biom-format definition}.
#'  
#' @param biom_file (Required). A character string indicating the 
#'  file location of the biom formatted file. This is a JSON formatted file
#'  specific to biological datasets. 
#'  The format is formally defined at 
#'  \href{http://biom-format.org/documentation/biom_format.html}{the biom-format definition}
#'
#' @details
#' \code{write_biom()} serialises the entire BIOM object to a single JSON
#' string via \code{\link[jsonlite]{toJSON}}.  Because R character strings are
#' limited to \eqn{2^{31}-1} bytes, this function will fail with a
#' \emph{"character strings are limited to 2^31-1 bytes"} error for very
#' large datasets.  If your dataset has more than a few thousand samples or
#' features, use \code{\link{write_hdf5_biom}} instead — the HDF5 format has
#' no such size constraint.
#' @return Nothing. The first argument, \code{x}, is written to a file.
#'
#' @seealso 
#' 
#' Function to create a biom object from R data,
#' \code{\link{make_biom}}.
#' 
#' Definition of the
#' \code{\link{biom-class}}. 
#' 
#' The \code{\link{read_biom}} import function.
#'
#' Accessor functions like \code{\link{header}}.
#'
#' @references \url{http://biom-format.org/}
#'
#' @export
#' @importFrom jsonlite toJSON
#' @examples
#' biom_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "biomformat")
#' x = read_biom(biom_file)
#' outfile = tempfile()
#' write_biom(x, outfile)
#' y = read_biom(outfile)
#' identical(x, y) 
write_biom <- function(x, biom_file){
	cat(toJSON(x, always_decimal=TRUE, auto_unbox=TRUE), file=biom_file)
}
################################################################################
#' Write a biom object to an HDF5 (BIOM v2) file.
#'
#' Serialises a \code{\link{biom-class}} object to the
#' \href{http://biom-format.org/documentation/biom_format.html}{BIOM v2 HDF5 format}.
#' Both the sample-major and observation-major compressed-sparse representations
#' required by the spec are written, along with sample and observation metadata.
#'
#' The \code{rhdf5} package is required. If it is not installed a clear error
#' is thrown.  In normal use you should prefer \code{\link{write_biom}} for
#' JSON (BIOM v1) output; use this function only when HDF5 output is explicitly
#' needed.
#'
#' @param x (Required). A \code{\link{biom-class}} object.
#' @param biom_file (Required). Character string path to the output file.
#'   Any existing file at that path is overwritten.
#'
#' @return \code{biom_file} (invisibly).
#'
#' @seealso
#' \code{\link{read_hdf5_biom}}, \code{\link{write_biom}},
#' \code{\link{biom-class}}, \code{\link{make_biom}}.
#'
#' @references \url{http://biom-format.org/}
#'
#' @export
#' @importFrom utils packageVersion
#' @examples
#' biom_file <- system.file("extdata", "rich_sparse_otu_table.biom",
#'                           package = "biomformat")
#' x <- read_biom(biom_file)
#' outfile <- tempfile(fileext = ".biom")
#' write_hdf5_biom(x, outfile)
#' y <- read_biom(outfile)
#' identical(biom_data(x), biom_data(y))
write_hdf5_biom <- function(x, biom_file) {
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop(
      "The 'rhdf5' package is required to write HDF5/BIOM-v2 files.\n",
      "Install it with: BiocManager::install(\"rhdf5\")"
    )
  }

  # -- Extract the count matrix as a plain base R matrix --
  mat    <- as.matrix(biom_data(x))
  n_obs  <- nrow(mat)
  n_samp <- ncol(mat)

  # -- Build sample-major CCS (outer = samples = columns) --
  # Iterate columns; collect (obs_idx, val) pairs in column order.
  # Indices within each column are kept in ascending row order because
  # Matrix::sparseMatrix / which() already returns them sorted.
  s_indptr  <- integer(n_samp + 1L)
  s_indices <- integer(0)
  s_data    <- numeric(0)
  for (j in seq_len(n_samp)) {
    nz <- which(mat[, j] != 0)
    s_indptr[j + 1L] <- s_indptr[j] + length(nz)
    s_indices <- c(s_indices, nz - 1L)   # 0-based
    s_data    <- c(s_data,    mat[nz, j])
  }

  # -- Build observation-major CSR (outer = obs = rows) --
  o_indptr  <- integer(n_obs + 1L)
  o_indices <- integer(0)
  o_data    <- numeric(0)
  for (i in seq_len(n_obs)) {
    nz <- which(mat[i, ] != 0)
    o_indptr[i + 1L] <- o_indptr[i] + length(nz)
    o_indices <- c(o_indices, nz - 1L)   # 0-based
    o_data    <- c(o_data,    mat[i, nz])
  }

  # -- Create HDF5 file and group hierarchy --
  if (file.exists(biom_file)) file.remove(biom_file)
  rhdf5::h5createFile(biom_file)
  for (grp in c("observation", "observation/matrix",
                "sample",      "sample/matrix")) {
    rhdf5::h5createGroup(biom_file, grp)
  }

  # -- Write matrix datasets --
  rhdf5::h5write(rownames(mat),     biom_file, "observation/ids")
  rhdf5::h5write(o_indptr,          biom_file, "observation/matrix/indptr")
  rhdf5::h5write(o_indices,         biom_file, "observation/matrix/indices")
  rhdf5::h5write(as.numeric(o_data),biom_file, "observation/matrix/data")

  rhdf5::h5write(colnames(mat),     biom_file, "sample/ids")
  rhdf5::h5write(s_indptr,          biom_file, "sample/matrix/indptr")
  rhdf5::h5write(s_indices,         biom_file, "sample/matrix/indices")
  rhdf5::h5write(as.numeric(s_data),biom_file, "sample/matrix/data")

  # -- Write metadata --
  # Access raw per-element metadata directly from x$rows / x$columns to
  # preserve the original data structure (e.g. vector-valued taxonomy).
  # Each named metadata field is written as a dataset under the group:
  #   vector fields  -> 1-D array of length n_obs / n_samp
  #   list/vector    -> 2-D matrix [n_levels x n_obs/n_samp]  (taxonomy style)
  write_meta_group <- function(elements, group, n_elem) {
    first_meta <- elements[[1L]]$metadata
    if (is.null(first_meta) || length(first_meta) == 0L) return(invisible(NULL))
    rhdf5::h5createGroup(biom_file, paste0(group, "/metadata"))
    for (field in names(first_meta)) {
      vals <- lapply(elements, function(e) e$metadata[[field]])
      # Scalar per element -> write as 1-D character/numeric array
      if (all(lengths(vals) == 1L)) {
        rhdf5::h5write(unlist(vals), biom_file,
                       paste0(group, "/metadata/", field))
      } else {
        # Vector per element (e.g. taxonomy) -> write as [n_levels x n_elem] matrix
        mat_meta <- do.call(cbind, lapply(vals, as.character))
        rhdf5::h5write(mat_meta, biom_file,
                       paste0(group, "/metadata/", field))
      }
    }
    invisible(NULL)
  }

  write_meta_group(x$rows,    "observation", n_obs)
  write_meta_group(x$columns, "sample",      n_samp)

  # -- Write top-level HDF5 attributes --
  fid <- rhdf5::H5Fopen(biom_file)
  on.exit(rhdf5::H5Fclose(fid), add = TRUE)
  rhdf5::h5writeAttribute(
    as.character(if (is.null(x$id)) "No Table ID" else x$id),
    fid, "id")
  rhdf5::h5writeAttribute("http://biom-format.org",      fid, "format-url")
  rhdf5::h5writeAttribute(c(2L, 1L),                     fid, "format-version")
  rhdf5::h5writeAttribute(
    paste0("biomformat ", packageVersion("biomformat")),  fid, "generated-by")
  rhdf5::h5writeAttribute(as.character(Sys.time()),       fid, "creation-date")
  rhdf5::h5writeAttribute("OTU table",                    fid, "type")
  rhdf5::h5writeAttribute(as.integer(length(s_data)),     fid, "nnz")
  rhdf5::h5writeAttribute(as.integer(c(n_obs, n_samp)),   fid, "shape")

  invisible(biom_file)
}
################################################################################
#' Read in a biom-format v2 (HDF5) file, returning a \code{list}.
#'
#' This function reads a BIOM v2 HDF5 file directly. In normal use, you should
#' call \code{\link{read_biom}}, which automatically routes to this function
#' when it detects the HDF5 magic bytes. Call this function directly only when
#' you are certain the file is HDF5 format.
#'
#' The \code{rhdf5} package is required at runtime. If it is not available, or
#' if the underlying HDF5 system libraries are absent, a graceful R-level
#' warning is thrown instead of a fatal C-level abort.
#'
#' @param biom_file (Required). A character string path to an HDF5-format
#'  BIOM file.
#'
#' @return A named list representing the BIOM data, suitable for passing to
#'  \code{\link{biom}()}.
#'
#' @seealso
#' \code{\link{read_biom}}, \code{\link{make_biom}},
#' \code{\link{biom-class}}, \code{\link{write_biom}},
#' \code{\link{write_hdf5_biom}}.
#'
#' @references \url{http://biom-format.org/}
#'
#' @export
#' @examples
#' biom_file <- system.file("extdata", "rich_sparse_otu_table_hdf5.biom", package = "biomformat")
#' x = read_hdf5_biom(biom_file)
#' x = biom(x)
#' outfile = tempfile()
#' write_biom(x, outfile)
#' y = read_biom(outfile)
#' identical(observation_metadata(x),observation_metadata(y))
#' identical(sample_metadata(x),sample_metadata(y))
#' identical(biom_data(x), biom_data(y))
read_hdf5_biom <- function(biom_file){
  # -- Defensive check: require rhdf5 at runtime --
  # If rhdf5 is not installed, or its underlying HDF5 system libraries are
  # missing (common on stripped BBS nodes), emit a clear R-level warning and
  # stop gracefully rather than letting the C layer produce a fatal abort.
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    warning(
      "The 'rhdf5' package is required to read HDF5/BIOM-v2 files but is not ",
      "available. Please install it with:\n",
      "  BiocManager::install(\"rhdf5\")"
    )
    return(invisible(NULL))
  }

  # -- Attempt h5read; catch any C-level / system-library errors gracefully --
  result <- tryCatch({
    x <- rhdf5::h5read(biom_file, "/", read.attributes = TRUE)

    data    <- generate_matrix(x)
    rows    <- generate_metadata(x$observation)
    columns <- generate_metadata(x$sample)
    shape   <- c(length(data), length(data[[1]]))

    id                 <- attr(x, "id")
    vs                 <- attr(x, "format-version")
    format             <- sprintf("Biological Observation Matrix %s.%s", vs[1], vs[2])
    format_url         <- attr(x, "format-url")
    type               <- "OTU table"
    generated_by       <- attr(x, "generated-by")
    date               <- attr(x, "creation-date")
    matrix_type        <- "dense"
    matrix_element_type <- "int"

    namedList(id, format, format_url, type, generated_by, date,
              matrix_type, matrix_element_type, rows, columns, shape, data)
  },
  error = function(e) {
    warning(
      "Failed to read HDF5/BIOM-v2 file: ", biom_file, "\n",
      "This may indicate missing HDF5 system libraries (libhdf5) or a ",
      "corrupted file.\n",
      "Original error: ", conditionMessage(e)
    )
    return(invisible(NULL))
  })

  return(result)
}
################################################################################
