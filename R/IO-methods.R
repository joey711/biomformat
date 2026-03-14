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
#' 
#' \code{\link{read_biom}}, \code{\link{make_biom}},
#' \code{\link{biom-class}}, \code{\link{write_biom}}.
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
