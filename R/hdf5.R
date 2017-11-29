################################################################################
################################################################################
# Read HDF5, e.g. biom-format v2.0 & v2.1
################################################################################
################################################################################
# Generic internal function for generating the count matrix.
#' @keywords internal
generate_matrix <- function(x){
  indptr  = x$sample$matrix$indptr+1
  indices = x$sample$matrix$indices+1
  data    = x$sample$matrix$data
  nr = length(x$observation$ids)
  
  counts = sapply(2:length(indptr),function(i){
    x = rep(0,nr)
    seq = indptr[i-1]:(indptr[i]-1)
    x[indices[seq]] = data[seq]
    x
  })
  rownames(counts) = x$observation$ids
  colnames(counts) = x$sample$ids
  # I wish this next line wasn't necessary
  lapply(1:nrow(counts),function(i){
    counts[i,]
  })
}
################################################################################
# Generic internal function for generating the metadata.
#' @keywords internal
generate_metadata <- function(x){
  metadata = x$metadata
  metadata = lapply(1:length(x$ids),function(i){
    id_metadata = lapply(metadata,function(j){
      if(length(dim(j))>1){ as.vector(j[,i,drop=FALSE]) }
      else{ j[i] }
    })
    list(id = x$ids[i],metadata=id_metadata)
  })
  return(metadata)
}
################################################################################
#' Read in a biom-format vs 2 file, returning a \code{list}.
#'
#' This function is meant only to be used if the user knows the file is
#' a particular version / hdf5 format. Otherwise, the `read_biom` file should be used.
#'
#' @param biom_file (Required). A biom object that is going to be written to file
#'  as a proper biom formatted file, adhering to 
#'  \href{http://biom-format.org/documentation/biom_format.html}{the biom-format definition}.
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
#' The \code{\link{read_hdf5_biom}} import function.
#'
#' Accessor functions like \code{\link{header}}.
#'
#' @references \url{http://biom-format.org/}
#'
#' @export
#' @importFrom rhdf5 h5read
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
read_hdf5_biom<-function(biom_file){
  x = h5read(biom_file,"/",read.attributes = TRUE)
  data = generate_matrix(x)
  rows = generate_metadata(x$observation)
  columns = generate_metadata(x$sample)
  shape = c(length(data),length(data[[1]]))
  
  id = attr(x,"id")
  vs = attr(x,"format-version")
  format = sprintf("Biological Observation Matrix %s.%s",vs[1],vs[2])
  format_url = attr(x,"format-url")
  type = "OTU table"
  generated_by = attr(x,"generated-by")
  date = attr(x,"creation-date")
  matrix_type = "dense"
  matrix_element_type = "int"
  
  namedList(id,format,format_url,type,generated_by,date,matrix_type,matrix_element_type,
            rows,columns,shape,data)
}
################################################################################
################################################################################
# Write HDF5
################################################################################
################################################################################
# Writing requires coercion to v2.1 matrix format.
# Oddly, this is defined redundantly as both transpositions of the matrix 
# in the format simultaneously, rather than just picking one as a convention.
#
#
#' Access compressed sparse column format (CSC) from a matrix
#' 
#' This relies heavily on the compressed sparse format 
#' defined in the [Matrix::Matrix] package [Matrix::CsparseMatrix-class].
#' 
#' See a
#' [python doc on CSC](http://www.scipy-lectures.org/advanced/scipy_sparse/csc_matrix.html)
#' 
#' @param M A [matrix], or matrix-like coercible to [Matrix::CsparseMatrix-class].
#' @param slotName The name ([chracter]) of the CsparseMatrix slot.
#' 
#' @return The values corresponding to the indicated compressed sparse column feature.
#' 
#' @importClassesFrom Matrix CsparseMatrix
#' 
#' @export
#' 
#' @examples 
#' mat1 = matrix(c(1,0,2,0,0,3,4,5,6), nrow = 3, byrow = TRUE)
#' access_sparse(mat1, "p")
#' sparse_indptr(mat1)
#' sparse_indices(mat1)
#' sparse_data(mat1)
access_sparse = function(M, slotName){
  slot(as(M, "CsparseMatrix"), slotName)
}
#' @describeIn access_sparse Access the index pointer, indptr
#' @export
sparse_indptr = function(M){
  access_sparse(M, "p")
}
#' @describeIn access_sparse Access the indices
#' @export
sparse_indices = function(M){
  access_sparse(M, "i")
}
#' @describeIn access_sparse Access the sparse data entries
#' @export
sparse_data = function(M){
  access_sparse(M, "x")
}
#' Transform matrix to biom-format list representation of Compressed Sparse Column matrix
#' 
#' @param m A base [matrix].
#' 
#' @seealso
#' A [python doc on CSC](http://www.scipy-lectures.org/advanced/scipy_sparse/csc_matrix.html)
#' 
#' @importClassesFrom Matrix CsparseMatrix
#' 
#' @export
#' 
matrix_to_biom_csc = function(m){
  # Coerce up-front to avoid repeating coercion at each step
  M = as(m, "CsparseMatrix")
  list(
    data = sparse_data(M),
    indices = sparse_indices(M),
    indptr = sparse_indptr(M)
  )
}
#' Transform matrix to biom format (sparse) matrix list
#' 
#' Assumes that rows are observations and columns are samples.
#' 
#' @param m The matrix that you want to coerce to a biom-like list.
#' 
#' @importClassesFrom Matrix CsparseMatrix
#' @importFrom Matrix t
#' 
#' @export
#' 
matrix_to_biom_matlist = function(m){
  # Coerce up-front to avoid repeating coercion at each step
  M = as(m, "CsparseMatrix")
  list(
    # biom format requires (redundantly) shape and nnz
    shape = dim(M),
    nnz = length(sparse_data(M)),
    # Rows are observations are taxa
    observation = list(
      ids = rownames(M),
      matrix = matrix_to_biom_csc(t(M)),
      # Add required metadata and group-metadata groups
      metadata = list(),
      `group-metadata` = list()
    ),
    # Columns are samples
    sample = list(
      ids = colnames(M),
      matrix = matrix_to_biom_csc(M),
      # Add required metadata and group-metadata groups
      metadata = list(),
      `group-metadata` = list()
    )
  )
}
#
#' Acceptable values for biom v2.1 top-level `type`
#' 
#' @export
#' 
biomTypes2.1 = c(
  "Taxon table",
  "OTU table",
  "Pathway table",
  "Function table",
  "Ortholog table",
  "Gene table",
  "Metabolite table")
#
#' Generate required top-level attributes for biom version 2.1
#' 
#' See [the format definition](http://biom-format.org/documentation/format_versions/biom-2.1.html).
#' Note that this excludes the two matrix-dependent top-level attributes,
#' `shape` and `nnz`.
#' 
#' @param id Optional table identifier string.
#' @param type String indicating the table type, from a controlled vocabulary, [`biomTypes2.1`].
#' @param format_url The URL string for the biom-format project.
#' @param format_version The biom-format version tuple.
#' @param generated_by The package and version that generated the corresponding biom file.
#' @param creation_date The date the file was built, in ISO 8601 format.
#' 
#' @export
#' 
biom_header_2.1 = function(
  id = "",
  type = biomTypes2.1,
  format_url = "http://biom-format.org",
  format_version = c(2L, 1L),
  generated_by = sprintf("biomformat %s", packageVersion("biomformat")),
  creation_date = strptime(x = Sys.time(), format = "%Y-%m-%d %H:%M:%S")
){
  list(
    id = id,
    type = type[1],
    `format-url` = format_url,
    `format-version` = format_version,
    `generated-by` = generated_by,
    `creation-date` = as.character(creation_date)
  )
}
# From:
# http://biom-format.org/documentation/format_versions/biom-2.1.html
#
# Required top-level attributes:
#   
# id                   : <string or null> a field that can be used to id a table (or null)
# type                 : <string> Table type (a controlled vocabulary)
# Acceptable values:
#   "OTU table"
# "Pathway table"
# "Function table"
# "Ortholog table"
# "Gene table"
# "Metabolite table"
# "Taxon table"
# format-url           : <url> A string with a static URL providing format details
# format-version       : <tuple> The version of the current biom format, major and minor
# generated-by         : <string> Package and revision that built the table
# creation-date        : <datetime> Date the table was built (ISO 8601 format)
# shape                : <list of ints>, the number of rows and number of columns in data
# nnz                  : <int> The number of non-zero elements in the table
# 
# Required groups:
#
# observation/               : The HDF5 group that contains observation specific information and an observation oriented view of the data
# observation/matrix         : The HDF5 group that contains matrix data oriented for observation-wise operations (e.g., in compressed sparse row format)
# observation/metadata       : The HDF5 group that contains observation specific metadata information
# observation/group-metadata : The HDF5 group that contains observation specific group metadata information (e.g., phylogenetic tree)
# sample/                    : The HDF5 group that contains sample specific information and a sample oriented data oriented view of the data
# sample/matrix              : The HDF5 group that contains matrix data oriented for sample-wise operations (e.g., in compressed sparse column format)
# sample/metadata            : The HDF5 group that contains sample specific metadata information
# sample/group-metadata      : The HDF5 group that contains sample specific group metadata information (e.g., relationships between samples)
#  
# Required datasets:
#   
# observation/ids            : <string> or <variable length string> A (N,) dataset of the observation IDs, where N is the total number of IDs
# observation/matrix/data    : <float64> A (nnz,) dataset containing the actual matrix data
# observation/matrix/indices : <int32> A (nnz,) dataset containing the column indices (e.g., maps into samples/ids)
# observation/matrix/indptr  : <int32> A (M+1,) dataset containing the compressed row offsets
# sample/ids                 : <string> or <variable length string> A (M,) dataset of the sample IDs, where M is the total number of IDs
# sample/matrix/data         : <float64> A (nnz,) dataset containing the actual matrix data
# sample/matrix/indices      : <int32> A (nnz,) dataset containing the row indices (e.g., maps into observation/ids)
# sample/matrix/indptr       : <int32> A (N+1,) dataset containing the compressed column offsets

#' The top-level attributes required according to biom-format 2.1 specification
#' 
#' See [the format definition](http://biom-format.org/documentation/format_versions/biom-2.1.html).
#' 
#' @export
#' 
requiredAttributes2.1 <-
  c("id", "type", "format-url", "format-version", "generated-by", "creation-date", "shape", "nnz")
#' Convert matrix to biom-format 2.1 as an R list
#' 
#' See [The version 2.1 specification](http://biom-format.org/documentation/format_versions/biom-2.1.html)
#' for more details.
#' 
#' @inheritParams access_sparse
#' @param ... Additional named arguments passed to [biom_header_2.1()].
#' 
#' @seealso 
#' 
#' [biom_header_2.1()]
#' 
#' [matrix_to_biom_matlist()]
#' 
#' @export
#' 
matrix_to_biom_2.1_list = function(M, ...){
  c(
    biom_header_2.1(...),
    matrix_to_biom_matlist(M)
  )
}
