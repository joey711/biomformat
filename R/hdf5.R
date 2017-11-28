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
