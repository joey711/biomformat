library(rhdf5)
outfile <- "inst/extdata/zero_col_hdf5.biom"
if (file.exists(outfile)) file.remove(outfile)
h5createFile(outfile)
for (grp in c("observation","observation/matrix","sample","sample/matrix")) {
  h5createGroup(outfile, grp)
}
# Observation IDs
h5write(c("OTU1","OTU2","OTU3"), outfile, "observation/ids")
# Sample IDs
h5write(c("Samp1","ZeroSamp","Samp3"), outfile, "sample/ids")

# Count matrix (Samp1=[1,3,5], ZeroSamp=[0,0,0], Samp3=[2,4,6])
# Sample-major CCS:
#   Samp1 has non-zeros at obs 0,1,2 with vals 1,3,5
#   ZeroSamp has no non-zeros
#   Samp3 has non-zeros at obs 0,1,2 with vals 2,4,6
# indptr: [0, 3, 3, 6]  (3 nnz in Samp1, 0 in ZeroSamp, 3 in Samp3)
# indices: [0,1,2, 0,1,2]  (0-based obs indices)
# data: [1,3,5, 2,4,6]
h5write(as.integer(c(0L,3L,3L,6L)), outfile, "sample/matrix/indptr")
h5write(as.integer(c(0L,1L,2L, 0L,1L,2L)), outfile, "sample/matrix/indices")
h5write(as.numeric(c(1,3,5, 2,4,6)), outfile, "sample/matrix/data")

# Observation-major CSR:
#   OTU1 has non-zeros at samples 0,2 with vals 1,2
#   OTU2 has non-zeros at samples 0,2 with vals 3,4
#   OTU3 has non-zeros at samples 0,2 with vals 5,6
# indptr: [0, 2, 4, 6]
# indices: [0,2, 0,2, 0,2]
# data: [1,2, 3,4, 5,6]
h5write(as.integer(c(0L,2L,4L,6L)), outfile, "observation/matrix/indptr")
h5write(as.integer(c(0L,2L, 0L,2L, 0L,2L)), outfile, "observation/matrix/indices")
h5write(as.numeric(c(1,2, 3,4, 5,6)), outfile, "observation/matrix/data")

# HDF5 attributes
fid <- H5Fopen(outfile)
h5writeAttribute("test_zero_col", fid, "id")
h5writeAttribute("http://biom-format.org", fid, "format-url")
h5writeAttribute(c(2L,1L), fid, "format-version")
h5writeAttribute("biomformat test fixture", fid, "generated-by")
h5writeAttribute(as.character(Sys.time()), fid, "creation-date")
h5writeAttribute("OTU table", fid, "type")
h5writeAttribute(6L, fid, "nnz")
h5writeAttribute(c(3L,3L), fid, "shape")
H5Fclose(fid)
cat("Created", outfile, "\n")
