################################################################################
# Use testthat to test file import and resulting class (and values)
################################################################################
library("biomformat"); library("testthat")
packageVersion("biomformat")
# # # # TESTS!
min_dense_file    = system.file("extdata", "min_dense_otu_table.biom",        package = "biomformat")
min_sparse_file   = system.file("extdata", "min_sparse_otu_table.biom",       package = "biomformat")
rich_dense_file   = system.file("extdata", "rich_dense_otu_table.biom",       package = "biomformat")
rich_sparse_file  = system.file("extdata", "rich_sparse_otu_table.biom",      package = "biomformat")
rich_dense_char   = system.file("extdata", "rich_dense_char.biom",            package = "biomformat")
rich_sparse_char  = system.file("extdata", "rich_sparse_char.biom",           package = "biomformat")
hdf5_file         = system.file("extdata", "rich_sparse_otu_table_hdf5.biom", package = "biomformat")
# Read the biom-format files
x1 = read_biom(min_dense_file)
x2 = read_biom(min_sparse_file)
x3 = read_biom(rich_dense_file)
x4 = read_biom(rich_sparse_file)
x5 = read_biom(rich_dense_char)
x6 = read_biom(rich_sparse_char)



# # # # TESTS!

test_that("Classes are all biom", {
	expect_is(x1, "biom")
	expect_is(x2, "biom")
	expect_is(x3, "biom")
	expect_is(x4, "biom")
	expect_is(x5, "biom")
	expect_is(x6, "biom")
})

test_that("min/rich files have same biom_data", {
	expect_identical(biom_data(x2), biom_data(x4))
	expect_identical(biom_data(x1), biom_data(x3))
})

test_that("biom_datas can be manipulated mathematically", {
	expect_identical(2*biom_data(x2), 4*biom_data(x4)/2)
	expect_identical(2*biom_data(x1)-biom_data(x1), biom_data(x3))
})

test_that("empty stuff is NULL", {
	expect_is(sample_metadata(x1), "NULL")
	expect_is(sample_metadata(x1, 2:4), "NULL")
	expect_is(observation_metadata(x1), "NULL")
})

test_that("Expected classes of non-empty components", {
	expect_is(observation_metadata(x3), "data.frame")
	expect_is(observation_metadata(x3, 2:4), "data.frame")
	expect_is(observation_metadata(x3, 3), "data.frame")
	expect_is(observation_metadata(x1), "NULL")
	expect_true(is(sample_metadata(x3), "data.frame"))
	expect_true(is(biom_data(x3), "Matrix"))
	expect_true(is(header(x3), "list"))
})

test_that("imported biom files are S4", {
	expect_true(isS4(x1))
	expect_true(isS4(x2))
	expect_true(isS4(x3))
	expect_true(isS4(x4))
})

test_that("show method output tests",{
	expect_output(print(x1), "biom object. \ntype:")
	expect_output(print(x4), "biom object. \ntype:")
})

# --- Milestone 3: drop=FALSE regression (PR #11 / PR #12) --------------------
# Sparse BIOM: single-row or single-column subset must retain matrix dimensions
# and must carry correct row/col names.
test_that("sparse single-row subset retains matrix dimensions (drop=FALSE, PR #11)", {
	# x2 is a sparse BIOM (min_sparse_otu_table): 5 rows x 6 cols
	r1 <- biom_data(x2, 1L, 1:6L)
	expect_false(is.null(dim(r1)),
	             info = "single-row sparse subset must have dim() (not a plain vector)")
	expect_equal(nrow(r1), 1L)
	expect_equal(ncol(r1), 6L)
	expect_equal(length(rownames(r1)), 1L)
	expect_equal(length(colnames(r1)), 6L)
})

test_that("sparse single-column subset retains matrix dimensions (drop=FALSE, PR #12)", {
	c1 <- biom_data(x2, 1:5L, 1L)
	expect_false(is.null(dim(c1)),
	             info = "single-col sparse subset must have dim() (not a plain vector)")
	expect_equal(nrow(c1), 5L)
	expect_equal(ncol(c1), 1L)
	expect_equal(length(rownames(c1)), 5L)
	expect_equal(length(colnames(c1)), 1L)
})

test_that("sparse single-cell subset retains matrix dimensions (drop=FALSE)", {
	s1 <- biom_data(x2, 1L, 1L)
	expect_false(is.null(dim(s1)),
	             info = "single-cell sparse subset must have dim()")
	expect_equal(nrow(s1), 1L)
	expect_equal(ncol(s1), 1L)
})

test_that("full sparse matrix still returns correctly after drop=FALSE fix", {
	m <- biom_data(x2)
	expect_false(is.null(dim(m)))
	expect_equal(nrow(m), 5L)
	expect_equal(ncol(m), 6L)
})

# --- Milestone 3: magic-byte HDF5 routing regression -------------------------
test_that("read_biom routes HDF5 file via magic bytes without jsonlite error", {
	# Must parse cleanly -- no 'lexical error' from jsonlite
	expect_error(read_biom(hdf5_file), NA)
	xh <- read_biom(hdf5_file)
	expect_is(xh, "biom")
})

test_that("read_biom correctly classifies JSON vs HDF5 fixtures via magic bytes", {
	# JSON fixtures must NOT be misidentified as HDF5
	expect_is(read_biom(min_dense_file),  "biom")
	expect_is(read_biom(min_sparse_file), "biom")
	expect_is(read_biom(rich_dense_file), "biom")
	# HDF5 fixture must parse without error
	expect_is(read_biom(hdf5_file),       "biom")
})
