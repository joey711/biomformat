################################################################################
# Tests for write_hdf5_biom() (M5 / v1.39.9)
################################################################################
library("biomformat"); library("testthat")

skip_if_not_installed("rhdf5")

json_file  <- system.file("extdata", "rich_sparse_otu_table.biom",       package = "biomformat")
hdf5_file  <- system.file("extdata", "rich_sparse_otu_table_hdf5.biom",  package = "biomformat")
min_file   <- system.file("extdata", "min_sparse_otu_table.biom",        package = "biomformat")
x_json     <- read_biom(json_file)
x_hdf5_ref <- read_biom(hdf5_file)
x_min      <- read_biom(min_file)

test_that("write_hdf5_biom() returns output path invisibly", {
  outf <- tempfile(fileext = ".biom")
  result <- write_hdf5_biom(x_min, outf)
  expect_equal(result, outf)
  expect_true(file.exists(outf))
})

test_that("JSON -> HDF5 -> biom round-trip preserves count matrix", {
  outf <- tempfile(fileext = ".biom")
  write_hdf5_biom(x_json, outf)
  y <- read_biom(outf)
  expect_equal(as.matrix(biom_data(x_json)),
               as.matrix(biom_data(y)))
})

test_that("JSON -> HDF5 -> biom round-trip preserves row/col names", {
  outf <- tempfile(fileext = ".biom")
  write_hdf5_biom(x_json, outf)
  y <- read_biom(outf)
  expect_equal(rownames(x_json), rownames(y))
  expect_equal(colnames(x_json), colnames(y))
})

test_that("min biom (no metadata) round-trips cleanly", {
  outf <- tempfile(fileext = ".biom")
  write_hdf5_biom(x_min, outf)
  y <- read_biom(outf)
  expect_equal(as.matrix(biom_data(x_min)),
               as.matrix(biom_data(y)))
  expect_null(observation_metadata(y))
  expect_null(sample_metadata(y))
})

test_that("output file is valid HDF5 (magic bytes)", {
  outf <- tempfile(fileext = ".biom")
  write_hdf5_biom(x_json, outf)
  magic <- readBin(outf, what = "raw", n = 4L)
  expect_equal(magic[1], as.raw(0x89))
  expect_equal(magic[2], as.raw(0x48))   # 'H'
  expect_equal(magic[3], as.raw(0x44))   # 'D'
  expect_equal(magic[4], as.raw(0x46))   # 'F'
})

test_that("write_hdf5_biom() overwrites existing file without error", {
  outf <- tempfile(fileext = ".biom")
  write_hdf5_biom(x_min, outf)
  expect_no_error(write_hdf5_biom(x_min, outf))
})

test_that("JSON sample metadata survives HDF5 round-trip", {
  outf <- tempfile(fileext = ".biom")
  write_hdf5_biom(x_json, outf)
  y <- read_biom(outf)
  smd_x <- sample_metadata(x_json)
  smd_y <- sample_metadata(y)
  expect_false(is.null(smd_y))
  expect_equal(rownames(smd_x), rownames(smd_y))
  # BODY_SITE column must be preserved
  expect_equal(as.character(smd_x[, "BODY_SITE"]),
               as.character(smd_y[, "BODY_SITE"]))
})

test_that("missing rhdf5 gives a clear error", {
  # Temporarily mock requireNamespace to pretend rhdf5 is absent
  local({
    outf <- tempfile(fileext = ".biom")
    with_mocked_bindings(
      requireNamespace = function(pkg, ...) if (pkg == "rhdf5") FALSE else TRUE,
      .package = "base",
      expect_error(write_hdf5_biom(x_min, outf), "rhdf5")
    )
  })
})
