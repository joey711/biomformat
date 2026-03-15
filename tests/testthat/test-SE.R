################################################################################
# Milestone 2: SummarizedExperiment interoperability tests
################################################################################
library("biomformat")
library("testthat")

# Skip entire file if SummarizedExperiment is not installed
skip_if_not_installed("SummarizedExperiment")
skip_if_not_installed("S4Vectors")

suppressPackageStartupMessages({
  library("SummarizedExperiment")
  library("S4Vectors")
})

# ── Fixtures ──────────────────────────────────────────────────────────────────
rich_dense_file <- system.file("extdata", "rich_dense_otu_table.biom",
                               package = "biomformat")
min_dense_file  <- system.file("extdata", "min_dense_otu_table.biom",
                               package = "biomformat")

x_rich <- read_biom(rich_dense_file)   # 5 obs x 6 samples, has both metadata tables
x_min  <- read_biom(min_dense_file)    # 5 obs x 6 samples, no metadata at all

# ── Tests ─────────────────────────────────────────────────────────────────────

test_that("biom_to_SummarizedExperiment() returns a SummarizedExperiment", {
  se <- biom_to_SummarizedExperiment(x_rich)
  expect_true(is(se, "SummarizedExperiment"))
})

test_that("assay 'counts' matches biom_data() exactly", {
  se  <- biom_to_SummarizedExperiment(x_rich)
  mat <- biom_data(x_rich)
  # Compare as plain matrices — storage class may differ (dgeMatrix vs dgCMatrix)
  expect_equal(as.matrix(assay(se, "counts")), as.matrix(mat))
})

test_that("colData rows match sample_metadata() rows", {
  se  <- biom_to_SummarizedExperiment(x_rich)
  smd <- sample_metadata(x_rich)
  # Row names of colData must match sample names
  expect_equal(rownames(colData(se)), rownames(smd))
  # Column names must match too
  expect_equal(colnames(colData(se)), colnames(smd))
  # Spot-check a value
  expect_equal(as.character(colData(se)[["BODY_SITE"]][2]),
               smd[2, "BODY_SITE"])
})

test_that("rowData rows match observation_metadata() rows", {
  se  <- biom_to_SummarizedExperiment(x_rich)
  omd <- observation_metadata(x_rich)
  expect_equal(rownames(rowData(se)), rownames(omd))
  expect_equal(colnames(rowData(se)), colnames(omd))
  # Spot-check a taxonomy value
  expect_equal(as.character(rowData(se)[[1]][1]),
               as.character(omd[1, 1]))
})

test_that("as(x, 'SummarizedExperiment') coercion gives the same result", {
  se1 <- biom_to_SummarizedExperiment(x_rich)
  se2 <- as(x_rich, "SummarizedExperiment")
  expect_equal(as.matrix(assay(se1, "counts")),
               as.matrix(assay(se2, "counts")))
  expect_equal(rownames(colData(se1)), rownames(colData(se2)))
  expect_equal(rownames(rowData(se1)), rownames(rowData(se2)))
})

test_that("NULL metadata produces empty DataFrame in colData and rowData", {
  se <- biom_to_SummarizedExperiment(x_min)
  # No metadata in min_dense — both should be empty DataFrames with correct names
  expect_equal(nrow(colData(se)), ncol(biom_data(x_min)))
  expect_equal(ncol(colData(se)), 0L)
  expect_equal(nrow(rowData(se)), nrow(biom_data(x_min)))
  expect_equal(ncol(rowData(se)), 0L)
  # Row names must still be set from the matrix dim names
  expect_equal(rownames(colData(se)), colnames(biom_data(x_min)))
  expect_equal(rownames(rowData(se)), rownames(biom_data(x_min)))
})

test_that("SE dimensions match original biom object", {
  se <- biom_to_SummarizedExperiment(x_rich)
  # biom nrow()/ncol() return *named* integers; SE returns plain integers.
  # Compare the numeric value only via as.integer().
  expect_equal(as.integer(nrow(se)), as.integer(nrow(x_rich)))
  expect_equal(as.integer(ncol(se)), as.integer(ncol(x_rich)))
  expect_equal(rownames(se), rownames(biom_data(x_rich)))
  expect_equal(colnames(se), colnames(biom_data(x_rich)))
})
