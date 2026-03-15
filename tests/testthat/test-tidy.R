################################################################################
# Milestone 3: tidy long-format interoperability tests
################################################################################
library("biomformat")
library("testthat")

# ── Fixtures ──────────────────────────────────────────────────────────────────
rich_dense_file <- system.file("extdata", "rich_dense_otu_table.biom",
                               package = "biomformat")
min_dense_file  <- system.file("extdata", "min_dense_otu_table.biom",
                               package = "biomformat")
rich_sparse_file <- system.file("extdata", "rich_sparse_otu_table.biom",
                                package = "biomformat")

x_rich   <- read_biom(rich_dense_file)   # 5 obs x 6 samples, full metadata
x_min    <- read_biom(min_dense_file)    # 5 obs x 6 samples, no metadata
x_sparse <- read_biom(rich_sparse_file)  # 5 obs x 6 samples, full metadata, sparse

# ── Tests ─────────────────────────────────────────────────────────────────────

test_that("as.data.frame.biom returns correct row count", {
  # 5 features x 6 samples = 30 rows in long form
  df <- as.data.frame(x_rich)
  expect_equal(nrow(df), 30L)
})

test_that("as.data.frame.biom has correct base column names", {
  df <- as.data.frame(x_rich)
  expect_true("feature_id" %in% colnames(df))
  expect_true("sample_id"  %in% colnames(df))
  expect_true("count"      %in% colnames(df))
})

test_that("as.data.frame.biom returns a plain data.frame", {
  df <- as.data.frame(x_rich)
  expect_true(is.data.frame(df))
})

test_that("count values match biom_data() matrix values", {
  df  <- as.data.frame(x_rich)
  mat <- as(biom_data(x_rich), "matrix")
  # Spot-check: GG_OTU_3 / Sample4 should equal mat[3,4] = 4
  cell <- df[df$feature_id == "GG_OTU_3" & df$sample_id == "Sample4", "count"]
  expect_equal(as.numeric(cell), as.numeric(mat["GG_OTU_3", "Sample4"]))
  # Zero cell: GG_OTU_1 / Sample1
  cell0 <- df[df$feature_id == "GG_OTU_1" & df$sample_id == "Sample1", "count"]
  expect_equal(as.numeric(cell0), as.numeric(mat["GG_OTU_1", "Sample1"]))
  # All 30 values collectively match the matrix (order-insensitive)
  mat_vals  <- sort(as.numeric(mat))
  long_vals <- sort(as.numeric(df$count))
  expect_equal(long_vals, mat_vals)
})

test_that("sample metadata columns are present when metadata exists", {
  df <- as.data.frame(x_rich)
  # rich_dense has 4 sample metadata columns
  expect_true("BarcodeSequence"    %in% colnames(df))
  expect_true("LinkerPrimerSequence" %in% colnames(df))
  expect_true("BODY_SITE"          %in% colnames(df))
  expect_true("Description"        %in% colnames(df))
  # Total columns: 3 base + 4 sample + 7 obs = 14
  expect_equal(ncol(df), 14L)
})

test_that("observation metadata columns are present when metadata exists", {
  df <- as.data.frame(x_rich)
  # rich_dense has 7 observation (taxonomy) metadata columns
  expect_true("taxonomy1" %in% colnames(df))
  expect_true("taxonomy7" %in% colnames(df))
})

test_that("NULL metadata works without error (min_dense has no metadata)", {
  # Must not error, must still produce 30-row data.frame with 3 columns
  df <- expect_no_error(as.data.frame(x_min))
  expect_equal(nrow(df), 30L)
  expect_equal(ncol(df), 3L)
  expect_equal(sort(colnames(df)), c("count", "feature_id", "sample_id"))
})

test_that("sparse biom long-format count values match its own biom_data()", {
  # The dense and sparse fixture files are not identical matrices (two cells
  # differ: GG_OTU_3/Sample5 and GG_OTU_3/Sample6 are swapped), so we
  # validate each fixture independently against its own biom_data() output.
  df_sparse <- as.data.frame(x_sparse)
  mat_s     <- as(biom_data(x_sparse), "matrix")
  expect_equal(nrow(df_sparse), 30L)
  # Build key -> count lookup and compare against the matrix
  key_s  <- paste(df_sparse$feature_id, df_sparse$sample_id, sep = "|")
  vals_s <- setNames(as.numeric(df_sparse$count), key_s)
  for (feat in rownames(mat_s)) {
    for (samp in colnames(mat_s)) {
      k <- paste(feat, samp, sep = "|")
      expect_equal(unname(vals_s[k]), as.numeric(mat_s[feat, samp]),
                   label = paste0("sparse[", feat, ",", samp, "]"))
    }
  }
})

test_that("as_tibble.biom works when tibble is installed", {
  skip_if_not_installed("tibble")
  # Use fully-qualified call so the generic is always found regardless of
  # whether tibble is attached; as_tibble.biom() also works when called directly.
  tbl <- tibble::as_tibble(x_rich)
  expect_true(inherits(tbl, "tbl_df"))
  expect_equal(nrow(tbl), 30L)
  expect_true("feature_id" %in% colnames(tbl))
  expect_true("count"      %in% colnames(tbl))
})
