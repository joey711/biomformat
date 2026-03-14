################################################################################
# Use testthat to test file import and resulting class (and values)
################################################################################
library("biomformat"); library("testthat")
packageVersion("biomformat")
# # # # TESTS!
min_dense_file   = system.file("extdata", "min_dense_otu_table.biom", package = "biomformat")
min_sparse_file  = system.file("extdata", "min_sparse_otu_table.biom", package = "biomformat")
rich_dense_file  = system.file("extdata", "rich_dense_otu_table.biom", package = "biomformat")
rich_sparse_file = system.file("extdata", "rich_sparse_otu_table.biom", package = "biomformat")
min_dense_file   = system.file("extdata", "min_dense_otu_table.biom", package = "biomformat")
rich_dense_char  = system.file("extdata", "rich_dense_char.biom", package = "biomformat")
rich_sparse_char  = system.file("extdata", "rich_sparse_char.biom", package = "biomformat")
# Test read biom
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
