################################################################################
# Use testthat to test file import and resulting class (and values)
################################################################################
library("biomformat"); library("testthat")
library("Matrix")

packageVersion("biomformat")

mat1 = matrix(c(1,0,2,0,0,3,4,5,6), nrow = 3, byrow = TRUE)
rownames(mat1) <- paste0("Taxa_", 1:nrow(mat1))
colnames(mat1) <- paste0("Sample_", 1:ncol(mat1))
mat2 = as(mat1, "CsparseMatrix")

test_that("Test biom-format v2.1 sparse CSC accessors", {
  # index pointer
  # 0 2 3 6
  expect_equal(
    c(0, 2, 3, 6),
    sparse_indptr(mat1)
  )
  expect_equal(
    c(0, 2, 3, 6),
    sparse_indptr(mat2)
  )
  # Indices
  # 0 2 2 0 1 2
  expect_equal(
    c(0, 2, 2, 0, 1, 2),
    sparse_indices(mat1)
  )
  expect_equal(
    c(0, 2, 2, 0, 1, 2),
    sparse_indices(mat2)
  )
  # data
  # 1 4 5 2 3 6
  expect_equal(
    c(1, 4, 5, 2, 3, 6),
    sparse_data(mat1)
  )
  expect_equal(
    c(1, 4, 5, 2, 3, 6),
    sparse_data(mat2)
  )
})

