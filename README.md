# The biomformat Package for R

[![Bioconductor release](https://img.shields.io/badge/Bioconductor-release-blue)](https://bioconductor.org/packages/release/bioc/html/biomformat.html)
[![Bioconductor devel](https://img.shields.io/badge/Bioconductor-devel-green)](https://bioconductor.org/packages/devel/bioc/html/biomformat.html)

## About

This is an R package for interfacing with both the JSON (BIOM v1) and HDF5
(BIOM v2) versions of the [BIOM file format](http://biom-format.org/).
It includes tools for reading BIOM files, accessing and subsetting the data
matrix and associated metadata from a biom object, as well as support for
writing a biom object back to a BIOM-format file.

## Installation

Install the current Bioconductor release version:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomformat")
```

Install the latest development version directly from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("joey711/biomformat")
```

## Quick start

```r
library(biomformat)

# Read a BIOM file (JSON v1 or HDF5 v2 — format is detected automatically)
biom_file <- system.file("extdata", "rich_sparse_otu_table.biom",
                         package = "biomformat")
x <- read_biom(biom_file)

# Inspect the object
x
biom_shape(x)   # rows (observations) x columns (samples)

# Extract the count matrix (returns a Matrix object)
biom_data(x)

# Extract sample and observation metadata
sample_metadata(x)
observation_metadata(x)

# Build a biom object from R data and write it to file
y <- make_biom(data = my_count_matrix,
               sample_metadata = my_sample_df,
               observation_metadata = my_tax_df)
write_biom(y, "output.biom")
```

## Recent changes (v1.39.2 – v1.39.4)

These versions rescue the package from Bioconductor deprecation and fix
several long-standing data-integrity bugs. See [NEWS](NEWS) for full details
and GitHub issue/PR cross-references.

| Version | Summary |
|---------|---------|
| **1.39.4** | `biom_data()`: `drop = FALSE` on dense and sparse paths — single-row/col subsets now always return a matrix, never a vector. Closes [PR #12](https://github.com/joey711/biomformat/pull/12). |
| **1.39.3** | `read_biom()`: deterministic HDF5 magic-bytes router replaces fragile JSON-first fallback, eliminating the `lexical error: invalid char` warning. Closes [Issue #14](https://github.com/joey711/biomformat/issues/14), supersedes [PR #16](https://github.com/joey711/biomformat/pull/16). `read_hdf5_biom()`: graceful `requireNamespace()` + `tryCatch()` guards replace fatal C-level aborts. Bumped `Matrix >= 1.7-0`. |
| **1.39.2** | Fix fatal `R CMD check` ERROR from deprecated `testthat` 2.x helpers (`expect_that`/`is_true()`). Addresses [Issue #17](https://github.com/joey711/biomformat/issues/17). |

## Design

The API is intended to match the Python biom-format API where appropriate,
while maintaining a decidedly "R flavour" — S4 classes and methods, and
extensions of common base-R generics (`nrow`, `ncol`, `rownames`,
`colnames`, `show`).

## Support

- Please post bug reports and feature requests on the
  [Issues Tracker](https://github.com/joey711/biomformat/issues).
- Issues related to the BIOM format specification (rather than the R
  interface) should go to the
  [biom-format team's tracker](https://github.com/biocore/biom-format/issues).
- For Bioconductor-specific questions, use the
  [Bioconductor support site](https://support.bioconductor.org/) or the
  [bioc-devel mailing list](mailto:bioc-devel@r-project.org).
- This package is available on
  [Bioconductor](https://bioconductor.org/packages/biomformat/), where
  the [rhdf5](http://bioconductor.org/packages/rhdf5/) dependency also lives.
