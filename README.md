
# The biomformat Package for R

[![Travis-CI Build Status](https://travis-ci.org/joey711/biomformat.svg?branch=master)](https://travis-ci.org/joey711/biomformat)

# About 

This is an R package for interfacing with both the JSON and HDF5 versions of the [biom file format](http://biom-format.org/). This package includes basic tools for reading biom-format files, accessing and subsetting data tables from a biom object, as well as limited support for writing a biom-object back to a biom-format file.

# Design

The design of this API is intended to match the python API where appropriate, while maintaining a decidedly "R flavor" that should be familiar to R users. This includes S4 classes and methods, as well as extensions of common core functions/methods.

# Installation

To install the latest stable release of the biom package enter the following in an R session (after official release).

```S
source("http://bioconductor.org/biocLite.R")
biocLite("biom")
```

Before official release, or to install the latest development version, enter the following into an R session.

```S
install.packages("devtools") # if not already installed
devtools::install_github("biomformat", "joey711")
```

# Support

 * Please post feature or support requests and bugs at the [Issues Tracker for the biomformat package](https://github.com/joey711/biomformat/issues) on GitHub. Issues related to the format itself and not the R interface should be posted on the [issues tracker for the biom format](https://github.com/biom-format/biom-format/issues).
 
 * Note that this is a separate (but friendly!) project from the biom-format team, and the software license is different between this package and much of the rest of the biom-format software, which has switched to BSD.

 * The official release version of this package will be submitted to Bioconductor project, where [the rhdf5 package](http://bioconductor.org/packages/release/bioc/html/rhdf5.html) lives. Other than that, CRAN would be a suitable destination.

 * The original release of this package was [made available through CRAN](http://cran.r-project.org/web/packages/biom/index.html), 
 as the "biom" package, supporting version 1.x (JSON) only.
 The current plan is to let "biom" remain on CRAN in maintenance-only mode until "biomformat" is released on Bioconductor. At which point "biom" will be considered deprecated.

