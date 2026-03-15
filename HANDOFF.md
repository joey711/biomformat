# biomformat — GitHub Issues + Test Coverage + Vignette Handoff

## Repository
- GitHub: https://github.com/joey711/biomformat
- Local clone: /Users/paulmcmurdie/github/biomformat
- Active branch: devel (also fully merged to master/origin)
- Current version: 1.39.10
- R CMD check status: 0 ERRORs, 0 NOTEs, 2 pre-existing vignette WARNINGs
  (inst/doc absent — built by Bioconductor BBS, not locally; not actionable)

---

## Package Context

biomformat is a Bioconductor I/O package for the BIOM file format
(http://biom-format.org), used widely in microbiome research.
Its core jobs are:
  - read_biom() / write_biom()        JSON (BIOM v1) round-trip
  - read_hdf5_biom() / write_hdf5_biom()  HDF5 (BIOM v2) round-trip
  - biom_data()                        extract the count matrix
  - sample_metadata() / observation_metadata()  extract metadata
  - make_biom()                        construct a biom from R data
  - as.data.frame.biom() / as_tibble.biom()  tidy long-format output
  - biom_to_SummarizedExperiment()     Bioconductor SE interop

Key source files:
  R/BIOM-class.R      — make_biom(), biom_data(), generate_matrix(), etc.
  R/IO-methods.R      — read_biom(), write_biom(), read_hdf5_biom(),
                        write_hdf5_biom()
  R/SE-methods.R      — biom_to_SummarizedExperiment(), as() coercion
  R/tidy-methods.R    — as.data.frame.biom(), as_tibble.biom()
  R/allClasses.R      — S4 class definition
  R/validity-methods.R

Test files:
  tests/testthat/test-IO.R           — read/write round-trips, HDF5 routing
  tests/testthat/test-accessors.R    — biom_data(), metadata accessors,
                                       subsetting, value regression tests
  tests/testthat/test-SE.R           — SummarizedExperiment coercion
  tests/testthat/test-tidy.R         — as.data.frame.biom(), as_tibble.biom()
  tests/testthat/test-hdf5-write.R   — write_hdf5_biom() round-trips

Fixtures (inst/extdata/):
  min_dense_otu_table.biom            5×6, no metadata, dense JSON
  min_sparse_otu_table.biom           5×6, no metadata, sparse JSON
  min_sparse_otu_table_hdf5.biom      5×6, no metadata, HDF5
  rich_dense_otu_table.biom           5×6, taxonomy + sample meta, dense JSON
  rich_sparse_otu_table.biom          5×6, taxonomy + sample meta, sparse JSON
  rich_sparse_otu_table_hdf5.biom     5×6, taxonomy + sample meta, HDF5
  rich_dense_char.biom                5×6, character values, dense JSON
  rich_sparse_char.biom               5×6, character values, sparse JSON

IMPORTANT FIXTURE NOTE: min_sparse and rich_sparse are INDEPENDENT datasets
that happen to share the same count matrix (verified by test-IO.R line 36).
They differ from rich_dense at [GG_OTU_3, Sample5] and [GG_OTU_3, Sample6] —
this is intentional, not a bug. Do not "fix" the fixtures to match each other.

---

## Open GitHub Issues — Triaged

### Issue #4: write_biom does not write valid BIOM 1.0 format
https://github.com/joey711/biomformat/issues/4
STATUS: Partially addressed by existing auto_unbox=TRUE in toJSON().
REMAINING PROBLEM: The issue reporter shows that field values like id,
format, format_url, type are written as JSON arrays (["None"]) instead of
strings ("None"), and row id fields are also arrays. This means
make_biom() is building some fields as length-1 lists rather than
plain scalars, so even with auto_unbox=TRUE they serialize incorrectly.
Root cause: make_biom() builds the list with fields that are already
length-1 lists at the point they are passed to biom(). When toJSON
serializes them, auto_unbox only applies at one level.
TASK: Reproduce the issue using the reporter's code pattern:
  library(biomformat)
  x <- read_biom(system.file("extdata","min_dense_otu_table.biom",
                               package="biomformat"))
  y <- make_biom(biom_data(x))
  write_biom(y, tmp <- tempfile())
  cat(readLines(tmp, n=1))
Examine the JSON output. If id/format/rows[i].id are arrays, trace back
through make_biom() to find where scalar fields are being wrapped in lists
and fix them. Add a regression test: after write_biom(make_biom(...)),
read_biom() must parse back without error AND the resulting biom object
must be valid (validObject() returns TRUE). Also verify against the BIOM
spec that the written JSON matches the expected structure.

### Issue #6: Exporting biom with taxonomy from dada2 — metadata is array vs object
https://github.com/joey711/biomformat/issues/6
STATUS: Open. Closely related to #4.
PROBLEM: When observation_metadata is a data.frame where each row is a
taxonomy rank (columns = ranks), make_biom() serializes the per-row
metadata as a JSON ARRAY  (["Bacteria","Firmicutes",...])  rather than
a named JSON OBJECT ({"taxonomy":["Bacteria","Firmicutes",...]}).
The BIOM spec requires observation metadata to be a JSON object keyed
by field name (e.g. "taxonomy"), not a bare array.
The phyloseq import_biom() function (and the biom Python library) expect
the named-object form; the bare-array form causes the
"$ operator is invalid for atomic vectors" error the reporter saw.
TASK: Confirm the bug by calling make_biom() with a data.frame
observation_metadata that has a column named "taxonomy" containing
multi-level character vectors, then inspect the written JSON.
The fix likely requires changing how make_biom() constructs the
per-row metadata list for observation_metadata: each row must produce
a named list  list(taxonomy = c("Bacteria",...))  not an unnamed vector.
Add a fixture (inst/extdata/make_biom_taxonomy.biom or similar) and a
regression test verifying the written JSON has the correct structure.

### Issue #7: Zero-count samples being imported with wrong counts from HDF5
https://github.com/joey711/biomformat/issues/7
STATUS: Open. Reporter provides a real QIIME2 HDF5 BIOM file where a
sample with 0 total counts (NTC4) is read back with nonzero counts assigned
to two features. This is a data-corruption bug on read.
ROOT CAUSE HYPOTHESIS: In generate_matrix() (R/BIOM-class.R), when a
sample has 0 non-zero entries (indptr[j] == indptr[j+1]), the rep() call
that builds j_vec must produce 0 elements for that sample.
diff(indptr) for a zero-count sample is 0, and rep(j, 0) should be fine,
but the original sapply-based code may have had an off-by-one or
indptr-indexing error for empty columns. The new sparseMatrix-based
generate_matrix() (introduced in v1.39.8) should handle this correctly,
but it has not been validated with a fixture that contains an all-zero
column.
TASK:
1. Create a minimal HDF5 BIOM fixture with at least one all-zero sample
   column (a sample with zero counts). This can be done in R with rhdf5
   by constructing the CCS data manually for a 3×3 matrix where one
   column has no non-zero entries, or by using the attached SV_minsamp.qza
   data from the issue (unzip and use the feature-table.biom inside it).
2. Add a test that reads this fixture and verifies biom_data()[, zero_col]
   is all zeros.
3. Also test the related error path: issue mentions some files fail to read
   entirely with the old "Both attempts to read..." error. This error no
   longer exists (read_biom() now uses magic-byte routing), but confirm
   this class of file now parses successfully.

### Issue #8: write_biom — "character strings are limited to 2^31-1 bytes"
https://github.com/joey711/biomformat/issues/8
STATUS: Open. This is a jsonlite / R limitation hit when serializing a
very large BIOM object to a single JSON string.
INVESTIGATION TASK: Determine whether this is actionable. Options:
  (a) Document the limitation in ?write_biom and suggest write_hdf5_biom()
      as the alternative for large datasets (HDF5 has no such limitation).
  (b) If a streaming JSON write is feasible (write chunk-by-chunk to a
      connection), implement it.
Most likely the correct fix is (a): add a note to write_biom() docs
pointing large-dataset users to write_hdf5_biom(), and add a note to the
vignette under the "Write BIOM format" section. Document in NEWS.

### Issue #9: make_biom + write_biom for dada2 output (user support)
https://github.com/joey711/biomformat/issues/9
STATUS: This is a user-support issue, not a distinct bug from #4 and #6.
Once #4 and #6 are fixed the user's workflow should work. No separate
code change needed; close with reference to #4/#6 fixes once those are
done, and add a vignette example showing the dada2 → biom → write pattern
(count table rows=ASVs, cols=samples; taxonomy as observation_metadata
data.frame with a "taxonomy" column containing character vectors of ranks).

### Issue #15: Loading group metadata (BIOM v2 feature request)
https://github.com/joey711/biomformat/issues/15
STATUS: Open feature request. BIOM format 2.1.0 spec defines
"group metadata" — additional per-sample or per-observation metadata
stored in HDF5 group attributes rather than the standard datasets.
This is a new feature, not a bug.
TASK: Investigate scope:
  1. Read the BIOM 2.1.0 spec section on group metadata:
     http://biom-format.org/documentation/format_versions/biom-2.1.html
  2. Check whether the existing HDF5 fixture (rich_sparse_otu_table_hdf5.biom)
     contains any group metadata, or if a new fixture would be needed.
  3. If the feature is straightforward (a few additional h5read calls
     in read_hdf5_biom() and h5write calls in write_hdf5_biom()), implement
     it. If it requires significant new API surface, file a scoped plan.

---

## Test Coverage Gaps — New Fixtures Needed

Current fixtures are ALL 5×6 matrices with the same OTU/sample IDs.
Missing edge cases:

1. ALL-ZERO SAMPLE COLUMN (critical — directly related to issue #7):
   A sparse or HDF5 BIOM where at least one sample has count=0 for all
   features. Tests: biom_data()[, zero_col] is all zeros; make_biom()
   roundtrip preserves the zero column; as.data.frame() includes that
   sample with count=0 for all features.

2. SINGLE-SAMPLE BIOM (1×N matrix):
   A BIOM with exactly one sample. Tests: ncol(x)==1, biom_data() returns
   a 1-column matrix (not a vector), write_biom/read_biom round-trip.

3. SINGLE-FEATURE BIOM (M×1 matrix):
   A BIOM with exactly one feature/OTU. Tests: nrow(x)==1, biom_data()
   returns a 1-row matrix (not a vector), round-trip.

4. LARGE TAXONOMY make_biom ROUND-TRIP (directly tests #6):
   A small fixture (3×3) created via make_biom() with observation_metadata
   that has a "taxonomy" column containing 7-level character vectors, then
   written with write_biom() and read back. The read-back metadata must
   match the input.

5. HDF5 BIOM WITH ALL-ZERO SAMPLE (for issue #7):
   See issue #7 above. Minimal HDF5 fixture, 3 obs × 3 samples, middle
   sample all zero. Can be created programmatically with rhdf5 in a
   helper script rather than committed as a binary if preferred.

6. FLOAT-VALUED BIOM:
   A BIOM with matrix_element_type="float" and non-integer values.
   Currently all numeric fixtures are "int". Tests: biom_data() returns
   numeric (not integer) values; write_biom/read_biom round-trip preserves
   float precision within jsonlite defaults.

---

## Vignette Gaps

The current vignette (vignettes/biomformat.Rmd) covers:
  - Read/write JSON BIOM (read_biom, write_biom)
  - biom_data(), observation_metadata(), sample_metadata() accessors
  - HDF5 read and write (write_hdf5_biom, read_biom)
  - Tidy long-format (as.data.frame.biom, as_tibble.biom)
  - purrr-style summarisation
  - SummarizedExperiment coercion

Missing vignette content:

1. CONSTRUCTING A BIOM FROM R DATA (make_biom + write_biom workflow):
   This is the most common user question (issues #4, #6, #9). Add a
   section showing:
     mat <- matrix(sample(0:100, 15), nrow=3, ncol=5,
                   dimnames=list(paste0("OTU",1:3), paste0("Samp",1:5)))
     tax <- data.frame(taxonomy=I(lapply(1:3, function(i)
               paste0("level",1:7,"_",i))), row.names=rownames(mat))
     x <- make_biom(data=mat, observation_metadata=tax,
                    matrix_element_type="int")
     tmp <- tempfile()
     write_biom(x, tmp)
     y <- read_biom(tmp)
     identical(biom_data(x), biom_data(y))
   This directly addresses issues #4, #6, #9 for users landing on the
   vignette. Include a note about write_hdf5_biom() for large datasets
   (addressing issue #8).

2. SUBSETTING biom_data() BY NAME:
   The biom_data() function accepts character row/column names for
   subsetting, but this is not shown in the vignette or examples.
   Add a short code block demonstrating named subsetting.

---

## Workflow for This Session

Recommended order of operations:

1. REPRODUCE AND FIX ISSUE #4 (write_biom invalid JSON).
   This is the highest-impact fix — it affects every user who calls
   make_biom() and the problem has been open for years.
   - Use the reproduction recipe above.
   - Trace the bug through make_biom() to find scalar-vs-list wrapping.
   - Fix and add regression test.
   - Bump version to v1.39.11.

2. FIX ISSUE #6 (taxonomy metadata serialized as array not object).
   Likely the same root cause as #4 or an adjacent one.
   - Fix the observation_metadata named-list construction in make_biom().
   - Add fixture and regression test.
   - Include in v1.39.11 or bump to v1.39.12.

3. CREATE MISSING FIXTURES AND TESTS (all-zero column, 1×N, N×1, float).
   - Write a helper script (inst/extdata/create_fixtures.R or similar)
     that generates the new fixtures programmatically and documents their
     contents.
   - Add corresponding tests.

4. INVESTIGATE ISSUE #7 (zero-count sample data corruption).
   - Validate that the new generate_matrix() (v1.39.8) correctly handles
     all-zero columns using the new fixture from step 3.
   - If a bug is found, fix it.

5. ADDRESS ISSUE #8 (write_biom size limit).
   - Add documentation note and vignette callout pointing to
     write_hdf5_biom() for large datasets.

6. ADD VIGNETTE SECTIONS (make_biom workflow, named subsetting).

7. INVESTIGATE ISSUE #15 (group metadata) — scope and implement or file
   a bounded plan.

8. Throughout: after each substantive fix, run R CMD check and confirm
   0 ERRORs / 0 NOTEs before committing.

---

## Key Invariants to Preserve

- Do NOT change the count values in any existing fixture file.
  min_sparse != rich_dense at [GG_OTU_3, Sample5/Sample6] — intentional.
- Matrix orientation in sparse JSON: triplets are [obs_idx, sample_idx, value]
  (0-based), matching shape=[n_obs, n_samples].
- Matrix orientation in HDF5: sample/matrix uses CCS (outer=samples),
  observation/matrix uses CSR (outer=obs). generate_matrix() reads
  sample/matrix only. Both must be written by write_hdf5_biom().
- The return type of generate_matrix() is a list-of-named-vectors (one
  per observation row). Do not change this — it feeds directly into the
  biom() constructor.
- R CMD check must pass with 0 ERRORs and 0 NOTEs before any merge.
  The 2 vignette WARNINGs (inst/doc absent) are pre-existing and expected.

---

## Commit Discipline

- One logical fix per commit.
- Commit message format:
    Fix #N: short description (vX.Y.Z)
    
    Longer explanation of root cause and fix.
    
    R CMD check: 0 ERRORs, 0 NOTEs, 2 pre-existing vignette WARNINGs.
- After each fix, bump the patch version in DESCRIPTION and prepend a
  NEWS entry.
- After completing all work, fast-forward origin/master to match devel:
    git checkout master
    git merge --ff-only devel
    git push origin master devel
