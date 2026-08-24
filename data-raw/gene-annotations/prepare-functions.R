# Source the package-owned, dependency-light preparation helpers. Keeping the
# implementation under R/ lets package tests exercise the same code while this
# data-raw wrapper keeps the maintainer workflow runnable with Rscript.

script_dir <- dirname(normalizePath(sys.frame(1)$ofile, mustWork = TRUE))
source(file.path(script_dir, "..", "..", "R", "gene-annotation-prep.R"))
