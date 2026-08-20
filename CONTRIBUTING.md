# Contributing

This repository accompanies a single methods paper and is not actively
maintained software in the usual sense — think of it as a citable,
reproducible snapshot rather than an ongoing project looking for
contributors. That said, questions, bug reports, and small fixes are welcome.

## Questions or issues

Please open a [GitHub issue](../../issues) rather than emailing directly, so
answers are visible to other users hitting the same thing.

## If you'd like to submit a change

- **Coding conventions**: paths via `here::here()` (never hardcoded absolute
  paths, never `setwd()`); all random operations fixed at `seed = 42`;
  functions documented with roxygen-style comments; scripts runnable from
  the command line via `Rscript`, not just interactively; tidyverse style.
- **Tests**: if you're changing anything in `R/`, please add or update a
  corresponding test in `tests/testthat/`. Run the full suite with
  `Rscript tests/testthat.R` before opening a pull request.
- **Dependencies**: if your change adds a new package dependency, run
  `renv::snapshot()` afterwards and include the updated `renv.lock` in your
  pull request.

## Scope

Bug fixes and clarity improvements to the existing feature extraction and
data-pathway scripts are the best fit. Feature additions that go beyond
what the accompanying paper describes may be politely declined, or asked to
be pursued as a fork, since this repo is meant to stay a faithful
implementation of the published method.
