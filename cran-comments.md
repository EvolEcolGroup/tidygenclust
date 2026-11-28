## Test environments

* local: Ubuntu 24.04, R 4.6.2 (release)
* GitHub Actions: Ubuntu (R-devel, release, oldrel-1), macOS and Windows (release)
* win-builder: R-devel (Windows)
* R-hub: linux, macos-arm64, windows and m1-san (ASan/UBSan), all R-devel

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* Possibly misspelled words in DESCRIPTION: 'natively' and 'pseudohaploid'
  are correct ('pseudohaploid' is a term of population genetics: one allele
  sampled at each locus of low-coverage genomes).

## Notes

* `gt_fastmixture()` and `gt_adamixture()` run optional Python programs
  through 'reticulate', which `tgc_tools_install()` installs in a conda
  environment only when the user calls it. Their tests are skipped on CRAN,
  and their examples are in `\donttest{}` and only run when the programs are
  installed. The examples of the installers are in `\dontrun{}`, as they
  install software.
