## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Notes

* `gt_fastmixture()` and `gt_adamixture()` run Python programs through
  `reticulate`, and `tgc_tools_install()` installs them in a conda
  environment only when the user calls it. The tests and examples that need
  these programs are skipped when they are not installed.
