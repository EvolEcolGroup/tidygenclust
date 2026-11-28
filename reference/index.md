# Package index

## ADMIXTURE and clumppling

Native estimation of ancestry proportions with the ADMIXTURE method, and
alignment of runs across K with clumppling, with their plotting
functions.

- [`gt_radmix()`](https://evolecolgroup.github.io/tidygenclust/reference/gt_radmix.md)
  : Estimate ancestry proportions with the ADMIXTURE method

- [`gt_clumppling()`](https://evolecolgroup.github.io/tidygenclust/reference/gt_clumppling.md)
  : Run clumppling

- [`subset_gt_clumppling()`](https://evolecolgroup.github.io/tidygenclust/reference/subset_gt_clumppling.md)
  :

  Subset a `gt_clumppling` object

- [`tidy(`*`<gt_clumppling>`*`)`](https://evolecolgroup.github.io/tidygenclust/reference/tidy_gt_clumppling.md)
  :

  Tidy a `gt_clumppling` object

- [`autoplot(`*`<gt_clumppling>`*`)`](https://evolecolgroup.github.io/tidygenclust/reference/autoplot_gt_clumppling.md)
  : autoplot for clumppling objects

## External tools

Optional tools, installed in the conda environment `ctidygenclust`
(Linux only).

### Managing the tools

- [`tgc_tools_install()`](https://evolecolgroup.github.io/tidygenclust/reference/tgc_tools_install.md)
  :

  Install optional tools for `tidygenclust`

- [`tgc_tools_version()`](https://evolecolgroup.github.io/tidygenclust/reference/tgc_tools_version.md)
  :

  Print version of the optional tools installed by `tidygenclust`

### fastmixture

- [`tgc_fastmixture_install()`](https://evolecolgroup.github.io/tidygenclust/reference/tgc_fastmixture_install.md)
  :

  Add `fastmixture` to the conda environment of `tidygenclust`

- [`gt_fastmixture()`](https://evolecolgroup.github.io/tidygenclust/reference/gt_fastmixture.md)
  : fastmixture algorithm for population genetics clustering

### ADAMIXTURE

- [`tgc_adamixture_install()`](https://evolecolgroup.github.io/tidygenclust/reference/tgc_adamixture_install.md)
  :

  Add `ADAMIXTURE` to the conda environment of `tidygenclust`

- [`gt_adamixture()`](https://evolecolgroup.github.io/tidygenclust/reference/gt_adamixture.md)
  : ADAMIXTURE algorithm for population genetics clustering

### ADMIXTURE

The ADMIXTURE binary is run by
[`tidypopgen::gt_admixture()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_admixture.html).

- [`tgc_admixture_install()`](https://evolecolgroup.github.io/tidygenclust/reference/tgc_admixture_install.md)
  :

  Add `ADMIXTURE` to the conda environment of `tidygenclust`

## Example data

- [`capeverde_pops`](https://evolecolgroup.github.io/tidygenclust/reference/capeverde_pops.md)
  : Population membership for the Capeverde dataset
