# Provide the data of gt_radmix() for a class (for package developers)

[`gt_radmix()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/gt_radmix.md)
reads the data of `x` through `gt_radmix_data()`, a generic with a
method for a
[`tidypopgen::gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.html).
Packages that store their data differently can provide methods for their
own classes: this is how
[`gt_radmix()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/gt_radmix.md)
fits genotype likelihoods, without tidygenclust depending on the
packages that store them. Users do not need to call it.

## Usage

``` r
gt_radmix_data(x, ...)

# Default S3 method
gt_radmix_data(x, ...)

# S3 method for class 'gen_tbl'
gt_radmix_data(x, use_imputed = FALSE, ...)
```

## Arguments

- x:

  the object given to
  [`gt_radmix()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/gt_radmix.md)

- ...:

  further arguments for methods

- use_imputed:

  boolean; for a `gen_tibble`, whether imputed genotypes should be used,
  as in
  [`gt_radmix()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/gt_radmix.md)

## Value

a list (see the details)

## Details

A method for genotype likelihoods returns a list with the elements:

- `kind`: `"gl_source"`;

- `version`: `1L`, the version of the interface;

- `fill`: an external pointer (made with `R_MakeExternalPtrFn()`) to a C
  function that fills the likelihoods of a block of individuals at one
  locus, tagged with the symbol `tidygenclust_gl_fill`;

- `ctx`: an external pointer to the data that the function reads (its
  context);

- `n_indiv`, `n_loci`: the numbers of individuals and loci, in the order
  of the rows of `x` and of
  [`tidypopgen::show_loci()`](https://evolecolgroup.github.io/tidypopgen/reference/show_loci.html)
  (which label the rows of Q and of P);

and any other element that must stay alive while
[`gt_radmix()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/gt_radmix.md)
uses the source (e.g. the file-backed matrix that the context points
into). The C function is specified in the header
`include/tidygenclust/radmix_gl_source.h` of the installed package,
which a package can copy (or use with `LinkingTo: tidygenclust`): it
gives the likelihoods of 1 and 2 copies of `allele2` relative to that of
0 copies, and whether each individual has data, and it is called from
several threads, so it must not call R. The file
`src/radmix_gl_source.cpp` of tidygenclust's sources is an example, used
by its tests. `P` then gives the frequency of `allele2`. Starts from the
SVD and imputed genotypes are not available for genotype likelihoods:
[`gt_radmix()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/gt_radmix.md)
checks them itself, and calls a method with `x` alone (`use_imputed` is
only given to the method for a `gen_tibble`, when it is `TRUE`).

The list for a `gen_tibble` is internal to tidygenclust, and may change.

## See also

[`gt_radmix()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/gt_radmix.md),
whose details describe the model of genotype likelihoods
