# tidygenclust 0.1.1

* Updated versions of fastmixture and clumppling
* Added experimental Windows support for Fastmixture, ADAMIXTURE and Clumppling
  through WSL2 while keeping R and RStudio on Windows
* Added a Windows setup vignette and Windows R/WSL2 continuous integration
* Added persistent WSL distribution selection with `TIDYGENCLUST_WSL_DISTRO`
  and checked prerequisites and micromamba downloads during installation
* Locked WSL Conda builds and Python runtime/build dependencies with download
  checksums, while preserving existing environments unless explicitly reset
* Fixed the pinned ADAMIXTURE 1.8.5 commit

# tidygenclust 0.1.0

* Initial public version
