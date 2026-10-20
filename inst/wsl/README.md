# WSL dependency locks

Each environment has an explicit Conda lock (package URLs, builds and SHA-256
hashes) and a pip lock (versions and compatible distribution SHA-256 hashes).
These are linux-64 locks for WSL Ubuntu, not Windows-native or macOS locks.
Tool source commits remain in the R installer's existing arguments; the shell
installer pins micromamba separately.

Fresh installations and `reset = TRUE` use these locks. Existing environments
are reused, not automatically rebuilt. Adding or changing a tool installs its
locked dependencies as needed. Source builds use those dependencies directly
(`--no-build-isolation --no-deps`), followed by `pip check`; custom source
revisions must remain compatible. Pip's cache is disabled during installation
to avoid reusing locally built wheels made with different build dependencies.
ADAMIXTURE is still optional, although its
small CPU dependency set is included in the shared core dependency lock.

To update the locks, first test the new dependencies in disposable environments.
Run `data-raw/wsl_locks.py` with each environment's Python, passing the
micromamba executable path. It prints JSON with `conda` and `pip` contents for
the corresponding lockfiles. Review the changes, including the source-build
backend versions specified in that script, and retain the two header comments.
The exporter retains Conda-forge packages, includes Python source-build
dependencies, and selects wheel hashes compatible with
Ubuntu 22.04 or newer, falling back to a source archive where necessary.
Rebuild from the updated locks and test all three tools, optional ADAMIXTURE
installation, repeated installation and explicit reset before accepting them.

These locks freeze the analysis environments, not Ubuntu's compiler/system
packages, Windows/R packages or hardware. They do not promise bitwise-identical
results across machines. Keep data, source versions, seeds and thread counts
with an analysis, and record `sessionInfo()` and `tgc_tools_version()`.

## Continuous integration

`.github/workflows/windows-wsl.yaml` adds a Windows R (release) check alongside
the unchanged Linux/macOS workflow. It runs on pull requests, pushes to `main`
and manual dispatch. Each run installs Ubuntu 24.04 in WSL2, creates an ordinary
Linux user, and calls `tgc_tools_install(adamixture = TRUE)` twice from Windows R.
The analysis environments are not cached, so this checks a fresh installation
from the locks as well as reuse. Only the distribution installer and R
dependencies use the setup actions' usual caches.

`TIDYGENCLUST_WSL_CI=true` makes missing tools and mismatched pip dependencies
fail the tests, rather than silently skip analyses. The usual `R CMD check`
builds both vignettes and tests all three tools, including Windows paths with
spaces and non-ASCII characters. The standard `r-lib/actions` check step
preserves check logs on failure; installer output is in the job log.

The `windows-2025` image includes WSL2, but GitHub does not officially guarantee
[nested virtualisation](https://docs.github.com/en/actions/concepts/runners/github-hosted-runners).
Confirm the first hosted run before treating this job as a required check.
A WSL setup failure must remain visible, not fall back to WSL1 or skip the job.
