test_that("WSL distribution selection respects the environment and R option", {
  old_env <- Sys.getenv("TIDYGENCLUST_WSL_DISTRO", unset = NA_character_)
  old_options <- options(tidygenclust.wsl_distro = NULL)
  tryCatch({
    Sys.unsetenv("TIDYGENCLUST_WSL_DISTRO")
    expect_identical(.tgc_wsl_distro(), "Ubuntu")
    Sys.setenv(TIDYGENCLUST_WSL_DISTRO = "Ubuntu-24.04")
    expect_identical(.tgc_wsl_distro(), "Ubuntu-24.04")
    options(tidygenclust.wsl_distro = "Ubuntu-22.04")
    expect_identical(.tgc_wsl_distro(), "Ubuntu-22.04")
  }, finally = {
    options(old_options)
    if (is.na(old_env)) {
      Sys.unsetenv("TIDYGENCLUST_WSL_DISTRO")
    } else {
      Sys.setenv(TIDYGENCLUST_WSL_DISTRO = old_env)
    }
  })
})

test_that("WSL dependency locks are portable and checksum-pinned", {
  for (envname in c("ctidygenclust", "cclumppling")) {
    for (type in c("conda", "pip")) {
      path <- system.file("wsl", paste0(envname, "-", type, ".lock"),
        package = "tidygenclust", mustWork = TRUE
      )
      entries <- readLines(path)
      entries <- entries[nzchar(entries) & !startsWith(entries, "#")]
      expect_true(length(entries) > 1L)
      expect_identical(anyDuplicated(entries), 0L)
      if (type == "conda") {
        expect_identical(entries[1], "@EXPLICIT")
        expect_true(all(grepl(
          paste0("^https://conda[.]anaconda[.]org/conda-forge/",
            "(linux-64|noarch)/[^/]+#[a-f0-9]{64}$"
          ),
          entries[-1]
        )))
      } else {
        expect_true(all(grepl(
          "^[a-z0-9-]+==[0-9A-Za-z.+]+( --hash=sha256:[a-f0-9]{64})+$",
          entries
        )))
      }
    }
  }
})

test_that("Windows CI has all required WSL tools", {
  skip_if(Sys.getenv("TIDYGENCLUST_WSL_CI") != "true")

  # A missing backend must fail CI, rather than skip every analysis test.
  expect_identical(.Platform$OS.type, "windows")
  expect_identical(.tgc_backend(), "wsl")
  # Ensure the path tests below exercise non-ASCII names too.
  expect_true(l10n_info()[["UTF-8"]])
  tools <- c("fastmixture", "adamixture", "clumppling")
  for (tool in tools) {
    expect_true(.tgc_tool_available(tool), info = tool)
  }
  expect_setequal(names(tgc_tools_version()), tools)
  expect_false("gt_admixture" %in% getNamespaceExports("tidygenclust"))
})

test_that("fresh Windows CI environments match the pip locks", {
  skip_if(Sys.getenv("TIDYGENCLUST_WSL_CI") != "true")

  # CI starts with a fresh WSL distribution; ordinary user environments may
  # have been retained from an earlier installation and are not checked here.
  for (envname in c("ctidygenclust", "cclumppling")) {
    path <- system.file("wsl", paste0(envname, "-pip.lock"),
      package = "tidygenclust", mustWork = TRUE
    )
    entries <- readLines(path)
    entries <- sub(" .*", "", entries[grepl("^[a-z]", entries)])
    requirements <- do.call(rbind, strsplit(entries, "==", fixed = TRUE))
    code <- paste(
      "import importlib.metadata as m, json, sys;",
      "print(json.dumps([m.version(n) for n in sys.argv[1:]]))"
    )
    versions <- .tgc_wsl_run_env(envname, "python",
      c("-c", code, requirements[, 1]), stdout = TRUE, stderr = TRUE
    )
    expect_null(attr(versions, "status"), info = envname)
    expect_identical(jsonlite::fromJSON(versions),
      unname(requirements[, 2]), info = envname
    )
  }
})

test_that("binary matrix bridge preserves row-major matrices", {
  values <- matrix(seq_len(12), nrow = 3, byrow = TRUE)
  path <- tempfile(fileext = ".bin")
  connection <- file(path, open = "wb")
  writeBin(
    as.double(t(values)),
    connection,
    size = 8L,
    endian = "little"
  )
  close(connection)

  observed <- .tgc_read_binary_matrix(path, dim(values))
  expect_identical(observed, unname(values * 1.0))
  expect_error(.tgc_read_binary_matrix(path, c(4L, 4L)), "truncated")
})

test_that("Windows paths can be translated for the WSL backend", {
  skip_if(.Platform$OS.type != "windows")
  skip_if(!.tgc_wsl_available())

  translated <- .tgc_wsl_path(getwd())
  expect_match(translated, "^/")
  expect_identical(
    .tgc_wsl_capture(c("/usr/bin/test", "-d", translated))$status, 0L
  )
  expect_true(nzchar(.tgc_wsl_root()))

  spaced_dir <- file.path(tempdir(), "tidygenclust path with spaces")
  dir.create(spaced_dir, showWarnings = FALSE)
  expect_match(
    .tgc_wsl_path(spaced_dir),
    "/tidygenclust path with spaces$",
    fixed = FALSE
  )
})

test_that("WSL preserves command arguments and reports command failures", {
  skip_if(.Platform$OS.type != "windows")
  skip_if(!.tgc_wsl_available())

  values <- c("", "space name", 'double"quote', "O'Brien", "a&b")
  if (l10n_info()[["UTF-8"]]) {
    values <- c(values, "caf\u00e9")
  }
  observed <- .tgc_wsl_capture(c("/usr/bin/printf", "%s\\n", values))
  expect_identical(observed$status, 0L)
  expect_identical(observed$output, values)
  expect_identical(.tgc_wsl_capture("/usr/bin/false")$status, 1L)

  skip_if(!.tgc_wsl_env_exists("ctidygenclust"))
  expect_error(
    .tgc_wsl_run_env("ctidygenclust", "/usr/bin/false"), "exit status 1"
  )
  expect_error(
    .tgc_wsl_run_env("missing-test-environment", "/usr/bin/true"),
    "does not exist"
  )
})

test_that("installed WSL tools are discoverable", {
  skip_if(.Platform$OS.type != "windows")
  skip_if(!.tgc_wsl_available())
  skip_if(!.tgc_wsl_env_exists("ctidygenclust"))
  skip_if(!.tgc_wsl_env_exists("cclumppling"))

  expect_true(.tgc_tool_available("fastmixture"))
  expect_true(.tgc_tool_available("clumppling"))
})

test_that("WSL analyses preserve special Windows paths", {
  skip_on_cran()
  skip_if(.Platform$OS.type != "windows")
  skip_if(!.tgc_tool_available("fastmixture"))
  skip_if(!.tgc_tool_available("adamixture"))

  # Use the bundled data, with both input and output paths on a Windows drive.
  accent <- if (l10n_info()[["UTF-8"]]) "caf\u00e9 " else ""
  path <- tempfile(paste0("tidygenclust spaces & apostrophe' ", accent))
  dir.create(path)
  vcf <- system.file(
    "extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
    package = "tidypopgen"
  )
  x <- tidypopgen::gen_tibble(
    vcf, quiet = TRUE, backingfile = file.path(path, "anolis data"),
    parser = "cpp"
  )
  prefix <- tidypopgen::gt_as_plink(
    x, file = file.path(path, "plink data"), chromosomes_as_int = TRUE
  )
  prefix <- sub("\\.bed$", "", prefix)

  for (run_model in list(gt_fastmixture, gt_adamixture)) {
    result <- run_model(prefix, k = 2, no_freqs = FALSE, cv = 3)
    expect_equal(dim(result$Q[[1]]), c(nrow(x), 2L))
    expect_equal(ncol(result$P[[1]]), 2L)
    expect_true(all(is.finite(result$Q[[1]])))
    expect_true(all(is.finite(result$P[[1]])))
    expect_true(all(is.finite(result$cv)))
    expect_equal(unname(rowSums(result$Q[[1]])), rep(1, nrow(x)),
      tolerance = 1e-6
    )
  }

  # Supervised labels are a second input path passed through the bridge.
  labels_path <- file.path(path, "population labels.txt")
  labels <- c(1L, 2L, rep(0L, nrow(x) - 2L))
  writeLines(as.character(labels), labels_path)
  supervised <- gt_fastmixture(
    prefix, k = 2, supervised = labels_path, threads = 2
  )
  expect_true(supervised$Q[[1]][1, 1] > 0.99)
  expect_true(supervised$Q[[1]][2, 2] > 0.99)
  expect_equal(unname(rowSums(supervised$Q[[1]])), rep(1, nrow(x)),
    tolerance = 1e-6
  )
})

test_that("WSL Clumppling translates paths without changing option values", {
  skip_on_cran()
  skip_if(.Platform$OS.type != "windows")
  skip_if(!.tgc_tool_available("clumppling"))

  accent <- if (l10n_info()[["UTF-8"]]) "caf\u00e9 " else ""
  path <- tempfile(paste0("Q matrices & apostrophe' ", accent))
  dir.create(path)
  input <- file.path(path, "admixture")
  dir.create(input)
  for (k in 2:3) {
    q <- matrix(seq_len(12 * k), ncol = k)
    q <- q / rowSums(q)
    utils::write.table(
      q, file.path(input, paste0("run", k, ".Q")),
      row.names = FALSE, col.names = FALSE, quote = FALSE
    )
  }
  # The relative input path deliberately matches the default input_format.
  old_wd <- getwd()
  result <- tryCatch({
    setwd(path)
    gt_clumppling("admixture", output_path = "clump output")
  }, finally = setwd(old_wd))
  expect_s3_class(result, "gt_clumppling")
  expect_equal(result$N, 12L)
  expect_equal(result$K_range, 2:3)
})
