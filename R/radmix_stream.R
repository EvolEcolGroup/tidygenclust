#' Whether gt_radmix() streams its data (internal)
#'
#' Adds the elements `stream`, `chunk_q` and `chunk_bytes` to the data of
#' `gt_radmix()` (the list returned by `gt_radmix_data()`), which tell the
#' C++ passes whether to walk the loci in chunks, prefetching the next one
#' and alternating the direction of the passes (see "Streaming
#' larger-than-memory data" in `src/radmix_geno_access.h`). Streaming is for
#' data larger than the memory, whose backing file cannot stay in the page
#' cache from one pass to the next. Genotypes that fit in memory get
#' `chunk_q`: the Q side alone runs in chunks, when it has one segment
#' (many individuals; "Chunks in memory" in the same file), which is faster
#' there. Neither ever changes the results, which are bit-identical either
#' way.
#'
#' The option `tidygenclust.radmix_stream` decides: `TRUE` (stream) or
#' `FALSE` (neither: the usual layout), or `NULL` (the default) to stream
#' when the bytes of the backing file read by a pass exceed `0.8` times the
#' memory available (`radmix_mem_available()`), and to use `chunk_q`
#' otherwise. The option `tidygenclust.radmix_chunk_mb` sets the size of
#' the chunks (in MB; by default, when streaming, a sixteenth of the memory
#' available, between 16 and 256 MB, and in memory
#' `radmix_chunk_q_mb()`).
#'
#' @param data the data of `gt_radmix()`
#' @returns `data`, with `stream` and `chunk_q` (`TRUE` or `FALSE`),
#'   `chunk_bytes` and `data_bytes` (the bytes read by a pass, or `NA` if
#'   unknown).
#' @keywords internal
#' @noRd
radmix_stream_mode <- function(data) {
  bytes <- radmix_data_bytes(data)
  mem <- radmix_mem_available()
  opt <- getOption("tidygenclust.radmix_stream", default = NULL)
  stream <- if (is.null(opt)) {
    is.finite(bytes) && bytes > 0.8 * mem
  } else {
    isTRUE(opt)
  }
  # genotypes in memory: the Q side in chunks (with one segment), unless
  # the option asks for the usual layout
  chunk_q <- !stream && !isFALSE(opt) && !identical(data$kind, "gl_source")
  chunk_mb <- getOption("tidygenclust.radmix_chunk_mb", default = NULL)
  if (is.null(chunk_mb)) {
    chunk_mb <- if (stream) {
      # a sixteenth of the memory available, between 16 and 256 MB: the
      # chunk being computed and the one being prefetched take memory from
      # the page cache, which keeps the data that the next pass reads first
      min(256, max(16, mem / 16 / 2^20))
    } else {
      radmix_chunk_q_mb()
    }
  }
  if (!(length(chunk_mb) == 1 && is.numeric(chunk_mb) && chunk_mb > 0)) {
    stop("option 'tidygenclust.radmix_chunk_mb' must be a positive number",
      call. = FALSE
    )
  }
  data$stream <- stream
  data$chunk_q <- chunk_q
  data$chunk_bytes <- chunk_mb * 2^20
  data$data_bytes <- bytes
  data
}

# The size of the chunks of the Q side of genotypes in memory, in MB (see
# radmix_stream_mode()). On manica01 (27 Nov 2026, K = 5, 8 cores), chunks of
# 8 to 64 MB were all 15-19% faster than the usual layout and within about 2%
# of each other: 64 MB was best with 250,000 individuals (smaller chunks hold
# a single tile of 64 loci there) and 1.7% behind 16 MB with 10,000
# (dev_notes/radmix_streaming.md, "Chunks in memory")
radmix_chunk_q_mb <- function() 64

# The bytes of the backing file that a pass reads: whole columns of the FBM
# for a gen_tibble, and locus_bytes per locus for a source of genotype
# likelihoods that gives it (NA otherwise)
radmix_data_bytes <- function(data) {
  if (identical(data$kind, "gl_source")) {
    if (is.null(data$locus_bytes)) {
      return(NA_real_)
    }
    return(as.numeric(data$locus_bytes) * data$n_loci)
  }
  as.numeric(data$fbm$nrow) * length(data$cols)
}

#' The memory available to this process (internal)
#'
#' On Linux, the smaller of `MemAvailable` in `/proc/meminfo` (free memory
#' plus the page cache and other memory that can be reclaimed) and the room
#' left under the memory limit of the process's control group (cgroup v2
#' `memory.max`, or v1 `memory.limit_in_bytes`, minus the memory used
#' outside the page cache). Elsewhere `Inf`, so that gt_radmix() does not
#' stream unless asked to.
#'
#' @returns the number of bytes, or `Inf` if unknown.
#' @keywords internal
#' @noRd
radmix_mem_available <- function() {
  meminfo <- "/proc/meminfo"
  if (!file.exists(meminfo)) {
    return(Inf)
  }
  lines <- readLines(meminfo, warn = FALSE)
  avail <- grep("^MemAvailable:", lines, value = TRUE)
  mem <- if (length(avail) == 1) {
    as.numeric(sub("^MemAvailable:\\s+(\\d+).*", "\\1", avail)) * 1024
  } else {
    Inf
  }
  min(mem, radmix_cgroup_available(), na.rm = TRUE)
}

# The room under the memory limit of the cgroup of this process (Inf if
# there is no limit or it cannot be read): the limit minus the memory
# charged to the cgroup, not counting its page cache (which the kernel
# reclaims under the limit)
radmix_cgroup_available <- function() {
  cg <- tryCatch(readLines("/proc/self/cgroup", warn = FALSE),
    error = function(e) character(0)
  )
  read_num <- function(f) {
    v <- tryCatch(readLines(f, n = 1, warn = FALSE),
      error = function(e) NA_character_
    )
    suppressWarnings(as.numeric(v))
  }
  read_stat <- function(f, key) {
    v <- tryCatch(readLines(f, warn = FALSE),
      error = function(e) character(0)
    )
    v <- grep(paste0("^", key, " "), v, value = TRUE)
    if (length(v) != 1) {
      return(0)
    }
    as.numeric(sub(".* ", "", v))
  }
  # cgroup v1: the memory controller
  v1 <- grep("^[0-9]+:([^:]*,)?memory(,[^:]*)?:", cg, value = TRUE)
  if (length(v1) == 1) {
    dir <- file.path("/sys/fs/cgroup/memory", sub("^[^:]*:[^:]*:", "", v1))
    limit <- read_num(file.path(dir, "memory.limit_in_bytes"))
    used <- read_num(file.path(dir, "memory.usage_in_bytes"))
    cache <- read_stat(file.path(dir, "memory.stat"), "cache")
    if (is.finite(limit) && is.finite(used) && limit < 2^60) {
      return(limit - (used - cache))
    }
    return(Inf)
  }
  # cgroup v2
  v2 <- grep("^0::", cg, value = TRUE)
  if (length(v2) == 1) {
    dir <- file.path("/sys/fs/cgroup", sub("^0::", "", v2))
    limit <- read_num(file.path(dir, "memory.max"))
    used <- read_num(file.path(dir, "memory.current"))
    cache <- read_stat(file.path(dir, "memory.stat"), "file")
    if (is.finite(limit) && is.finite(used)) {
      return(limit - (used - cache))
    }
  }
  Inf
}
