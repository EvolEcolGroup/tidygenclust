# Plots of the multi-core scaling (phase E of the verification): speed-up
# against the number of cores, for the kernels (scaling_kernels.R) and for
# whole fits (bench_fits.R run on several core counts), base against comb.
#
#   Rscript plot_scaling.R <kernels dir> <fits dir> <out png>
#
# e.g. Rscript plot_scaling.R desktop_results/E_scaling desktop_results/E_fits \
#        desktop_results/E_scaling/scaling.png
# Speed-ups are relative to the 1-core time of the same build, binding and
# case. Kernels: the threads bound to cores (OMP_PROC_BIND=close,
# OMP_PLACES=cores) if those runs exist, otherwise all runs.
suppressMessages({
  library(ggplot2)
  library(patchwork)
})
args <- commandArgs(TRUE)
kern <- read.csv(file.path(args[1], "scaling_kernels.csv"))
fits <- read.csv(file.path(args[2], "fits.csv"))
out <- args[3]

if (any(grepl("close", kern$binding))) kern <- kern[grepl("close", kern$binding), ]
kern$case <- paste0(kern$data, ", K = ", kern$k)
kern$kernel <- factor(kern$kernel, c(
  "update_q", "update_p", "em", "loglik", "update_q_cv", "update_p_cv"
))
cores <- sort(unique(c(kern$n_cores, fits$n_cores)))

theme_set(theme_bw(base_size = 10))
ideal <- geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50")
p1 <- ggplot(kern, aes(n_cores, speedup, colour = label)) +
  ideal +
  geom_line() +
  geom_point(size = 1) +
  facet_grid(case ~ kernel) +
  scale_x_continuous(breaks = cores) +
  scale_y_continuous(breaks = cores) +
  coord_equal() +
  labs(
    x = "cores", y = "speed-up over 1 core", colour = "build",
    title = "Kernels (one pass, fixed Q and P)"
  )

fm <- aggregate(seconds ~ label + data + k + n_cores, fits, median)
fm$t1 <- fm$seconds[match(
  paste(fm$label, fm$data, fm$k, 1),
  paste(fm$label, fm$data, fm$k, fm$n_cores)
)]
fm$speedup <- fm$t1 / fm$seconds
fm$case <- paste0(fm$data, ", K = ", fm$k)
p2 <- ggplot(fm, aes(n_cores, speedup, colour = label)) +
  ideal +
  geom_line() +
  geom_point(size = 1) +
  facet_wrap(~case, nrow = 1) +
  scale_x_continuous(breaks = cores) +
  scale_y_continuous(breaks = cores) +
  coord_equal() +
  labs(
    x = "cores", y = "speed-up over 1 core", colour = "build",
    title = "Whole fits (gt_radmix(), seed 1)"
  )
p3 <- ggplot(fm, aes(n_cores, seconds, colour = label)) +
  geom_line() +
  geom_point(size = 1) +
  facet_wrap(~case, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = cores) +
  expand_limits(y = 0) +
  labs(x = "cores", y = "seconds", colour = "build", title = "Whole fits: time")
ggsave(out, p1 / p2 / p3 + plot_layout(heights = c(3, 1.2, 1.2)),
  width = 11, height = 13, dpi = 110
)
