# Summary of run_passes.sh: for each dataset, memory cap and disk, and
# version, the medians of each kind of pass (update_q, update_q_loglik:
# Q side; update_p: P side), from the third pass on (the first two, EM and
# a fused Q update, start from a cold cache): time, reads of the file, read
# throughput and CPU use.
# Usage: Rscript summarise_passes.R <passes.csv>
x <- utils::read.csv(commandArgs(TRUE)[1])
x <- x[x$i >= 3, ]
x$config <- sub(":[^:]*$", "", sub("^[^:]*:", "", x$label))
x$data <- sub(":.*", "", x$label)
x$disk <- ifelse(x$throttle_mbps == 0, "as is", paste(x$throttle_mbps, "MB/s"))
agg <- stats::aggregate(
  cbind(elapsed_s, file_reads, read_mbps, cpu_pct) ~ data + config + disk +
    mode + pass,
  x, stats::median
)
agg <- agg[order(agg$data, agg$config, agg$disk, agg$pass, agg$mode), ]
agg$elapsed_s <- round(agg$elapsed_s, 2)
agg$file_reads <- round(agg$file_reads, 2)
agg$read_mbps <- round(agg$read_mbps)
agg$cpu_pct <- round(agg$cpu_pct)
print(agg, row.names = FALSE)
