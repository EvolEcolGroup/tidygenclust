# OpenADMIXTURE.jl (Ko et al. 2023), the Julia implementation of
# OpenADMIXTURE, run from the command line by admixture_benchmark_cli.sh (in
# the same directory as this script) for the benchmark shown in the article
# admixture_benchmark.Rmd, as:
#
#   julia --threads=<cores>,0 admixture_benchmark_openadmixture.jl \
#     <bed file> <K> <seed> <output prefix> [<warm-up bed file>]
#
# It runs OpenADMIXTURE.jl on all the loci of the bed file (without its
# selection of ancestry informative markers), with its default settings,
# which are those of gt_radmix(). It writes the Q and P matrices to
# <output prefix>.Q and <output prefix>.P, as ADMIXTURE does, and prints its
# progress, the log-likelihood and the time taken by the fit (including
# reading the bed file and writing Q and P).
#
# Julia compiles the code of a function the first time it is called. With a
# warm-up bed file, OpenADMIXTURE.jl is first run on it, for one iteration,
# so that the fit is timed with the code already compiled, as in a Julia
# session in which OpenADMIXTURE.jl has already been used.
#
# (--threads=<cores>,0 starts Julia without an interactive thread:
# OpenADMIXTURE.jl indexes its per-thread buffers with threadid(), which
# breaks with the interactive thread added by default since Julia 1.12.)

## ---- julia-run ----
using OpenADMIXTURE, Printf, Random

bed, k, seed, prefix = ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3]), ARGS[4]
warmup = length(ARGS) >= 5 ? ARGS[5] : nothing

# a matrix written with one row per line, as ADMIXTURE writes its Q and P
# matrices (6 decimals)
function write_matrix(file, m)
    open(file, "w") do io
        for i in axes(m, 1)
            println(io, join([@sprintf("%.6f", x) for x in m[i, :]], " "))
        end
    end
end

if warmup !== nothing
    OpenADMIXTURE.run_admixture(warmup, k; rng=Xoshiro(seed), admix_n_iter=1,
        progress_bar=false)
end

seconds = @elapsed begin
    # sparsity = nothing (the default): all the loci are used
    d, _, _ = OpenADMIXTURE.run_admixture(bed, k; rng=Xoshiro(seed), progress_bar=false)
    # d.q and d.p have one column per individual and per locus
    write_matrix(prefix * ".Q", permutedims(d.q))
    write_matrix(prefix * ".P", permutedims(d.p))
end
println("Loglikelihood: ", @sprintf("%.6f", d.ll_new))
println("Fit time (s): ", @sprintf("%.3f", seconds))
