# Reference computations with OpenADMIXTURE.jl, run by
# data-raw/radmix/regenerate_julia_fixtures.R to create the fixtures used by
# tests/testthat/test_radmix_vs_julia.R.
#
# Usage:
#   julia --threads=1,0 openadmixture_reference.jl <bed file> <q0 csv> <p0 csv> \
#     <output dir> <n_iter>
#
# (--threads=1,0 starts Julia with a single thread and no interactive thread
# pool: OpenADMIXTURE.jl indexes its per-thread buffers with threadid(), which
# breaks with the interactive thread added by default since Julia 1.12.)
#
# q0 (K x I) and p0 (K x J) are the starting values, written by R as plain
# comma-separated matrices. Starting from them, we call the internal functions
# of OpenADMIXTURE.jl that gt_radmix() reimplements, and write each result
# to <output dir> for R to compare:
#   loglik.csv        log-likelihood at (q0, p0)
#   em_q.csv/em_p.csv one FRAPPE EM step (em_q!/em_p! in :base mode)
#   upd_q.csv         one Newton update of Q (update_q!)
#   upd_p.csv         one Newton update of P given the new Q (update_p!)
#   run_q.csv/run_p.csv/run_ll.csv
#                     the full quasi-Newton algorithm (admixture_qn!) started
#                     from (q0, p0), with the log-likelihood of each iteration
#
# Note that SnpArrays counts copies of the second allele of the bim file,
# whereas gt_radmix() counts the first (allele_alt in tidypopgen), so the
# P matrices exchanged with R are frequencies of the second allele (R flips
# them as needed).

using OpenADMIXTURE
using SnpArrays
using DelimitedFiles
using LinearAlgebra
using Random

const OA = OpenADMIXTURE

bed_file, q0_file, p0_file, out_dir = ARGS[1:4]
n_iter = parse(Int, ARGS[5])

# genotypes, read exactly as OpenADMIXTURE.jl's driver does
g = SnpArray(bed_file)
g_la = SnpLinAlg{Float64}(g)
I, J = size(g_la)

q0 = readdlm(q0_file, ',', Float64)
p0 = readdlm(p0_file, ',', Float64)
K = size(q0, 1)
@assert size(q0) == (K, I)
@assert size(p0) == (K, J)

# skipmissing = true, as in run_admixture(); the rng is only used to fill the
# buffers with random values, which we overwrite
d = OA.AdmixData{Float64}(I, J, K, 3; skipmissing = true,
                          rng = MersenneTwister(1))

function reset!(d)
    d.q .= q0
    d.p .= p0
    d.ll_new = NaN
end

writecsv(name, x) = writedlm(joinpath(out_dir, name), x, ',')

# log-likelihood
reset!(d)
ll = OA.loglikelihood(g_la, d.q, d.p, d.qp_small, K, d.skipmissing)
writecsv("loglik.csv", [ll])

# one EM step: both em_q! and em_p! read (q, p) and write (q_next, p_next)
reset!(d)
OA.em_q!(d, g_la)
OA.em_p!(d, g_la)
writecsv("em_q.csv", d.q_next)
writecsv("em_p.csv", d.p_next)

# one round of block updates: update_q! writes q_next from (q, p), then
# update_p! writes p_next from (q_next, p)
reset!(d)
OA.update_q!(d, g_la)
OA.update_p!(d, g_la)
writecsv("upd_q.csv", d.q_next)
writecsv("upd_p.csv", d.p_next)

# the full algorithm from (q0, p0), as run_admixture() does when both init_q
# and init_p are given (no EM priming). admixture_qn! prints the
# log-likelihood of each iteration, which we capture and parse.
reset!(d)
lls = Float64[]
out = mktemp() do path, io
    redirect_stdout(io) do
        OA.admixture_qn!(d, g_la, n_iter, 1e-7; mode = :ZAL)
    end
    flush(io)
    read(path, String)
end
for m in eachmatch(r"LogLikelihood = (-?[0-9.eE+-]+)", out)
    push!(lls, parse(Float64, m.captures[1]))
end
writecsv("run_q.csv", d.q)
writecsv("run_p.csv", d.p)
writecsv("run_ll.csv", lls)
println("done")
