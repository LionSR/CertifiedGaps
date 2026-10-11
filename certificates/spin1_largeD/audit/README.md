# Independent review and fresh spectral replay

`INDEPENDENT_REVIEW.md` records a separately delegated mathematical/code review.
`reviewer_checks.jl` adds an exact Cartesian-to-normalized-Sz bridge and 256-bit
residual diagnostics. The diagnostics supplement the norm proof, not replace it.

Historical `check_pinned.jl` and `test/runtests.jl` consume stored spectral radii.
They do not recompute eigenvectors, residuals, or orthogonality defects. This
separate replay rebuilds all model blocks at L=6 and L=12, recomputes every dense
eigensystem, and computes new construction, residual and orthogonality bounds.
It uses current `CertifiedGaps`, without repeating the numerical search.

From repository root with Julia 1.12.4 and dependencies installed:

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=julia certificates/spin1_largeD/audit/reviewer_checks.jl
julia --project=julia certificates/spin1_largeD/audit/check.jl
```

The second command is a fast check of this new stored evidence, including source
and block-file hashes, complete multiplicity accounting, partition enclosures and
rational recurrence inputs. Like any stored-evidence checker it does not itself
repeat dense ED. To repeat that expensive stage, run these four independent
processes (each with one BLAS thread), then assemble:

```sh
julia --project=julia certificates/spin1_largeD/audit/replay_worker.jl 4//1 heavy
julia --project=julia certificates/spin1_largeD/audit/replay_worker.jl 4//1 mixed
julia --project=julia certificates/spin1_largeD/audit/replay_worker.jl 1603//500 heavy
julia --project=julia certificates/spin1_largeD/audit/replay_worker.jl 1603//500 mixed
julia --project=julia certificates/spin1_largeD/audit/assemble_replay.jl
```

Use a disposable checkout to retain this run's evidence. The scripts write only
`audit/results/`; they never overwrite original `../results/`. The worker reads
no historical spectrum or radius. The assembler reads historical data only for
comparison after computing a fresh certified bound.

The charge groups are {1,4,…,12} and {0,2,3}; the latter also recomputes L=6.
Charge and momentum multiplicities recover all 531441 levels from 92 L=12
blocks. Four processes were selected to fit the four-CPU, 16 GiB environment
while keeping each invocation below approximately 45 minutes.

`initial-driver-failure/` preserves an aborted first launch: adding string
metadata to the eigensolver's real-valued dictionary raised a type error on the
first L=12 block. No L=12 evidence was accepted from that launch. The corrected
driver converts the metadata dictionary to `Dict{String,Any}`; all four workers
restart from model inputs in an empty output directory.
