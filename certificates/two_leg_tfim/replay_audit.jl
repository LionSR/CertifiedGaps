# Read-only replay of the fixed seed from model inputs; the old certificate is used
# ONLY after fresh dense spectra, spectral radii, partition intervals and seeds exist.
include("src/Ladder.jl")
using .Ladder, ThermalMPOGap, ThermalMPOGap.Results, LinearAlgebra, Dates
BLAS.set_num_threads(1)
setprecision(BigFloat,256)
@assert VERSION == v"1.12.4"
@assert rounding(Float64) == RoundNearest && !get_zero_subnormals()
root=@__DIR__
original=joinpath(root,"results","certificate_g3_1.json")
output=joinpath(root,get(ARGS,1,"results/audit_2026-10-11_replay.json"))
@assert !isfile(output) "Audit outputs are immutable: choose a new output filename"
original_hash=sha256_file(original)
base_commit=strip(read(`git -C $root rev-parse HEAD`,String))
started=string(now(UTC)); t=time()
g=3//1; n=4; beta=BigInt(3)//BigInt(2)
@assert rp_identity(g) && iseven(n) && g>0
fresh=[full_spectrum(L,g;verified=true) for L in (n,2n)]
ln=logz_interval(fresh[1],beta); l2=logz_interval(fresh[2],beta)
l22=logz_interval(fresh[2],2beta)
ls=Ladder.down(()->l2[1]-2ln[2]); lt=Ladder.down(()->l22[1]-2l2[2])
es=Rational{BigInt}(-Ladder.down(()->expm1(ls)))
et=Rational{BigInt}(-Ladder.down(()->expm1(lt)))
@assert 0<=es<1 && 0<=et<1
result=certified_gap_bound(es,et,beta)
@assert result.bound !== nothing
# The old numerical values first enter here, exclusively for comparison.
old=read_exact_json(original)
@assert parse(Int,old["n0"])==n && parse(Rational{BigInt},old["beta"])==beta
@assert parse(Rational{BigInt},old["g"])==g
comparisons=[]
for (sp,osp) in zip(fresh,old["spectra"])
    @assert sp.L==parse(Int,osp["L"])
    @assert length(sp.sectors)==length(osp["sectors"])
    for (sec,osec) in zip(sp.sectors,osp["sectors"])
        @assert sec.parity==parse(Int,osec["parity"])
        @assert sec.character==parse(Int,osec["character"])
        @assert sec.dimension==parse(Int,osec["dimension"])
        newrad=parse(Rational{BigInt},sec.verification.radius)
        oldrad=parse(Rational{BigInt},osec["verification"]["radius"])
        push!(comparisons,(;L=sp.L,parity=sec.parity,character=sec.character,
            dimension=sec.dimension,eigenvalues_identical=sec.eigenvalues==Float64.(osec["eigenvalues"]),
            radius_identical=newrad==oldrad,radius_difference=string(newrad-oldrad),
            residual_identical=sec.verification.residual==osec["verification"]["residual"],
            orthogonality_identical=sec.verification.orthogonality==osec["verification"]["orthogonality"],
            construction_identical=sec.verification.construction==osec["verification"]["construction"]))
    end
end
oldes=parse(Rational{BigInt},old["epsilon_s"]); oldet=parse(Rational{BigInt},old["epsilon_t"])
hashes=merge(source_hashes(root;dirs=("src",),files=["replay_audit.jl","certify.jl","test/runtests.jl"]),
    Dict("../../julia/"*k=>v for (k,v) in source_hashes(joinpath(root,"../../julia"))),
    file_hashes(root,["../../julia/Project.toml","../../julia/Manifest.toml",
        "../../outputs/thermal-mpo-gap-scope.tex"]))
@assert all(sha256_file(joinpath(root,p))==h for (p,h) in old["metadata"]["source_sha256"])
@assert sha256_file(original)==original_hash
agreement=(;all_eigenvalues_identical=all(c.eigenvalues_identical for c in comparisons),
    all_radii_identical=all(c.radius_identical for c in comparisons),
    all_residual_bounds_identical=all(c.residual_identical for c in comparisons),
    all_orthogonality_bounds_identical=all(c.orthogonality_identical for c in comparisons),
    all_construction_bounds_identical=all(c.construction_identical for c in comparisons),
    epsilon_s_identical=es==oldes,epsilon_t_identical=et==oldet,
    epsilon_s_difference=string(es-oldes),epsilon_t_difference=string(et-oldet),
    bound_identical=result.bound==old["bound"],original_bound_supported=result.bound>=old["bound"])
assumptions=[
    "Pauli normalization; one rung cell d=4; periodic legs and onsite rung coupling of unit strength.",
    "Exact rational RP: C1=Z⊗I, C2=I⊗Z, c1=c2=1; real symmetric onsite term.",
    "Even n0=4 and beta=3//2; full dense character blocks exhaust 4^L at L=4,8.",
    "IEEE binary64 round-to-nearest conventional BLAS dot products with gradual underflow, no overflow.",
    "MPFR correctly rounded scalar operations at 256 bits; gamma_(2m+8) norm corrections and additive underflow term.",
    "Polar/Weyl whole-spectrum enclosure; directed lower numerator and upper denominator log-partition bounds.",
    "Connected negative single-spin-flip graph for g>0 gives finite-ring uniqueness by Perron–Frobenius.",
    "Uniform conclusion uses repository two-purity theorem and certified_gap_bound; no phase-diagram estimate is a premise."]
payload=(;status="fresh spectral verification complete",started_utc=started,finished_utc=string(now(UTC)),
    seconds=time()-t,base_commit,original_certificate_sha256=original_hash,
    original_artifact_preserved=true,source_hashes=hashes,
    g=string(g),n0=n,beta=string(beta),epsilon_s=string(es),epsilon_t=string(et),
    bound=result.bound,bound_exact_binary=string(Rational{BigInt}(result.bound)),recurrence=result,
    logZ=[map(x->string(Rational{BigInt}(x)),v) for v in (ln,l2,l22)],
    agreement,sector_comparisons=comparisons,spectra=fresh,theorem_assumptions=assumptions,
    metadata=run_metadata((;mode="fixed-seed replay; no search",bigfloat_precision=precision(BigFloat),
        blas=string(BLAS.get_config()),rounding="nearest",gradual_underflow=true);
        source_sha256=hashes,packages=(ThermalMPOGap,),status="certified replay",certified=true),
    reproduction="OPENBLAS_NUM_THREADS=1 julia --project=julia certificates/two_leg_tfim/replay_audit.jl results/replay_review.json",
    witness_scope="All spectra and outward scalar bounds retained. Eigenvectors regenerated by the reproduction command; no old radii used in fresh computations.")
write_result(output,payload)
println("REPLAY AGREEMENT: ",agreement)
println("REPLAY bound=$(result.bound) elapsed=$(time()-t) output=$output")
