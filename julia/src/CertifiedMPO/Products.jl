"""
    CertifiedGaps.CertifiedMPO.Products

Certified zip-up compression of an MPO product A·B with unequal auxiliary dimensions (in
rational or interval candidate arithmetic), its replay, its witness encoding, and the
three-layer overlap ⟨C, A·B⟩ in Arb. Moved verbatim from spin2 `Spin2MPOProducts`; `MS` is
`CertifiedMPO.Squaring` instead of a private include copy.
"""
module Products

using LinearAlgebra
import Nemo
using ..CertifiedMPO: nemo_array
import ..Squaring
const MS=Squaring
const MN=MS.MN
const F=MN.F

"Three-layer overlap <C,A*B>, without constructing product cores or a full operator."
function product_inner_ball(C, A, B; weights = nothing, precision = 192)
    @assert length(C)==length(A)==length(B) && !isempty(C)
    RR=Nemo.ArbField(precision)
    environment=reshape([RR(1)], 1, 1, 1)
    for (c, a, b) in zip(C, A, B)
        lc, d, di, rc=size(c)
        la, da, mid, ra=size(a)
        lb, midb, dib, rb=size(b)
        @assert d==da && di==dib && mid==midb && size(environment)==(lc, la, lb)
        blocks(X) = [
            any(!iszero, X[:, s, t, :]) ? Nemo.matrix(RR, X[:, s, t, :]) : nothing for
            s in axes(X, 2), t in axes(X, 3)
        ]
        cb, ab, bb=blocks(c), blocks(a), blocks(b)
        incoming=Nemo.matrix(RR, reshape(permutedims(environment, (1, 3, 2)), lc*lb, la))
        outgoing=Nemo.zero_matrix(RR, rc*ra, rb)
        for s in 1:d, u in 1:mid
            isnothing(ab[s, u]) && continue
            contracted=incoming*ab[s, u]
            first=Nemo.matrix(RR, reshape(nemo_array(contracted, RR), lc, lb*ra))
            for t in 1:di
                (isnothing(cb[s, t]) || isnothing(bb[u, t])) && continue
                second=transpose(cb[s, t])*first
                raw=reshape(nemo_array(second, RR), rc, lb, ra)
                third=Nemo.matrix(RR, reshape(permutedims(raw, (1, 3, 2)), rc*ra, lb))*bb[u, t]
                weight=isnothing(weights) ? F(1) : F(weights[t])/F(weights[s])
                outgoing+=RR(weight)*third
            end
        end
        environment=reshape(nemo_array(outgoing, RR), rc, ra, rb)
    end
    @assert size(environment)==(1, 1, 1)
    return only(environment)
end

"Certified zip-up compression of A*B, including unequal auxiliary dimensions."
function product_candidate(
    A,
    B;
    maxdim = 16,
    bits = 48,
    precision = 192,
    weights = nothing,
    first_gauges = nothing,
    candidate_arithmetic = :rational,
    physical_charges = nothing,
    progress = nothing,
)
    @assert length(A)==length(B) && !isempty(A)
    @assert size(A[1], 1)==size(B[1], 1)==size(A[end], 4)==size(B[end], 4)==1
    @assert candidate_arithmetic in (:rational, :interval)
    RR = Nemo.ArbField(precision)
    first=isnothing(first_gauges) ? A :
          MN.gauge_transform(
        A,
        first_gauges;
        interval_precision = candidate_arithmetic == :interval ? precision : nothing,
    )
    charge_a = isnothing(physical_charges) ? nothing : MS.CB.infer_charges(first, physical_charges)
    charge_b = isnothing(physical_charges) ? nothing : MS.CB.infer_charges(B, physical_charges)
    prefix_charges = [0]
    operator_bounds=[MN.core_norm_upper(c, :operator; weights, precision) for c in first]
    hs_bounds=[MN.core_norm_upper(c, :right_hs; weights, precision) for c in B]
    n=length(A)
    tails=ones(F, n+1)
    for j in n:-1:1
        tails[j]=tails[j+1]*operator_bounds[j]*hs_bounds[j]
    end
    bridge=ones(F, 1, 1, 1)
    cores=Array{F,4}[]
    bridges=Array{F,3}[]
    residual_bounds, prefix_bounds=F[], F[]
    error, prefix=F(0), F(1)
    for j in 1:n
        # Interior products may use outward-rounded enclosures. The terminal
        # equality remains rational, closing the exact telescoping identity.
        field = candidate_arithmetic == :interval && j < n ? RR : Nemo.QQ
        Z =
            isnothing(physical_charges) ? MS.product_bridge(bridge, first[j], B[j]; field) :
            MS.CB.product_bridge(
                bridge,
                first[j],
                B[j],
                prefix_charges,
                charge_a[j],
                charge_a[j+1],
                charge_b[j],
                charge_b[j+1],
                physical_charges;
                field,
            )
        chi, d, di, ra, rb=size(Z)
        @assert d==di
        M=reshape(Z, chi*d*d, ra*rb)
        if j==n
            push!(cores, reshape(M, chi, d, d, 1))
            isnothing(progress) || progress((; stage = :compression, site = j, total = n))
            break
        end
        metric=MN.hs_row_metric(chi, d, weights)
        root=sqrt.(Float64.(metric))
        labels =
            isnothing(physical_charges) ? nothing :
            MS.CB.product_labels(prefix_charges, charge_a[j+1], charge_b[j+1], physical_charges)
        decomposition =
            isnothing(labels) ? MN.numerical_svd(Float64.(M) .* root) :
            MS.CB.truncated_svd(
                Float64.(M) .* root,
                labels...;
                maxdim,
                svd_function = MN.numerical_svd,
            )
        isnothing(labels) || (prefix_charges = decomposition.charges)
        rank=min(maxdim, length(decomposition.S))
        U=MS.relative_quantize(decomposition.U[:, 1:rank] ./ root; bits)
        R=MS.relative_quantize(Diagonal(decomposition.S[1:rank])*decomposition.Vt[1:rank, :]; bits)
        residual=M-MS.matrix_product(U, R, field)
        delta=MS.frobenius_upper(residual, field; metric)
        C=reshape(U, chi, d, d, rank)
        alpha=MN.core_norm_upper(C, :left_hs; weights, precision)
        error+=prefix*delta*tails[j+1]
        prefix*=alpha
        bridge=reshape(R, rank, ra, rb)
        push!(cores, C)
        push!(bridges, bridge)
        push!(residual_bounds, delta)
        push!(prefix_bounds, alpha)
        isnothing(progress) || progress((; stage = :compression, site = j, total = n))
    end
    return (;
        cores,
        bridges,
        residual_bounds,
        prefix_bounds,
        operator_bounds,
        hs_bounds,
        error,
        weights,
        first_gauges,
    )
end

"Replay the local product identities and all norm inequalities, without an SVD."
function verify_product(
    A,
    B,
    proof;
    precision = 192,
    residual_arithmetic = :rational,
    progress = nothing,
)
    @assert residual_arithmetic in (:rational, :interval)
    RR = Nemo.ArbField(precision)
    n=length(A)
    @assert length(B)==length(proof.cores)==length(proof.operator_bounds)==length(proof.hs_bounds)==n
    @assert length(proof.bridges)==length(proof.residual_bounds)==length(proof.prefix_bounds)==n-1
    weights=proof.weights
    first=isnothing(proof.first_gauges) ? A :
          MN.gauge_transform(
        A,
        proof.first_gauges;
        interval_precision = residual_arithmetic == :interval ? precision : nothing,
    )
    tails=ones(F, n+1)
    for j in n:-1:1
        @assert MN.verify_core_bound(
            first[j],
            :operator,
            proof.operator_bounds[j];
            weights,
            precision,
        )
        @assert MN.verify_core_bound(B[j], :right_hs, proof.hs_bounds[j]; weights, precision)
        tails[j]=tails[j+1]*proof.operator_bounds[j]*proof.hs_bounds[j]
    end
    bridge=ones(F, 1, 1, 1)
    prefix, error=F(1), F(0)
    for j in 1:n
        field = residual_arithmetic == :interval && j < n ? RR : Nemo.QQ
        Z=MS.product_bridge(bridge, first[j], B[j]; field)
        chi, d, di, ra, rb=size(Z)
        @assert d==di
        M=reshape(Z, chi*d*d, ra*rb)
        C=proof.cores[j]
        @assert size(C)[1:3]==(chi, d, d)
        if j==n
            @assert ra==rb==size(C, 4)==1 && reshape(C, size(M))==M
            isnothing(progress) || progress((; stage = :residuals, site = j, total = n))
            break
        end
        rank=size(C, 4)
        next=proof.bridges[j]
        @assert size(next)==(rank, ra, rb)
        U, R = reshape(C, chi*d*d, rank), reshape(next, rank, ra*rb)
        residual=M-MS.matrix_product(U, R, field)
        metric=MN.hs_row_metric(chi, d, weights)
        delta, alpha=proof.residual_bounds[j], proof.prefix_bounds[j]
        @assert delta >= 0
        if field === Nemo.QQ
            @assert delta^2 >= sum(abs2.(residual) .* metric)
        elseif !(RR(delta)^2 > sum(abs2.(residual) .* RR.(metric)))
            # Inconclusive intervals, including an exactly zero residual,
            # require the rational identity; they are never accepted directly.
            exact_core =
                isnothing(proof.first_gauges) ? A[j] :
                MN.exact_gauge_core(A[j], proof.first_gauges[j], proof.first_gauges[j+1])
            exact = reshape(MS.product_bridge(bridge, exact_core, B[j]), size(M))
            exact_residual = exact-MS.rational_product(U, R)
            @assert delta^2 >= sum(abs2.(exact_residual) .* metric)
        end
        @assert MN.verify_core_bound(C, :left_hs, alpha; weights, precision)
        error+=prefix*delta*tails[j+1]
        prefix*=alpha
        bridge=next
        isnothing(progress) || progress((; stage = :residuals, site = j, total = n))
    end
    @assert error==proof.error
    return true
end

function encode(A, B, proof)
    return Dict(
        "first_input"=>MN.encode_cores(A),
        "second_input"=>MN.encode_cores(B),
        "output"=>MN.encode_cores(proof.cores),
        "bridges"=>[
            Dict("dimensions"=>collect(size(R)), "entries"=>string.(vec(R))) for R in proof.bridges
        ],
        "residual_bounds"=>string.(proof.residual_bounds),
        "prefix_bounds"=>string.(proof.prefix_bounds),
        "operator_bounds"=>string.(proof.operator_bounds),
        "hs_bounds"=>string.(proof.hs_bounds),
        "weights"=>isnothing(proof.weights) ? nothing : string.(proof.weights),
        "first_gauges"=>isnothing(proof.first_gauges) ? nothing :
                        MN.encode_gauges(proof.first_gauges),
        "error"=>string(proof.error),
    )
end

function decode(data)
    A, B=MN.decode_cores(data["first_input"]), MN.decode_cores(data["second_input"])
    vectors(key) = MN.parse_fraction.(data[key])
    proof=(;
        cores = MN.decode_cores(data["output"]),
        bridges = [
            reshape(MN.parse_fraction.(R["entries"]), Int.(R["dimensions"])...) for
            R in data["bridges"]
        ],
        residual_bounds = vectors("residual_bounds"),
        prefix_bounds = vectors("prefix_bounds"),
        operator_bounds = vectors("operator_bounds"),
        hs_bounds = vectors("hs_bounds"),
        weights = isnothing(data["weights"]) ? nothing : vectors("weights"),
        first_gauges = isnothing(data["first_gauges"]) ? nothing :
                       MN.decode_gauges(data["first_gauges"]),
        error = MN.parse_fraction(data["error"]),
    )
    return A, B, proof
end

end
