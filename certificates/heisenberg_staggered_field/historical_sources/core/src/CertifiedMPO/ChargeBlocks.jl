"""
Exact U(1) support and sectorwise numerical proposals for neutral MPOs in
(left, output, input, right) order. No tolerance is used to discard entries.
The generic residual verifiers remain independent of this optimization.
"""
module ChargeBlocks

using LinearAlgebra
import Nemo

"Infer integer bond charges from every nonzero, including disconnected support."
function infer_charges(cores, physical)
    @assert !isempty(cores) && all(x -> x isa Integer, physical)
    dims = vcat(size(first(cores), 1), [size(c, 4) for c in cores])
    @assert first(dims) == last(dims) == 1
    @assert all(
        size(c) == (dims[j], length(physical), length(physical), dims[j+1]) for
        (j, c) in enumerate(cores)
    )
    offsets = cumsum(vcat(0, dims))
    parent = collect(1:sum(dims))
    potential = zeros(Int, length(parent)) # q_node - q_parent
    function root(i)
        parent[i] == i && return i
        p = parent[i]
        parent[i] = root(p)
        potential[i] += potential[p]
        return parent[i]
    end
    function constrain(a, b, difference)
        ra, rb = root(a), root(b)
        if ra == rb
            potential[b]-potential[a] == difference || throw(
                ArgumentError("MPO support is not charge adapted; use a charge-preserving gauge"),
            )
        else
            parent[rb] = ra
            potential[rb] = difference+potential[a]-potential[b]
        end
    end
    constrain(1, offsets[end-1]+1, 0)
    for (j, c) in enumerate(cores), index in CartesianIndices(c)
        iszero(c[index]) && continue
        l, s, t, r = Tuple(index)
        constrain(offsets[j]+l, offsets[j+1]+r, physical[s]-physical[t])
    end
    roots = root.(eachindex(parent))
    origin = potential[1]
    charges = [potential[i]-(roots[i] == roots[1] ? origin : 0) for i in eachindex(parent)]
    return [charges[(offsets[j]+1):offsets[j+1]] for j in eachindex(dims)]
end

groups(labels) = Dict(q => findall(==(q), labels) for q in sort(unique(labels)))

"Check exact zeros between different sectors before skipping their arithmetic."
function check_support(A, rows, cols)
    @assert size(A) == (length(rows), length(cols))
    for j in axes(A, 2), i in axes(A, 1)
        rows[i] == cols[j] ||
            iszero(A[i, j]) ||
            throw(ArgumentError("Nonzero outside the declared charge sectors"))
    end
    return true
end

"Exact or outward-rounded product, evaluated only inside matching charge blocks."
function matrix_product(A, B, rows, inner, cols; field = Nemo.QQ)
    check_support(A, rows, inner)
    check_support(B, inner, cols)
    rg, kg, cg = groups(rows), groups(inner), groups(cols)
    zero_entry = field === Nemo.QQ ? big(0)//big(1) : field(0)
    C = fill(zero_entry, length(rows), length(cols))
    for q in sort(collect(intersect(keys(rg), keys(kg), keys(cg))))
        ri, ki, ci = rg[q], kg[q], cg[q]
        block = Nemo.matrix(field, A[ri, ki])*Nemo.matrix(field, B[ki, ci])
        for (b, j) in enumerate(ci), (a, i) in enumerate(ri)
            C[i, j] = field === Nemo.QQ ? Rational{BigInt}(block[a, b]) : block[a, b]
        end
    end
    return C
end

"Contract A·B without multiplying charge-forbidden matrix blocks."
function product_bridge(
    bridge,
    A,
    B,
    prefix,
    leftA,
    rightA,
    leftB,
    rightB,
    physical;
    field = Nemo.QQ,
)
    l, d, mid, r = size(A)
    lb, mid2, di, rb = size(B)
    chi = length(prefix)
    @assert d == mid == mid2 == di == length(physical)
    @assert size(bridge) == (chi, l, lb)
    X = reshape(permutedims(bridge, (1, 3, 2)), chi*lb, l)
    Y = reshape(A, l, d*mid*r)
    rows = vec([p-b for p in prefix, b in leftB])
    cols = vec([a-s+u for s in physical, u in physical, a in rightA])
    first = reshape(matrix_product(X, Y, rows, leftA, cols; field), chi, lb, d, mid, r)
    X2 = reshape(permutedims(first, (1, 3, 5, 2, 4)), chi*d*r, lb*mid)
    Y2 = reshape(B, lb*mid, di*rb)
    rows2 = vec([p+s-a for p in prefix, s in physical, a in rightA])
    inner2 = vec([b+u for b in leftB, u in physical])
    cols2 = vec([b+t for t in physical, b in rightB])
    second = reshape(matrix_product(X2, Y2, rows2, inner2, cols2; field), chi, d, r, di, rb)
    return permutedims(second, (1, 2, 4, 3, 5))
end

product_labels(prefix, rightA, rightB, physical) = (
    vec([p+s-t for p in prefix, s in physical, t in physical]),
    vec([a+b for a in rightA, b in rightB]),
)

"Block SVD with one global rank budget, rather than a separate budget per sector."
function truncated_svd(M, rows, cols; maxdim, svd_function = svd)
    @assert maxdim >= 1
    check_support(M, rows, cols)
    rg, cg = groups(rows), groups(cols)
    decompositions = Dict(
        q => svd_function(M[rg[q], cg[q]]) for q in sort(collect(intersect(keys(rg), keys(cg))))
    )
    choices = [
        (value = S, charge = q, index = k) for q in sort(collect(keys(decompositions))) for
        (k, S) in enumerate(decompositions[q].S)
    ]
    sort!(choices; by = x -> (-x.value, x.charge, x.index))
    if isempty(choices)
        # A zero matrix with disjoint row/column charges still has a rank-one
        # zero factorization. Choose a legal left vector; the bridge is zero.
        U = zeros(Float64, size(M, 1), 1)
        U[1, 1] = 1
        return (; U, S = [0.0], Vt = zeros(Float64, 1, size(M, 2)), charges = [first(rows)])
    end
    keep = choices[1:min(maxdim, length(choices))]
    U, Vt = zeros(Float64, size(M, 1), length(keep)), zeros(Float64, length(keep), size(M, 2))
    for (j, item) in enumerate(keep)
        q, k = item.charge, item.index
        U[rg[q], j] = decompositions[q].U[:, k]
        Vt[j, cg[q]] = decompositions[q].Vt[k, :]
    end
    return (; U, S = [x.value for x in keep], Vt, charges = [x.charge for x in keep])
end

"Least-squares bridge fit constrained to the same exact sectors as the rounded left factor."
function fit_bridge(U, M, rows, middle, cols)
    check_support(U, rows, middle)
    check_support(M, rows, cols)
    rg, mg, cg = groups(rows), groups(middle), groups(cols)
    R = zeros(Float64, length(middle), length(cols))
    for q in sort(collect(intersect(keys(rg), keys(mg), keys(cg))))
        R[mg[q], cg[q]] = U[rg[q], mg[q]] \ M[rg[q], cg[q]]
    end
    return R
end

"Schmidt balancing within sectors, including a regularized full right nullspace."
function schmidt_gauges(cores, physical; floor_ratio = 1e-18, svd_function = svd)
    @assert 0 < floor_ratio < 1
    charges = infer_charges(cores, physical)
    gauges = [ones(1, 1)]
    prefix, prefix_charges = ones(1, 1), [0]
    for (j, core) in enumerate(cores[1:(end-1)])
        l, d, di, r = size(core)
        columns = reshape(prefix*reshape(core, l, d*di*r), size(prefix, 1)*d*di, r)
        rows = vec([p+s-t for p in prefix_charges, s in physical, t in physical])
        check_support(columns, rows, charges[j+1])
        rg, cg = groups(rows), groups(charges[j+1])
        sectors = sort(collect(keys(cg)))
        blocks = Dict(q => Matrix(qr(columns[get(rg, q, Int[]), cg[q]]).R) for q in sectors)
        decompositions = Dict(q => svd_function(blocks[q]; full = true) for q in sectors)
        largest = maximum((maximum(s.S; init = 0.0) for s in values(decompositions)); init = 0.0)
        reference = iszero(largest) ? 1.0 : largest
        gauge = zeros(r, r)
        next = zeros(sum(size(b, 1) for b in values(blocks)), r)
        next_charges, offset = Int[], 0
        for q in sectors
            block, decomposition, indices = blocks[q], decompositions[q], cg[q]
            singular = vcat(decomposition.S, zeros(length(indices)-length(decomposition.S)))
            scales = sqrt.(max.(singular, floor_ratio*reference))
            gauge[indices, indices] = Diagonal(scales)*decomposition.Vt
            nr = size(block, 1)
            next[(offset+1):(offset+nr), indices] = block
            append!(next_charges, fill(q, nr))
            offset += nr
        end
        push!(gauges, gauge)
        prefix, prefix_charges = next, next_charges
    end
    push!(gauges, ones(1, 1))
    return gauges
end

end
