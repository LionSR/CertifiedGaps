"""
    CertifiedGaps.RingModels

A translation-invariant ring H_L = Σ_x h_{x,x+1} + Σ_x o_x given by exact matrices, and the
exact reflection-positivity certificate −h = Σ_α w_α C_α ⊗ R conj(C_α) R (w_α > 0, o real) that
the criterion needs. A c-site unit cell is one site of dimension dᶜ, with the terms inside the
cell in `onsite`. Bond matrices act on `kron(site x, site x+1)`.

Everything here is exact (`Rational{BigInt}` entries, possibly complex); the floating-point
`ReflectionPositivity.rp_test` is only a screening step.
"""
module RingModels

using LinearAlgebra
using ..ExactArithmetic: RationalBig, parse_rational
using ..ReflectionPositivity: realign

export RingModel, cell_model, isreal_model, conserves_charge, RPCertificate, rp_certificate,
    verify_rp, exact_matrix, encode_matrix, decode_matrix, encode_model, decode_model, encode_rp,
    decode_rp

const ExactEntry = Complex{RationalBig}

"Convert a matrix of integers, rationals or exact complex rationals to `Complex{RationalBig}`."
exact_matrix(A::AbstractMatrix) = map(_exact, A)
_exact(x::Union{Integer,Rational}) = ExactEntry(RationalBig(x))
_exact(x::Complex{<:Union{Integer,Rational}}) =
    ExactEntry(RationalBig(real(x)), RationalBig(imag(x)))
_exact(x) = throw(ArgumentError("refusing inexact matrix entry $(repr(x)) of type $(typeof(x))"))

"""
    RingModel(d, bond, onsite = zeros(d, d); charge = nothing, modulus = 0, flip = nothing,
        name = "")

Ring with local dimension `d`, Hermitian bond `bond` (d²×d², on `kron(x, x+1)`) and Hermitian
on-site term `onsite`, all with exact entries. `charge[a]` is an optional additive quantum
number of local state `a`, conserved by H (checked exactly), summed modulo `modulus`
(0 for a U(1) charge); exact diagonalisation splits the ring into its sectors.

`flip` is an optional exact local unitary F that reverses the charge (F[a, b] ≠ 0 only when
`charge[a] ≡ −charge[b]`) and leaves h and o invariant, such as a spin flip. Then F^{⊗L} commutes
with H and with translations and maps the sector q onto −q, so the two have the same spectrum
and only one of them is diagonalised. All conditions are checked exactly.
"""
struct RingModel
    name::String
    d::Int
    bond::Matrix{ExactEntry}
    onsite::Matrix{ExactEntry}
    charge::Vector{Int}
    modulus::Int
    flip::Union{Nothing,Matrix{ExactEntry}}
    function RingModel(d::Integer, bond::AbstractMatrix, onsite::AbstractMatrix = zeros(Int, d, d);
        charge = nothing, modulus::Integer = 0, flip = nothing, name = "")
        B, O = exact_matrix(bond), exact_matrix(onsite)
        size(B) == (d^2, d^2) || throw(DimensionMismatch("bond must be $(d^2)×$(d^2)"))
        size(O) == (d, d) || throw(DimensionMismatch("onsite must be $d×$d"))
        B == B' || throw(ArgumentError("bond is not Hermitian"))
        O == O' || throw(ArgumentError("onsite term is not Hermitian"))
        q = charge === nothing ? zeros(Int, d) : collect(Int, charge)
        length(q) == d || throw(DimensionMismatch("charge must have length $d"))
        F = flip === nothing ? nothing : exact_matrix(flip)
        model = new(String(name), Int(d), B, O, q, Int(modulus), F)
        conserves_charge(model) || throw(ArgumentError("H does not conserve the given charge"))
        F === nothing || _valid_flip(model) ||
            throw(ArgumentError("flip is not a charge-reversing local symmetry of H"))
        return model
    end
end

"""
    cell_model(model, c; charge = nothing, kwargs...) -> RingModel

The same Hamiltonian as a ring of c-site cells, built exactly: the cell is `kron` of its c sites,
the on-site term holds the c on-site terms and the c − 1 bonds inside the cell, and the bond
joins the last site of a cell to the first site of the next, so the ring of L cells is the ring
of cL sites. `charge` (length dᶜ) and the other keywords (`modulus`, `flip`, …) are those of
`RingModel` and are checked exactly for the cell model. A cell may conserve a charge that the
site model does not (a sublattice charge, say).
"""
function cell_model(m::RingModel, c::Integer; charge = nothing, kwargs...)
    c >= 1 || throw(ArgumentError("a cell has at least one site"))
    d = m.d
    id(k) = Matrix{ExactEntry}(I, d^k, d^k)
    onsite = sum(kron(id(k), m.onsite, id(c - 1 - k)) for k in 0:(c-1)) +
             sum((kron(id(k), m.bond, id(c - 2 - k)) for k in 0:(c-2));
                 init = zeros(ExactEntry, d^c, d^c))
    return RingModel(d^c, kron(id(c - 1), m.bond, id(c - 1)), onsite; charge, kwargs...,
        name = "$(m.name), $c-site cells")
end

"True when every entry of the bond and the on-site term is real."
isreal_model(m::RingModel) = all(iszero ∘ imag, m.bond) && all(iszero ∘ imag, m.onsite)

_reduce(m::RingModel, q) = m.modulus == 0 ? q : mod(q, m.modulus)

# Exact check of the conditions on `flip` stated in the `RingModel` docstring.
function _valid_flip(m::RingModel)
    F, d = m.flip, m.d
    size(F) == (d, d) && F * F' == I || return false
    all(iszero(F[a, b]) || _reduce(m, m.charge[a] + m.charge[b]) == 0 for a in 1:d, b in 1:d) ||
        return false
    FF = kron(F, F)
    return FF * m.bond * FF' == m.bond && F * m.onsite * F' == m.onsite
end

"Exact check that every nonzero matrix element of h and o preserves the total charge."
function conserves_charge(m::RingModel)
    d, q = m.d, m.charge
    for j in 1:d^2, i in 1:d^2
        iszero(m.bond[i, j]) && continue
        (a, b), (c, e) = divrem(i - 1, d), divrem(j - 1, d)
        _reduce(m, q[a+1] + q[b+1] - q[c+1] - q[e+1]) == 0 || return false
    end
    return all(iszero(m.onsite[i, j]) || _reduce(m, q[i] - q[j]) == 0 for i in 1:d, j in 1:d)
end

"""
    RPCertificate(weights, factors, reflection = nothing)

Exact reflection-positivity data −h = Σ_α w_α C_α ⊗ R conj(C_α) R with rationals w_α > 0, where
R is `reflection` (a real symmetric involution commuting with the on-site term) or the
identity. C_α acts on the left site and R conj(C_α) R on the right one. Three forms, as in the
Lean ring model: "real" (R = I, real C_α: −h = Σ w C ⊗ C), "conj" (R = I, some C_α complex)
and "reflected" (R given, e.g. the swap of the two sites of a cell).
"""
struct RPCertificate
    weights::Vector{RationalBig}
    factors::Vector{Matrix{ExactEntry}}
    reflection::Union{Nothing,Matrix{ExactEntry}}
end
RPCertificate(weights, factors) = RPCertificate(weights, factors, nothing)

"""
    rp_certificate(model; factors = nothing, weights = nothing, reflection = nothing)

Exact reflection-positivity certificate for the bond of `model`, in the form of
`RPCertificate`. Given `factors` (and optional `weights`, default 1), they are checked with
`verify_rp`. Otherwise the Hermitian realignment of −(1 ⊗ R) h (1 ⊗ R) is factored exactly as
Σ_k p_k u_k u_k† by rational LDL† with diagonal pivoting, which succeeds iff it is positive
semidefinite; then C_k is u_k reshaped row-major. Throws if
the bond is not reflection positive in this basis or the on-site term is not real.
"""
function rp_certificate(model::RingModel; factors = nothing, weights = nothing,
    reflection = nothing)
    R = reflection === nothing ? nothing : exact_matrix(reflection)
    if factors !== nothing
        C = [exact_matrix(F) for F in factors]
        w = weights === nothing ? ones(RationalBig, length(C)) : RationalBig.(weights)
        cert = RPCertificate(w, C, R)
    else
        cert = _ldl_certificate(model, R)
    end
    verify_rp(model, cert) || throw(ArgumentError("the RP certificate does not verify"))
    return cert
end

function _ldl_certificate(model::RingModel, reflection)
    d = model.d
    # −h = Σ w C ⊗ R conj(C) R  iff  −(1 ⊗ R) h (1 ⊗ R) = Σ w C ⊗ conj(C).
    flip = reflection === nothing ? I : kron(Matrix{ExactEntry}(I, d, d), reflection)
    R = realign(-(flip * model.bond * flip), d)
    R == R' || throw(ArgumentError("R[−h] is not Hermitian: the bond is not RP in this basis"))
    R = copy(R)
    weights, factors = RationalBig[], Matrix{ExactEntry}[]
    for k in 1:d^2
        p = real(R[k, k])
        if iszero(p)
            # A PSD matrix with a zero diagonal entry has a zero row and column there.
            all(iszero, view(R, :, k)) ||
                throw(ArgumentError("R[−h] is not positive semidefinite"))
            continue
        end
        p > 0 || throw(ArgumentError("R[−h] is not positive semidefinite"))
        u = R[:, k] / p
        R .-= p .* (u * u')
        push!(weights, p)
        # R[A ⊗ conj(A)] = vec_row(A) vec_row(A)†, with vec_row the row-major vectorisation.
        push!(factors, permutedims(reshape(u, d, d)))
    end
    return RPCertificate(weights, factors, reflection)
end

"""
    verify_rp(model, cert) -> Bool

Exact identity −h = Σ_α w_α C_α ⊗ R conj(C_α) R, positivity of the weights, reality of the
on-site term o, and for a reflection R: R real, symmetric, R² = 1 and Ro = oR. All exact.
"""
function verify_rp(model::RingModel, cert::RPCertificate)
    all(>(0), cert.weights) || return false
    length(cert.weights) == length(cert.factors) || return false
    all(iszero ∘ imag, model.onsite) || return false
    R = cert.reflection === nothing ? Matrix{ExactEntry}(I, model.d, model.d) : cert.reflection
    size(R) == (model.d, model.d) || return false
    all(iszero ∘ imag, R) && R == transpose(R) && R * R == I || return false
    R * model.onsite == model.onsite * R || return false
    d2 = model.d^2
    total = zeros(ExactEntry, d2, d2)
    for (w, C) in zip(cert.weights, cert.factors)
        size(C) == (model.d, model.d) || return false
        total .+= w .* kron(C, R * conj.(C) * R)
    end
    return total == -model.bond
end

# JSON encoding: every exact number is a string "num//den"; a complex matrix is {re, im}.
_token(x::RationalBig) = denominator(x) == 1 ? string(numerator(x)) :
                         string(numerator(x), "//", denominator(x))

"Encode an exact matrix as `{\"re\": rows, \"im\": rows}` (`im` omitted when zero)."
function encode_matrix(A::AbstractMatrix)
    E = exact_matrix(A)
    rows(f) = [[_token(f(E[i, j])) for j in axes(E, 2)] for i in axes(E, 1)]
    record = Dict{String,Any}("re" => rows(real))
    all(iszero ∘ imag, E) || (record["im"] = rows(imag))
    return record
end

"Inverse of `encode_matrix`; refuses any inexact entry."
function decode_matrix(record)
    rows(key) = permutedims(reduce(hcat, [[parse_rational(x) for x in row] for row in record[key]]))
    re = rows("re")
    im = haskey(record, "im") ? rows("im") : zero(re)
    return ExactEntry.(re, im)
end

"JSON record of a model (exact strings); `charge_flip` only when the model has one."
function encode_model(m::RingModel)
    record = Dict{String,Any}("name" => m.name, "d" => m.d, "bond" => encode_matrix(m.bond),
        "onsite" => encode_matrix(m.onsite), "charge" => m.charge, "charge_modulus" => m.modulus)
    m.flip === nothing || (record["charge_flip"] = encode_matrix(m.flip))
    return record
end

"Rebuild a `RingModel` from `encode_model` output, re-running all exact checks."
decode_model(r) = RingModel(_int(r["d"]), decode_matrix(r["bond"]), decode_matrix(r["onsite"]);
    charge = _int.(r["charge"]), modulus = _int(r["charge_modulus"]),
    flip = haskey(r, "charge_flip") ? decode_matrix(r["charge_flip"]) : nothing,
    name = String(r["name"]))

# `parse_exact_json` reads every integer as a string.
_int(x::AbstractString) = parse(Int, x)
_int(x::Integer) = Int(x)

"""
    encode_rp(cert) -> Dict

JSON record of an RP certificate. Each factor carries its weight and whether it is real and
symmetric or antisymmetric, the case covered by the Lean ring model.
"""
function encode_rp(cert::RPCertificate)
    function factor(w, C)
        kind = C == transpose(C) ? "symmetric" : C == -transpose(C) ? "antisymmetric" : "neither"
        Dict{String,Any}("weight" => _token(w), "matrix" => encode_matrix(C),
            "real" => all(iszero ∘ imag, C), "transpose_symmetry" => kind)
    end
    factors = [factor(w, C) for (w, C) in zip(cert.weights, cert.factors)]
    if cert.reflection !== nothing
        return Dict{String,Any}("form" => "reflected", "factors" => factors,
            "reflection" => encode_matrix(cert.reflection),
            "identity" => "-h = sum_a weight_a * kron(C_a, R * conj(C_a) * R)")
    end
    form = all(f["real"] for f in factors) ? "real" : "conj"
    return Dict{String,Any}("form" => form, "factors" => factors,
        "identity" => "-h = sum_a weight_a * kron(C_a, conj(C_a))")
end

"Inverse of `encode_rp`."
decode_rp(r) = RPCertificate([parse_rational(f["weight"]) for f in r["factors"]],
    [decode_matrix(f["matrix"]) for f in r["factors"]],
    haskey(r, "reflection") ? decode_matrix(r["reflection"]) : nothing)

end
