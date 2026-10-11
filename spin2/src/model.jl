"""Spin-2 operators in descending m = 2, 1, 0, -1, -2 order; hbar = 1 (`ThermalMPOGap.Spins`)."""
spin_matrices() = Spins.spin_matrices(2)

function ITensors.space(::SiteType"S=2"; conserve_qns = false)
    return conserve_qns ? [QN("Sz", 2m) => 1 for m in 2:-1:-2] : 5
end
ITensors.op(::OpName"Sz", ::SiteType"S=2") = spin_matrices().sz
ITensors.op(::OpName"S+", ::SiteType"S=2") = spin_matrices().sp
ITensors.op(::OpName"S-", ::SiteType"S=2") = spin_matrices().sm
ITensors.op(::OpName"Sx", ::SiteType"S=2") = spin_matrices().sx
ITensors.op(::OpName"Sy", ::SiteType"S=2") = spin_matrices().sy
ITensors.op(::OpName"Id", ::SiteType"S=2") = Matrix{Float64}(I, 5, 5)

"""Real G_alpha = i S_alpha on symmetric traceless Cartesian tensors.

The generators act as Q -> J_alpha Q - Q J_alpha. They have norm 2,
not the norm 1 of the spin-1 generators in the existing notes.
Implemented by `ThermalMPOGap.Spins.cartesian_generators(2)`.
"""
real_generators() = Spins.cartesian_generators(2)

function bonds(n; periodic = true)
    n >= 2 || throw(ArgumentError("at least two sites are required"))
    periodic && (n < 4 || isodd(n)) && throw(ArgumentError("periodic runs use even n >= 4"))
    pairs = [(j, j+1) for j in 1:(n-1)]
    periodic && push!(pairs, (n, 1))
    return pairs
end

function hamiltonian(sites; periodic = true)
    os = OpSum()
    for (j, k) in bonds(length(sites); periodic)
        os += 1.0, "Sz", j, "Sz", k
        os += 0.5, "S+", j, "S-", k
        os += 0.5, "S-", j, "S+", k
    end
    return MPO(os, sites)
end

function bond_operator(sites, j, k)
    return op("Sz", sites, j)*op("Sz", sites, k) +
           0.5*op("S+", sites, j)*op("S-", sites, k) +
           0.5*op("S-", sites, j)*op("S+", sites, k)
end
