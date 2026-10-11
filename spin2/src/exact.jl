"""Independent Sz-sector Hamiltonian, assembled in the product basis.

Limited to n <= 6: full diagonalization is a reference calculation, not
the scalable tensor-network method. 'Exact' means untruncated Hilbert
space here, NOT interval-certified floating-point arithmetic.

Built by `ThermalMPOGap.ExactDiagonalization` from the two-site bond
Sz Sz + (S+ S- + S- S+)/2 on the pairs of `bonds(n; periodic)`; the local
charge is 2m, so the sector of magnetization M has charge 2M. Basis states
follow the core convention (site 1 most significant).
"""
function sector_hamiltonian(n, magnetization; periodic = true)
    n <= 6 || throw(ArgumentError("full diagonalization is limited to n <= 6"))
    chain = spin2_chain(n; periodic)
    states = ED.sector_states(chain, n, 2magnetization)
    isempty(states) && throw(ArgumentError("empty magnetization sector"))
    return Symmetric(Matrix(ED.sector_hamiltonian(chain, n, states)))
end

"Two-site Heisenberg bond Sz Sz + (S+ S- + S- S+)/2 on kron(site i, site j)."
function bond_matrix()
    s = spin_matrices()
    return kron(s.sz, s.sz)+(kron(s.sp, s.sm)+kron(s.sm, s.sp))/2
end

"The spin-2 Heisenberg chain on `bonds(n; periodic)`, with local charge 2m."
function spin2_chain(n; periodic = true)
    bonds(n; periodic) # validate the size
    h = bond_matrix()
    edges = periodic ? ED.ring_bonds(h) : ED.open_bonds(h)
    return ED.Chain(5, edges; charge = [4, 2, 0, -2, -4])
end

"Full sorted spectrum with multiplicities; sectors M < 0 are the spin flips of M > 0."
function exact_spectrum(n; periodic = true)
    n <= 6 || throw(ArgumentError("full diagonalization is limited to n <= 6"))
    return ED.spectrum(spin2_chain(n; periodic), n; spin_flip = true, dense_max = 5^6).E
end
