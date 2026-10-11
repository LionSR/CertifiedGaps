"""Bond projector P(total spin=3) from its polynomial in S.S, independently of the CG formula.

The eigenvalues of S.S on two spin-3/2 sites are -15/4, -11/4, -3/4 and 9/4 (total spin
0, 1, 2, 3), so P = (x+15/4)(x+11/4)(x+3/4)/90 with x = S.S. P is invariant under exchange of
the two sites, so it is the same on `kron(site i, site j)` and `kron(site j, site i)`.
"""
function polynomial_projector()
    s = spin_matrices(3//2)
    x = kron(s.sz, s.sz)+(kron(s.sp, s.sm)+kron(s.sm, s.sp))/2
    return (x+15I/4)*(x+11I/4)*(x+3I/4)/90
end

"""Full spectrum by sector-wise dense diagonalization in every integer 2Sz sector; n<=8.

Uses the core product-basis assembly with the polynomial projector on every bond. The dense
eigensolve is untruncated but uses Float64, not interval arithmetic. Sectors with 2Sz<0 are
counted through spin-flip symmetry.
"""
function full_spectrum(n, bonds; progress = true, sector_callback = nothing)
    n <= 8 || throw(ArgumentError("reference diagonalization is limited to eight sites"))
    chain = Chain(4, graph_bonds(bonds, polynomial_projector()); charge = [3, 1, -1, -3])
    energies = Float64[]
    sectors = Dict{String,Any}[]
    function on_sector(sector)
        for _ in 1:sector.multiplicity
            append!(energies, sector.levels)
        end
        row = Dict(
            "twice_sz"=>sector.charge,
            "dimension"=>sector.dimension,
            "spin_flip_multiplicity"=>sector.multiplicity,
            "minimum_energy"=>first(sector.levels),
            "seconds"=>sector.seconds,
        )
        push!(sectors, row)
        progress && println(
            "ED n=$n 2Sz=$(sector.charge) dim=$(sector.dimension) " *
            "E0=$(first(sector.levels)) seconds=$(sector.seconds)",
        )
        flush(stdout)
        !isnothing(sector_callback) && sector_callback(energies, sectors)
    end
    spec = spectrum(chain, n; spin_flip = true, dense_max = 4^n, on_sector)
    spec.exact && spec.N==4^n && length(spec.E)==4^n || error("incomplete full spectrum")
    return (; energies = spec.E, sectors)
end
