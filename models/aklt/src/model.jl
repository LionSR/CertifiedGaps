"""Periodic honeycomb quotient with two spin-3/2 sites per primitive cell.

A(x,y) joins B(x,y), B(x-1,y), B(x,y-1). Both coordinates are periodic.
Periods >=2 ensure three distinct neighbors and no duplicate bonds.
Small quotients may have short winding cycles; they are not open patches.
The site ordering snakes along rows; this does not change the Hamiltonian.
"""
function honeycomb(nx::Int, ny::Int)
    nx >= 2 && ny >= 2 || throw(ArgumentError("both periods must be at least two"))
    coordinates = Tuple{Int,Int,Int}[]
    for y in 0:(ny-1), x in (iseven(y) ? (0:(nx-1)) : ((nx-1):-1:0)), a in 0:1
        push!(coordinates, (x, y, a))
    end
    lookup = Dict(c=>j for (j, c) in enumerate(coordinates))
    groups = [Tuple{Int,Int}[] for _ in 1:3]
    for y in 0:(ny-1), x in 0:(nx-1)
        a = lookup[(x, y, 0)]
        for (color, (dx, dy)) in enumerate(((0, 0), (-1, 0), (0, -1)))
            b = lookup[(mod(x+dx, nx), mod(y+dy, ny), 1)]
            push!(groups[color], minmax(a, b))
        end
    end
    bonds = vcat(groups...)
    nsites = 2nx*ny
    length(unique(bonds)) == 3nsites÷2 || error("duplicate honeycomb bonds")
    degrees = zeros(Int, nsites)
    for (a, b) in bonds
        degrees[a] += 1
        degrees[b] += 1
    end
    all(==(3), degrees) || error("honeycomb graph is not trivalent")
    return (; nx, ny, nsites, bonds, groups, coordinates)
end

"""Bond projector P(total spin=3), using maximal-spin Clebsch-Gordan coefficients.

Local digits a,b=0:3 label m=3/2-a. Linear index is 1+a+4b.
Each fixed a+b block is a rank-one projector. All entries are real.
"""
function bond_projector()
    p = zeros(16, 16)
    for a in 0:3, b in 0:3, c in 0:3, d in 0:3
        a+b == c+d || continue
        p[1+a+4b, 1+c+4d] =
            sqrt(binomial(3, a)*binomial(3, b)*binomial(3, c)*binomial(3, d)) / binomial(6, a+b)
    end
    return p
end

function ITensors.space(::SiteType"AKLT3/2"; conserve_qns = false)
    return conserve_qns ? [QN("Sz", m)=>1 for m in (3, 1, -1, -3)] : 4
end
ITensors.op(::OpName"Id", ::SiteType"AKLT3/2") = Matrix{Float64}(I, 4, 4)

function bond_gate(sites, a, b, tau)
    matrix = Matrix{Float64}(I, 16, 16) + expm1(-tau)*bond_projector()
    return ITensor(
        reshape(matrix, 4, 4, 4, 4),
        prime(sites[a]),
        prime(sites[b]),
        dag(sites[a]),
        dag(sites[b]),
    )
end
