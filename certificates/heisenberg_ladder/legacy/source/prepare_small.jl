# Cache the small rings for the two remaining parameters while long rings run.
include("ladder.jl")
for J in (2//1,3//2), L in (2,3,4,5,6,8)
    tag=replace(string(Float64(J)),"."=>"p")
    path=joinpath(@__DIR__,"results","spectrum_J$(tag)_L$(L).json")
    isfile(path) && continue
    s=full_spectrum(L,J); es=s.energies
    write_result(path,Dict("metadata"=>provenance(Dict("J"=>string(J),"L"=>L)),
        "J"=>string(J),"L"=>L,"dimension"=>s.dimension,"sectors"=>s.rows,
        "seconds"=>s.seconds,"energies"=>es,"gap"=>es[2]-es[1]))
end
