using JSON3

"Read a result JSON, including losslessly gzip-compressed complete spectra."
function read_result_json(path)
    raw=endswith(path,".gz") ? read(`gzip -dc -- $path`,String) : read(path,String)
    JSON3.read(raw)
end

using CertifiedGaps.ExactDiagonalization: PartialSpectrum
import CertifiedGaps.ExactDiagonalization: logZ_bounds

"A complete numerical spectrum represented by levels and exact integer multiplicities."
struct WeightedSpectrum
    E::Vector{Float64}
    weights::Vector{Int}
    N::Int
    exact::Bool
end
function logZ_bounds(s::WeightedSpectrum,beta::Real)
    e0=s.E[1]
    z=sum(w*exp(-beta*(e-e0)) for (e,w) in zip(s.E,s.weights))
    value=log(z)-beta*e0
    value,value
end
function result_spectrum(result,t)
    row=result.spectra[t]
    e=collect(Float64,row.energies)
    if haskey(row,:weights)
        weights=collect(Int,row.weights)
        @assert length(e)==length(weights)
        @assert issorted(e) && all(>(0),weights)
        @assert sum(weights)==3^result.L
        return WeightedSpectrum(e,weights,sum(weights),true)
    end
    PartialSpectrum(e,3^result.L,true)
end
spectrum_gap(s::PartialSpectrum)=s.E[2]-s.E[1]
spectrum_gap(s::WeightedSpectrum)=s.weights[1]>1 ? 0.0 : s.E[2]-s.E[1]
spectrum_moment(s::PartialSpectrum,p)=sum(e->e^p,s.E)/s.N
spectrum_moment(s::WeightedSpectrum,p)=
    sum(w*e^p for (e,w) in zip(s.E,s.weights))/s.N
