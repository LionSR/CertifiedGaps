# Spectra of stored results. Reading (plain or gzip-compressed JSON) is
# `CertifiedGaps.Results.read_result(path; exact=false)`; complete multiplet-weighted spectra
# and their numerical log Z are the core's `WeightedSpectrum` and `logZ_bounds`.
using CertifiedGaps.ExactDiagonalization: PartialSpectrum, WeightedSpectrum, logZ_bounds
"Spectrum at angle index `t` of a complete-spectrum result (weighted if SU(2)-resolved)."
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
"First excitation gap; zero if the lowest level is a degenerate multiplet."
spectrum_gap(s::PartialSpectrum)=s.E[2]-s.E[1]
spectrum_gap(s::WeightedSpectrum)=s.weights[1]>1 ? 0.0 : s.E[2]-s.E[1]
"`p`-th moment of the energy over the complete spectrum (infinite temperature)."
spectrum_moment(s::PartialSpectrum,p)=sum(e->e^p,s.E)/s.N
spectrum_moment(s::WeightedSpectrum,p)=
    sum(w*e^p for (e,w) in zip(s.E,s.weights))/s.N
