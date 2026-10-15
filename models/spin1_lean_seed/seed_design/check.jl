include("certify.jl")

function check_certificate(path=joinpath(@__DIR__,"certificate.json"))
    c=JSON3.read(read(path,String),Dict{String,Any})
    for row in c["designs"]
        inp=row["input"];seed=row["seed"];chain=row["chain"]
        eps=parse_rational(inp["eps"])
        @assert audit_stored(inp["source"],eps)==row["input_sha256"]
        moments,errors=exact_sector_data(inp["source"],inp["mode"],eps)
        @assert parse_rational.(inp["masses"])==[m[inp["nmax"]]+e for (m,e) in zip(moments,errors)]
        for (m,e,f) in zip(moments,errors,inp["filters"])
            isnothing(f) && continue
            q=parse_rational.(f["q"])
            @assert exact_filter(m,e,f["shift"],q)==f
        end
        @assert verify_trial(seed["trial"])==row["trial"]
        calculated=exact_seed(inp,seed["trial"],seed["k"],parse_rational(seed["r"]))
        @assert calculated==seed
        @assert parse_rational(chain["beta"])==
            parse_rational(inp["base_beta"])*parse_rational(seed["r"])
        @assert parse_rational(first(chain["s"]))==parse_rational(seed["p"])
        @assert parse_rational(first(chain["t"]))==parse_rational(seed["q"])
        @assert chain_valid(chain)
        ratio=frac(784)/parse_rational(chain["beta"])
        @assert parse_rational(chain["r"])^numerator(ratio)>frac(20)^denominator(ratio)
        @assert 0<parse_rational(chain["gap_strict_lower"])<
            first(log_enclosure(parse_rational(chain["r"])))/parse_rational(chain["beta"])
        println("PASS ",row["name"]," exact scalar inequalities and ChainCertificate.Valid")
    end
    true
end
if abspath(PROGRAM_FILE) == @__FILE__
    check_certificate(isempty(ARGS) ? joinpath(@__DIR__,"certificate.json") : only(ARGS))
end
