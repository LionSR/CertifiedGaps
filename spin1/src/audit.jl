"SHA-256 of the environment files and of every `.jl` file under src/, scripts/ and test/."
implementation_fingerprint() = source_hashes(
    normpath(joinpath(@__DIR__, ".."));
    files = ("Project.toml", "Manifest.toml", ".JuliaFormatter.toml"),
)

"Audit stored Julia outputs; never describes reading a cache as fresh computation."
function audit_results_impl(directory; execution_mode)
    provenance = verify_source()
    records = Dict{String,Any}()
    checks = 0
    upstream_sectors = [
        parse_exact_json(line) for
        line in eachline(joinpath(SOURCE, "evidence/source1/thermal-run/sectors.jsonl"))
    ]
    for source in (1, 2)
        path = joinpath(directory, "source$source.json")
        data = read_exact_json(path)
        @assert data["status"] == "reproduced"
        expected_keys =
            source == 1 ? Set((n, g) for n in 4:12 for g in 0:1) :
            Set((n, g) for n in (6, 8, 10) for g in 0:2)
        actual_keys =
            Set((Int(integer(c["n"])), Int(integer(c["twist"]))) for c in data["certificates"])
        @assert actual_keys == expected_keys &&
                length(data["certificates"]) == length(expected_keys)
        for record in data["certificates"]
            n, g = Int(integer(record["n"])), Int(integer(record["twist"]))
            total = integer(record["total"])
            # Rebuild interval endpoints and compare full integers to frozen evidence.
            rebuilt = thermal_certificate(source, n, g, total)
            for (key, value) in rebuilt
                if value isa Integer
                    @assert integer(record[key]) == value
                else
                    @assert record[key] == value
                end
            end
            sectors = filter(x -> integer(x["n"]) == n && integer(x["twist"]) == g, data["sectors"])
            @assert length(sectors) == n + 1
            @assert Set(integer(s["magnetization"]) for s in sectors) == Set(big.(0:n))
            @assert sum(integer(s["total"]) for s in sectors) == total
            @assert sum(integer(s["multiplicity"]) for s in sectors) == 3^n
            if source == 1
                for sector in sectors
                    reference = only(
                        filter(
                            s ->
                                integer(s["n"]) == n &&
                                integer(s["tw"]) == g &&
                                integer(s["magnetization"]) == integer(sector["magnetization"]),
                            upstream_sectors,
                        ),
                    )
                    for (actual_key, reference_key) in (
                        ("total", "W"),
                        ("dimension", "dimension"),
                        ("representatives", "representatives"),
                        ("multiplicity", "column_multiplicity"),
                    )
                        @assert integer(sector[actual_key]) == integer(reference[reference_key])
                    end
                end
            end
            checks += 1
        end
        records["source$source"] = Dict(
            "sha256" => bytes2hex(sha256(read(path))),
            "certificates" => length(expected_keys),
            "original_wall_seconds" => data["wall_seconds"],
            "original_completed_utc" => data["completed_utc"],
        )
    end
    trials_path = joinpath(directory, "trials.json")
    trials = read_exact_json(trials_path)
    @assert trials["status"] == "reproduced" && length(trials["certificates"]) == 2
    @assert Set(integer(c["length"]) for c in trials["certificates"]) == Set(big.((72, 120)))
    for c in trials["certificates"]
        n, u, v = integer(c["length"]), integer(c["U"]), integer(c["V"])
        expected = read_exact_json(joinpath(SOURCE, "evidence/trials/certificate$n.json"))
        @assert u == integer(expected[n == 72 ? "one_bond_energy_numerator" : "U"])
        @assert v == integer(expected[n == 72 ? "norm" : "V"])
        @assert v > 0 && -700741v - 500000u > 0
        @assert integer(c["positive_margin"]) == -700741v - 500000u
        @assert integer(c["floor_1e12_U_over_V"]) == fld(big(10)^12 * u, v)
        @assert sum(integer(b["U"]) for b in c["blocks"]) == u
        @assert sum(integer(b["V"]) for b in c["blocks"]) == v
    end
    for name in ("trials", "scalars", "arithmetic", "itensors")
        path = joinpath(directory, "$name.json")
        data = read_exact_json(path)
        @assert data["status"] ==
                (name == "itensors" ? "numerical-crosscheck-passed" : "reproduced")
        records[name] = Dict(
            "sha256" => bytes2hex(sha256(read(path))),
            "original_completed_utc" => data["completed_utc"],
        )
    end
    # Recompute all inexpensive exact scalar and operator checks in the audit too.
    scalar_data = read_exact_json(joinpath(directory, "scalars.json"))
    for (name, rebuilt) in (
        ("source1", source1_initialization()),
        ("source2", source2_initialization()),
        ("bootstrap", bootstrap_checks()),
    )
        @assert scalar_data[name] == parse_exact_json(JSON3.write(rebuilt))
    end
    arithmetic = read_exact_json(joinpath(directory, "arithmetic.json"))
    for (name, rebuilt) in (
        ("matrices", matrix_and_symmetry_checks()),
        ("triad", triad_checks()),
        ("bounds", error_bound_checks()),
    )
        @assert arithmetic[name] == parse_exact_json(JSON3.write(rebuilt))
    end
    result = Dict(
        "status" => "passed",
        "execution_mode" => execution_mode,
        "source_commit" => UPSTREAM_COMMIT,
        "source_integrity" => provenance,
        "thermal_certificates" => checks,
        "trial_certificates" => 2,
        "components" => records,
        "implementation_sha256" => implementation_fingerprint(),
        "julia" => string(VERSION),
        "audit_utc" => string(now(UTC)),
        "scope" => "Periodic spin-1 finite certificates. Analytic transfer and induction arguments are not formalized here. ITensors results are numerical controls.",
    )
    println(
        "PASS aggregate audit: 27 thermal certificates, two integer MPS trials, two purity seeds and bootstrap constants",
    )
    return write_result(joinpath(directory, "verification.json"), result)
end

function audit_results(directory; execution_mode = "stored-output-audit")
    path = joinpath(directory, "verification.json")
    write_result(path, Dict("status" => "running", "execution_mode" => execution_mode))
    try
        return audit_results_impl(directory; execution_mode)
    catch exception
        write_result(
            path,
            Dict(
                "status" => "failed",
                "execution_mode" => execution_mode,
                "error" => sprint(showerror, exception),
                "audit_utc" => string(now(UTC)),
            ),
        )
        rethrow()
    end
end
