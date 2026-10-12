# Compact, byte-exact storage of the complete L=12 enclosed spectra.
#
# certify.jl / certify_parallel.jl write ~30 MB JSON files whose SHA-256 is recorded in each
# certificate (`enclosure_hashes`). The sorted level list contains every sector eigenvalue
# repeated by its symmetry multiplicity, so it run-length encodes to ~1.8 MB. We store
#   <name>.skeleton.json : the original JSON text with the two large arrays replaced by markers,
#   <name>.levels.bin    : radius table, run energies, run radius indices and run counts.
# `enclosure_text` rebuilds the original bytes exactly, so the checker re-verifies the
# certificate's recorded hash instead of trusting the compaction.
using SHA

const SPECTRA_MAGIC = b"SPIN1LD1"
const ARRAY_INDENT = "        "

"Text of one array exactly as JSON3.pretty wrote it (one element per line)."
array_text(xs) = "[\n" * join((ARRAY_INDENT * string(x) for x in xs), ",\n") * "\n    ]"

"Split the original pretty JSON into a skeleton and the energies/radii arrays."
function split_enclosure(text)
    arrays = Dict{String,Vector{Float64}}()
    skeleton = text
    for key in ("energies", "radii")
        head = "\"$key\": "
        start = last(findfirst(head * "[\n", skeleton)) - 1        # index of '['
        stop = last(findnext("\n    ]", skeleton, start))          # index of ']'
        body = skeleton[start:stop]
        values = [parse(Float64, rstrip(strip(line), ',')) for line in
                  split(body, '\n')[2:end-1]]
        array_text(values) == body || error("array $key is not in canonical form")
        arrays[key] = values
        skeleton = skeleton[1:start-1] * "\"@$key@\"" * skeleton[stop+1:end]
    end
    return skeleton, arrays["energies"], arrays["radii"]
end

"Run-length encode the sorted (energy, radius) list against a table of distinct radii."
function write_levels(path, energies, radii)
    table = sort(unique(radii))
    length(table) <= 256 || error("too many distinct radii")
    index = Dict(r => UInt8(i - 1) for (i, r) in enumerate(table))
    runE = Float64[]; runR = UInt8[]; runC = UInt8[]
    for (e, r) in zip(energies, radii)
        if !isempty(runE) && runE[end] === e && table[runR[end]+1] === r && runC[end] < 0xff
            runC[end] += 0x01
        else
            push!(runE, e); push!(runR, index[r]); push!(runC, 0x01)
        end
    end
    open(path, "w") do io
        write(io, SPECTRA_MAGIC)
        write(io, htol(UInt64(length(runE))), htol(UInt64(length(table))))
        write(io, htol.(table)); write(io, htol.(runE)); write(io, runR); write(io, runC)
    end
    return path
end

"Inverse of `write_levels`: the full sorted energies and radii."
function read_levels(path)
    open(path) do io
        read(io, length(SPECTRA_MAGIC)) == SPECTRA_MAGIC || error("bad magic in $path")
        nrun = Int(ltoh(read(io, UInt64))); ntab = Int(ltoh(read(io, UInt64)))
        table = ltoh.(read!(io, Vector{Float64}(undef, ntab)))
        runE = ltoh.(read!(io, Vector{Float64}(undef, nrun)))
        runR = read!(io, Vector{UInt8}(undef, nrun))
        runC = read!(io, Vector{UInt8}(undef, nrun))
        eof(io) || error("trailing bytes in $path")
        energies = Float64[]; radii = Float64[]
        for i in 1:nrun, _ in 1:runC[i]
            push!(energies, runE[i]); push!(radii, table[runR[i]+1])
        end
        return energies, radii
    end
end

"Compact `results/<name>.json` into skeleton + levels files, checking an exact round trip."
function compact_enclosure(dir, name)
    text = read(joinpath(dir, name * ".json"), String)
    skeleton, energies, radii = split_enclosure(text)
    write(joinpath(dir, name * ".skeleton.json"), skeleton)
    write_levels(joinpath(dir, name * ".levels.bin"), energies, radii)
    enclosure_text(dir, name; prefer_full = false) == text ||
        error("round trip of $name is not byte-exact")
    return text_sha256(text)
end

# sha256(::String) hashes the string object id per block (quadratic); hash a byte copy.
"Hex SHA-256 of the UTF-8 bytes of a string."
text_sha256(text) = bytes2hex(sha256(Vector{UInt8}(text)))

"Original JSON text of an enclosure: the full file if present, else rebuilt from compact form."
function enclosure_text(dir, name; prefer_full = true)
    full = joinpath(dir, name * ".json")
    prefer_full && isfile(full) && return read(full, String)
    skeleton = read(joinpath(dir, name * ".skeleton.json"), String)
    energies, radii = read_levels(joinpath(dir, name * ".levels.bin"))
    # Splice by positions found in the small skeleton; the arrays are streamed into a buffer.
    pieces = sort([(first(findfirst("\"@$key@\"", skeleton)), key)
                   for key in ("energies", "radii")])
    io = IOBuffer(); cursor = 1
    for (position, key) in pieces
        write(io, SubString(skeleton, cursor, position - 1))
        write(io, array_text(key == "energies" ? energies : radii))
        cursor = position + length("\"@$key@\"")
    end
    write(io, SubString(skeleton, cursor))
    return String(take!(io))
end
