using CertifiedGaps.Results, LinearAlgebra, JSON3
BLAS.set_num_threads(1)
root=@__DIR__
record=joinpath(root,"results/archive.json")
rows=isfile(record) ? Any[JSON3.read(read(record,String)).rows...] : Any[]
for name in readdir(joinpath(root,"results"))
    startswith(name,"spectrum_L") && endswith(name,".json") || continue
    path=joinpath(root,"results",name)
    filesize(path)>1_000_000 || continue
    original=sha256_file(path);original_bytes=filesize(path)
    output=path*".gz"
    run(pipeline(`gzip -n -c -- $path`,stdout=output*".tmp"))
    @assert read(`gzip -dc -- $(output*".tmp")`)==read(path)
    mv(output*".tmp",output;force=true)
    # Keep the original locally, outside versioned deliverables.
    mkpath(joinpath(root,"runtime/archive"))
    mv(path,joinpath(root,"runtime/archive",name);force=true)
    push!(rows,(;file=name*".gz",original_sha256=original,
        compressed_sha256=sha256_file(output),original_bytes,compressed_bytes=filesize(output)))
end
write_result(joinpath(root,"results/archive.json"),(;
    metadata=run_metadata((;method="lossless gzip, original JSON written by Results.write_result");
        source_sha256=source_hashes(root;dirs=(),files=("archive_spectra.jl","result_io.jl")),
        status="lossless storage conversion"),rows))
