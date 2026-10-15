# Verify the inventory only. Lean remains the certificate checker.
using SHA
root=dirname(@__DIR__)
text=read(joinpath(root,"certificate.json"),String)
block=match(r"\"source_hashes\"\s*:\s*\{([^}]*)\}",text)
block===nothing && error("missing source hashes")
for m in eachmatch(r"\"([^\"]+)\"\s*:\s*\"([0-9a-f]{64})\"",block[1])
    bytes2hex(sha256(read(joinpath(root,m[1]))))==m[2] || error("source hash mismatch: $(m[1])")
end
println("Source hashes match; this does not establish the open physical interpretation.")
