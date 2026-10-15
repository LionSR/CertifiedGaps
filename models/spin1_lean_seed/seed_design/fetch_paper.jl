using Downloads, SHA
const COMMIT="adc7f1241b42e322a6451854ab7e4b4c146bf78a"
const BASE="https://raw.githubusercontent.com/openai/math/$COMMIT/preprints/" *
    "The-periodic-spin-one-Haldane-gap-September-24-2026"
mkpath(joinpath(@__DIR__,"upstream"))
for name in ("initialization.tex","transfer.tex","thermal72.tex","trials.tex","bootstrap.tex")
    path=joinpath(@__DIR__,"upstream",name)
    Downloads.download("$BASE/build/$name",path)
    println(name," ",bytes2hex(sha256(read(path))))
end
