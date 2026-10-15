# Copy the original Lean sources into this project's writable build tree. Never build in lean/.
using SHA, JSON
const HERE=@__DIR__
const SOURCE=normpath(joinpath(HERE,"../../lean"))
const VENDOR=joinpath(HERE,"vendor")
mkpath(VENDOR)
hashes=Dict{String,String}()
for (dir, _, files) in walkdir(joinpath(SOURCE,"CertifiedGaps")), f in files
    endswith(f,".lean") || continue
    p=joinpath(dir,f); rel=relpath(p,SOURCE); dst=joinpath(VENDOR,rel)
    mkpath(dirname(dst)); cp(p,dst;force=true)
    hashes[rel]=bytes2hex(sha256(read(p)))
end
cp(joinpath(SOURCE,"CertifiedGaps.lean"),joinpath(VENDOR,"CertifiedGaps.lean");force=true)
# ProductLimit imports all of Mathlib in the original. Narrow ONLY the copied import header
# to avoid compiling all 8,943 cache modules when the cache host is blocked. Proof bodies and
# definitions remain byte-for-byte identical. This change is not applied to the repository's lean/.
p=joinpath(VENDOR,"CertifiedGaps/Transfer/ProductLimit.lean")
s=read(p,String)
imports="""
import Mathlib.Analysis.SpecialFunctions.Exponential
import Mathlib.Analysis.Calculus.Deriv.Mul
import Mathlib.Analysis.Calculus.Deriv.Comp
import Mathlib.Analysis.Calculus.Deriv.Slope
import Mathlib.Analysis.Normed.Algebra.MatrixExponential
import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Tactic.Abel
import Mathlib.Tactic.NoncommRing
import Mathlib.Tactic.GCongr
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.FinCases
"""
@assert count("import Mathlib\n",s)==1
write(p,replace(s,"import Mathlib\n"=>imports))
@assert split(read(p,String),"/-!";limit=2)[2]==split(s,"/-!";limit=2)[2]
open(joinpath(HERE,"source_hashes.json"),"w") do io
    JSON.print(io,Dict("lean_sources"=>hashes,"import_only_patch"=>"ProductLimit.lean",
        "toolchain"=>strip(read(joinpath(SOURCE,"lean-toolchain"),String))),2)
end
println("Copied original Lean source; only ProductLimit's local import list was narrowed.")
