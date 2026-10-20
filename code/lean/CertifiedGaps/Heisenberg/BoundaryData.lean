/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.WordBand

/-!
# The integer recurrence of the boundary-field thermal appendix

The appendix "Certified short-chain thermal moments" of OpenAI's *A boundary-field gap for the
spin-one Heisenberg chain* (2026) computes fifteen integers `k_{n,g}` (eq.
(appth:integer-table)) by an exact integer Horner recurrence (eqs. (appth:coefficients),
(appth:scalar), (appth:quotient)). This file implements that recurrence as computable
definitions, in the packed form of `Trace/PackedHorner.lean`, and states the hypothesis that it
reproduces the printed integers (`BoundaryTotals`).

## The recurrence

Words `w ∈ {-1, 0, 1}ⁿ` are numbered `k = ∑ⱼ (wⱼ + 1) 3ʲ` (site `0` least significant), so that
`wⱼ = digit k j - 1`. With `p_w = (∑ⱼ wⱼ) mod 2` (`par`) and
`C_n = 3⌊(n-1)/2⌋ + 2((n-1) mod 2) + 6/5` (`c10 n = 10 C_n`), the integer matrix
`D = 80 J A J⁻¹` has

* the diagonal `η(w) = 80 - 10 C_n - 10 ∑ⱼ wⱼ wⱼ₊₁ + 6 w₀ - 6 g w_{n-1}` (`eta`);
* the entries `-10` of the bond hops `w ↔ w + σ(eⱼ - eⱼ₊₁)` (`hopBand n j (j + 1) 10`,
  `hopBand n (j + 1) j 10` of `Heisenberg/WordBand.lean`; in the numbering, offsets `∓ 2·3ʲ`);
* for `g = 0`, the entries `-3(1 + p_w)` of the end hops `w ↔ w + σ e_{n-1}` (`endUp`,
  `endDown`; offsets `± 3^{n-1}`).

The coefficients are `tᵢ = ⌊Q sᵢ / S_K⌋` with `sᵢ = 49ⁱ K!/i!`, `K = 192`, `Q = 2⁵⁰`
(`coeffs`), and the recurrence `𝒳ᵢ = ⌊D 𝒳ᵢ₊₁ / 80⌋ + tᵢ I` from `𝒳_{K+1} = 0` runs on all `3ⁿ`
columns at once (`runBoundary`, through `Trace.Packed.run`: lanes of `128` bits with bias
`2⁵²`). The table entry is `k_{n,g} = ⌊10¹⁶ ∑_{w,r} x(w,r)² (2 - p_w)(1 + p_r) / (2 Q²)⌋`
(`quotient`).
-/

namespace CertifiedGaps.SpinOneHeisenberg

namespace BoundaryThermal

open Trace.Packed

/-- The chain lengths of the table. -/
def lengths : List ℕ := [2, 3, 4, 5, 6]

/-- The right fields of the table. -/
def fields : List ℤ := [1, -1, 0]

/-- `10 C_n = 30⌊(n-1)/2⌋ + 20((n-1) mod 2) + 12`. -/
def c10 (n : ℕ) : ℕ := 30 * ((n - 1) / 2) + 20 * ((n - 1) % 2) + 12

/-- The parity `p_w = (∑ⱼ wⱼ) mod 2 = (∑ⱼ digit k j + n) mod 2` of the word `k` of length
`n`. -/
def par (n k : ℕ) : ℕ := (((List.range n).map (digit k)).sum + n) % 2

/-- The diagonal entry `η(w)` of `D` (eq. (appth:diagonal)). -/
def eta (n : ℕ) (g : ℤ) (k : ℕ) : ℤ :=
  80 - (c10 n : ℤ) - 10 * ((List.range (n - 1)).map fun j ↦ letter k j * letter k (j + 1)).sum +
    6 * letter k 0 - (if g = 0 then 0 else 6 * g * letter k (n - 1))

/-- The end hop `w_{n-1} ↦ w_{n-1} + 1`: column `k + 3^{n-1}`, entry `-3(1 + p_w)`. -/
def endUp (n : ℕ) : Hop :=
  ⟨true, 3 ^ (n - 1), (List.range (3 ^ n)).map fun k ↦
    if digit k (n - 1) ≤ 1 then 3 * (1 + par n k) else 0⟩

/-- The end hop `w_{n-1} ↦ w_{n-1} - 1`: column `k - 3^{n-1}`, entry `-3(1 + p_w)`. -/
def endDown (n : ℕ) : Hop :=
  ⟨false, 3 ^ (n - 1), (List.range (3 ^ n)).map fun k ↦
    if 1 ≤ digit k (n - 1) then 3 * (1 + par n k) else 0⟩

/-- The off-diagonal bands of `D`: the two hops of every bond with entry `-10`, and the end
hops. -/
def hops (n : ℕ) (g : ℤ) : List Hop :=
  (List.range (n - 1)).flatMap (fun j ↦ [hopBand n j (j + 1) 10, hopBand n (j + 1) j 10]) ++
    (if g = 0 then [endUp n, endDown n] else [])

/-- The integer matrix `D = 80 J A J⁻¹` of eq. (appth:diagonal) as a banded matrix. -/
def band (n : ℕ) (g : ℤ) : Band where
  pos := (List.range (3 ^ n)).map fun k ↦ (eta n g k).toNat
  neg := (List.range (3 ^ n)).map fun k ↦ (-eta n g k).toNat
  hops := hops n g

/-- The packing parameters: lanes of `128` bits, bias `2⁵²`, scale `80`, and the magic
division `⌊u / 80⌋ = ⌊u ⌈2⁶⁸/80⌉ / 2⁶⁸⌋`. -/
def params : Params := ⟨128, 52, 80, 68, (2 ^ 68 + 79) / 80⟩

/-- The numerators `sᵢ = 49ⁱ K!/i!`, `K = 192`, of the Poisson weights. -/
def poissonNum (i : ℕ) : ℕ := 49 ^ i * (Nat.factorial 192 / Nat.factorial i)

/-- `S_K = ∑_{i ≤ K} sᵢ`. -/
def poissonDen : ℕ := ((List.range 193).map poissonNum).sum

/-- The coefficients `tᵢ = ⌊Q sᵢ / S_K⌋`, `Q = 2⁵⁰`, `i = 0, …, 192` (eq.
(appth:coefficients)). -/
def coeffs : List ℕ := (List.range 193).map fun i ↦ 2 ^ 50 * poissonNum i / poissonDen

/-- The row weights `2 - p_w`. -/
def rowWeights (n : ℕ) : List ℕ := (List.range (3 ^ n)).map fun k ↦ 2 - par n k

/-- The column weights `1 + p_r`. -/
def colWeights (n : ℕ) : List ℕ := (List.range (3 ^ n)).map fun k ↦ 1 + par n k

/-- The packed recurrence for `(n, g)`, all `3ⁿ` columns in one run. -/
def runBoundary (n : ℕ) (g : ℤ) : Option (List ℕ) :=
  run params (band n g) ((band n g).cst params (3 ^ n)) (3 ^ n) 0 (3 ^ n) coeffs

/-- **The table entry** `k_{n,g} = ⌊10¹⁶ ∑_{w,r} x(w,r)² (2 - p_w)(1 + p_r) / (2 Q²)⌋` (eq.
(appth:quotient)), or `none` if a range test of the packed recurrence fails. -/
def quotient (n : ℕ) (g : ℤ) : Option ℕ :=
  (runBoundary n g).map fun ys ↦
    10 ^ 16 * numerator params (rowWeights n) (colWeights n) ys / (2 * (2 ^ 50) ^ 2)

/-- The printed integers `k_{n,g}` of eq. (appth:integer-table). -/
def paperK : ℕ → ℤ → ℕ
  | 2, g => if g = 1 then 10185296075298 else if g = -1 then 4160017987 else 583578389575
  | 3, g => if g = 1 then 259254948723 else if g = -1 then 571111660440989 else 47796720049308
  | 4, g => if g = 1 then 318278534893 else if g = -1 then 661347985 else 43971049724
  | 5, g => if g = 1 then 260908126900 else if g = -1 then 29917222640283 else 5659378873682
  | 6, g => if g = 1 then 16467340956 else if g = -1 then 328219854 else 4086161507
  | _, _ => 0

/-- **The hypothesis on the integer recurrence**: for every `(n, g)` of the table, the packed
recurrence succeeds and reproduces the printed integer `k_{n,g}`. Each instance is a closed
computation (`decide +kernel`). -/
def BoundaryTotals : Prop := ∀ n ∈ lengths, ∀ g ∈ fields, quotient n g = some (paperK n g)

end BoundaryThermal

end CertifiedGaps.SpinOneHeisenberg
