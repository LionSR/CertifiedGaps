# Leaderboard: certified spectral gaps with Lean proofs

This directory holds the Lean library `Leaderboard`, declared in the root `lakefile.toml`. Each
leaderboard entry with a Lean proof is a theorem in a module `Leaderboard.Records.<Model>`. CI
builds these theorems, audits them, and re-checks them in the Lean kernel. The audit output
feeds the leaderboard page.

## Records

A record is a theorem whose statement has one of four forms, defined in
`Leaderboard/Record.lean`:

```lean
Leaderboard.Record.Unconditional spectrum L₀ r β
Leaderboard.Record.Seeded spectrum L₀ r β n₀ s t
Leaderboard.Record.SeededOn family a b L₀ r β n₀ s t
Leaderboard.Record.TwoSeeded spectrum L₀ r β n₀ s t n₁ β₁ s₁ t₁
```

`spectrum L : (Fin L → d) → ℝ` lists the eigenvalues of a periodic ring of `L` unit cells, with
multiplicity. The numbers `L₀`, `n₀` are natural numbers, and `r`, `β`, `s`, `t` are rationals.

* `Unconditional` states `0 < β`, `1 < r`, and `Δ_L > log r / β` for every even `L ≥ L₀`.
* `Seeded` states `0 < β`, `1 < r`, and the same gap, under two seed inequalities on the same
  spectrum: `P_s(n₀, β) ≥ 1 - s` and `P_τ(2n₀, β) ≥ 1 - t`.
* `SeededOn` states `Seeded (family x) L₀ r β n₀ s t` for every real parameter `x` with
  `a ≤ x ≤ b`. It is the form for a computer certificate that bounds the two defects on a whole
  parameter interval, for example by a perturbation bound.
* `TwoSeeded` states the same gap under the seed inequalities at `(n₀, β)` and at a second seed
  `(n₁, β₁)` with defects `s₁`, `t₁`. The second seed covers the rings from `n₁` to `n₀` in one
  doubling step, as for the temperature-lifted spin-1 seed (`Leaderboard/TwoSeeded.lean`).

`Δ_L > c` means `CertifiedGaps.GapGT (spectrum L) c`. The ground level is nondegenerate, and
every other level lies more than `c` above it. The purities are `CertifiedGaps.spatialPurity`
and `CertifiedGaps.thermalPurity` of the spectrum.

The definition fixes the hypothesis of a seeded record. A record cannot assume anything other
than the seed inequalities, and the numbers `n₀`, `s`, `t` (and `n₁`, `β₁`, `s₁`, `t₁`) determine
them. The seeds are
certified by computer in `models/<slug>/`, and a validator compares them with the certificate.
Once the seeds are proved in Lean, `Record.Seeded.unconditional` makes the record unconditional.

All chain certificates in this repository give bounds of the form `log r / β`. Records cover
even ring lengths only. Odd rings and open chains would need new record forms.

## Adding a record

1. Put the theorem in `leaderboard/Leaderboard/Records/<Model>.lean`, in the namespace
   `Leaderboard.Records.<Model>`. Import the file from `leaderboard/Leaderboard.lean`.
2. Give the theorem a one-line docstring that names the model and its parameters.
3. Add the row to `home_page/_data/records.yml` with a `lean:` entry naming the theorem, and run

   ```sh
   julia --project=tools/leaderboard tools/leaderboard/lean_chains.jl
   ```

   It reads the exact seed defects from the row's certificate, computes a chain certificate
   (14 doublings, 8 significant digits) and writes it to
   `leaderboard/Leaderboard/Chains/<Model>.lean` as `Leaderboard.Chains.<Model>.<point>`, with a
   theorem `<point>_valid` that the kernel checks by `decide +kernel`. It also fills in the
   row's `lean.r` and `lean.decimal`. The script only searches: Lean checks the chain, and the
   validator checks the seeds against the certificate.
4. Write the parameters and the numbers as literals. For a model in `lean/CertifiedGaps/Models`,
   `Leaderboard.RPRing.seeded` proves the record in one line from the chain:

   ```lean
   /-- The two-leg transverse-field Ising ladder, `g = 3`. -/
   theorem g_3 :
       Record.Seeded (model 3).spectrum 4 25974347e-7 (3 / 2) 4 11240837e-8 51615239e-9 :=
     RPRing.seeded (model 3) _ Chains.TwoLegIsing.g_3_valid
   ```

   The literals are the `r` and the first entries of `s` and `t` of the chain, copied as they
   are (`rfl` checks them); the chain itself is inferred from its validity theorem. For another
   spectrum with the moment inequalities, use `Record.Seeded.of_chain`, and for two seeds
   `Record.TwoSeeded.of_chain`, as in `Records/SpinOneHeisenberg.lean`.

   For a physical frame, add `|>.congr` with the spectrum equality, as in
   `Records/StaggeredHeisenberg.lean`.
5. Run the CI commands below.

Proofs may use only the axioms `propext`, `Classical.choice` and `Quot.sound`. Do not use
`sorry`, a new `axiom`, or `native_decide`, which trusts the compiler. For large finite checks,
use `decide +kernel`, as `lean/SeedKernel` does. Do not edit `lean/`; put new generic lemmas in
`leaderboard/Leaderboard/`.

## What CI checks

```sh
lake build
lake exe leaderboard_audit home_page/_data/lean_records.json
LEAN_NUM_THREADS=1 lake env leanchecker Leaderboard CertifiedGaps SeedKernel.Ising SeedKernel.Packed
```

`lake build` builds every default target, `Leaderboard` included. Run it first, because the
audit reads the compiled `.olean` files.

The audit lists every theorem of a `Leaderboard.Records.*` module whose statement is a record.
For each record it does three things.

* It collects the axioms of the proof, and fails if any is outside the three above.
* It reduces `L₀`, `r`, `β`, `n₀`, `s`, `t` to exact literals in the kernel.
* It reads the model parameters from the spectrum term.

If anything fails, the audit writes no file and exits 1. Otherwise it writes a JSON array,
sorted by name, with one object per record:

```json
{"name": "Leaderboard.Records.TwoLegIsing.g_3",
 "module": "Leaderboard.Records.TwoLegIsing",
 "kind": "conditional",
 "spectrum": {"model": "CertifiedGaps.TwoLegIsing.model", "frame": "spectrum",
              "params": ["3"], "printed": "CertifiedGaps.RPRing.spectrum (CertifiedGaps.TwoLegIsing.model 3)"},
 "L0": 4, "r": "25974347/10000000", "beta": "3/2", "n0": 4, "s": "11240837/100000000",
 "t": "51615239/1000000000", "interval": null, "seed2": null,
 "axioms": ["Classical.choice", "Quot.sound", "propext"]}
```

* `kind` is `unconditional` or `conditional` (any of the seeded forms). For an unconditional
  record, `n0`, `s`, `t` are `null`.
* `interval` is `[a, b]` for a `SeededOn` record, and `seed2` is `{n0, beta, s, t}` of the
  second seed of a `TwoSeeded` record; otherwise both are `null`.
* In `spectrum`, `model` is the constant that takes the parameters, and `params` are its
  explicit arguments as exact rationals. `frame` is `spectrum` for `R.spectrum` of a
  reflection-positive ring `R`; otherwise it is the last part of the name, such as
  `physicalSpectrum`. For the family `fun x ↦ …` of a `SeededOn` record, the parameter is
  listed as `"x"`. `printed` is the full term with full names.

`leanchecker` comes with the Lean toolchain; it is the upstreamed `lean4checker`. It replays
every declaration of the listed modules through the kernel. This catches a file that switched
the kernel off with `set_option debug.skipKernelTC true`. The listed prefixes cover the
Leaderboard library and the repository modules it imports. With one thread, leanchecker checks
one module at a time, so its memory peak is that of the largest module, about 9 GB for
`SeedKernel.Ising.Seed`. Parallel checks need much more memory. CI runs this check in
`.github/workflows/kernel-check.yml`, apart from the deploy. It runs when `leaderboard/` or
`lean/` changes, and weekly.

## What a reviewer must still check

The audit cannot check that a record is about the right Hamiltonian. A reviewer must check that
the model and frame are the claimed Hamiltonian. That includes the unit cell and what `L`
counts.
