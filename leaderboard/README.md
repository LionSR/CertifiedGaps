# Leaderboard: certified spectral gaps with Lean proofs

This directory is the Lean library `Leaderboard` (declared in the root `lakefile.toml`). Each entry
of the leaderboard is a Lean theorem tagged `@[gap_record]` that states a uniform spectral gap.
Continuous integration builds the theorems, audits their statements and axioms, and re-checks
every declaration in the kernel.

## Record statements

`Leaderboard/Record.lean` defines two statement forms. In both, `spectrum L : (Fin L → d) → ℝ` is
the spectrum (with multiplicity) of a periodic ring of `L` unit cells, and `r`, `β` are
rationals:

```lean
-- 0 < β, 1 < r, and Δ_L > log r / β for every even L ≥ L₀
Leaderboard.Record.Unconditional spectrum L₀ r β
-- 0 < β, 1 < r, and if H holds, then Δ_L > log r / β for every even L ≥ L₀
Leaderboard.Record.Conditional H spectrum L₀ r β
```

Here `Δ_L > c` means `CertifiedGaps.GapGT (spectrum L) c`: the ground level is nondegenerate, and
every other level lies more than `c` above it. All chain certificates in this repository give a
bound of the form `log r / β₀`, and any positive gap bound can be weakened to that form. For a
model whose seeds are still proved outside Lean, the hypothesis `H` is the pair of seed
inequalities. Records cover even ring lengths only; odd rings or open chains would need new
record forms.

## Records

| file | records |
|---|---|
| `Records/Ising.lean` | **unconditional**: the Ising chain `-ZZ - (1/4)X - (1/20)Z`, `Δ_L > log 13480 / (97/4)` for even `L ≥ 6` |
| `Records/SpinOneHeisenberg.lean` | the spin-1 Heisenberg ring, `log 20 / 784` for even `L ≥ 2304` |
| `Records/StaggeredHeisenberg.lean` | `h = 1, 3/4, 1/2, 7/20, 1/4` |
| `Records/DimerizedHeisenberg.lean` | `δ = 1/2, 3/10, 1/4` |
| `Records/SpinOneDimerized.lean` | `δ = 9/10` |
| `Records/SpinOneLargeD.lean` | `D = 4`, `D = 1603/500` |
| `Records/HeisenbergLadder.lean` | `J⊥ = 4` |
| `Records/TwoLegIsing.lean` | `g = 3` |
| `Records/IsingLongitudinal.lean` | the other 29 grid points `(g, h)` |

Every record except the first is conditional. Its hypothesis is the model's seed lemma (the `sorry`
in `lean/`) at that point, with the certified defects written out. Its proof uses only the proved
lemmas and the kernel-checked chain certificate, so a proof of the seed lemma would make the
record unconditional through `Record.Conditional.unconditional`. Three certified points have no
record because `lean/` has no chain certificate for them: the ladder at `J⊥ = 7/2, 3, 5/2`, and
the spin-1 dimerized chain at `δ = 7/10`, whose seed length is odd.

## Adding a record

1. Create a file `leaderboard/Leaderboard/Records/<Model>.lean` with the namespace
   `Leaderboard.Records.<Model>`, and import it from `leaderboard/Leaderboard.lean`.
2. State the record as a `theorem` tagged `@[gap_record]`, with a docstring that names the model
   and its parameters in words:

   ```lean
   /-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 1/4`, `h = 1/20`, unconditionally. -/
   @[gap_record] theorem g_1_4_h_1_20 :
       Record.Unconditional (IsingLongitudinal.model (1 / 4) (1 / 20)).spectrum 6 13480 (97 / 4) :=
     ...
   ```

   The arguments `L₀`, `r` and `β` must be closed terms that reduce to literals; the audit
   evaluates them in the kernel. For a reflection-positive ring (`CertifiedGaps.RPRing`) with a
   chain certificate, `Leaderboard.RPRing.conditional` produces the conditional record directly
   (see `Records/StaggeredHeisenberg.lean`).
3. Check the record locally with the three CI commands below.

Rules:

* A proof may use only the axioms `propext`, `Classical.choice` and `Quot.sound`. It may not use
  `sorry`, a custom `axiom`, or `native_decide`, which trusts the compiler. For large finite
  checks, use `decide +kernel` as `lean/SeedKernel` does.
* Leaderboard modules may not define notation, unexpanders or delaborators, because reviewers
  read the printed statements.
* Do not edit `lean/`. Generic lemmas belong in `leaderboard/Leaderboard/`.

## What CI checks

```sh
lake build Leaderboard
lake exe leaderboard_audit home_page/_data/lean_records.json
lake exe leaderboard_recheck
```

* `lake build Leaderboard` elaborates every record. Run it first, because both executables read
  the compiled `.olean` files.
* `leaderboard_audit` finds every `@[gap_record]` theorem. For each one it checks that the
  statement has the form `Record.Unconditional …` or `Record.Conditional …`, reduces `L₀`, `r`
  and `β` to literals in the kernel, and collects the axioms. It fails, and writes no file, if a
  record uses an axiom outside the three above, if a theorem states a record without the tag, or
  if a Leaderboard module changes how terms print. On success it writes one JSON object per
  record, sorted by name:

  ```json
  {"name": "Leaderboard.Records.IsingLongitudinal.g_1_4_h_1_20", "kind": "unconditional",
   "spectrum": "(CertifiedGaps.IsingLongitudinal.model (1 / 4) (1 / 20)).spectrum",
   "L0": 6, "r": "13480", "beta": "97/4", "hypothesis": null,
   "axioms": ["Classical.choice", "Quot.sound", "propext"],
   "module": "Leaderboard.Records.Ising"}
  ```

* `leaderboard_recheck` re-checks in the kernel every declaration of `Leaderboard` and of every
  repository module it imports (`CertifiedGaps.*`, `SeedKernel.*`). A file can skip the kernel
  check with `set_option debug.skipKernelTC true`; the recheck catches the resulting ill-typed
  declarations. It is the check that the toolchain's `leanchecker` performs (`leanchecker` is the
  upstreamed `lean4checker`), done over a single import of Mathlib. The option `--leanchecker`
  runs the toolchain's own `leanchecker` on each module instead; this is slower.

On this machine, `lake build Leaderboard` takes about 2.5 min from a built `lean/`, the audit about
15 s, and the recheck about 3 min.

## What a reviewer must still check

The audit cannot check that a record is about the right Hamiltonian. The printed `spectrum` must
be the spectrum of the claimed model, and its definition (for example,
`CertifiedGaps.IsingLongitudinal.model`) must be the claimed Hamiltonian, with the stated cell
structure and the stated meaning of `L`. A conditional record is only as strong as its
hypothesis. The hypothesis must be the intended seed statement, and the record does not prove
it.
