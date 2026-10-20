/-
SeedKernel: kernel-checked seed arithmetic (an opt-in library: `lake build SeedKernel`).

Every computational statement is checked by `decide +kernel` (no `native_decide`, no `sorry`);
see `SeedKernel/Audit.lean`. The spin-1 witnesses are arithmetic facts about packed integer
vectors. The Ising files prove two gaps outright (no `sorry`, standard axioms only):
`SeedKernel.Ising.ising_gap` for the rings of `IsingLongitudinal.model (1/4) (1/20)` and
`SeedKernel.IsingOpen.ising_free_gap` for its chains with free ends. lean/README.md has the file
map.
-/
import SeedKernel.Packed
import SeedKernel.PackingProofs
import SeedKernel.Obligations
import SeedKernel.PackedTests
import SeedKernel.Spin1.Packed
import SeedKernel.Spin1.Optimized
import SeedKernel.Spin1.Obligations
import SeedKernel.Spin1.Counting
import SeedKernel.Spin1.Tests
import SeedKernel.Spin1.H6
import SeedKernel.Spin1.H8
import SeedKernel.Spin1.H10
import SeedKernel.Ising.Generated
import SeedKernel.Ising.GeneratedChunks
import SeedKernel.Ising.ChunkCoverage
import SeedKernel.Ising.NegativeChecks
import SeedKernel.Spectral
import SeedKernel.Comparison
import SeedKernel.LowSpectrum
import SeedKernel.ExpBounds
import SeedKernel.Seeds
import SeedKernel.Ising.Hamiltonian
import SeedKernel.Ising.Bridge
import SeedKernel.Ising.Seed
import SeedKernel.Ising.Checks
import SeedKernel.IsingOpen.Generated
import SeedKernel.IsingOpen.GeneratedChunks
import SeedKernel.IsingOpen.Hamiltonian
import SeedKernel.IsingOpen.Bridge
import SeedKernel.IsingOpen.Seed
