/-
SeedKernel: kernel-checked seed arithmetic (an opt-in library: `lake build SeedKernel`).

Every computational statement is checked by `decide +kernel` (no `native_decide`, no `sorry`);
see `SeedKernel/Audit.lean`. The spin-1 witnesses are arithmetic facts about packed integer
vectors. The Ising files prove one seed outright: `SeedKernel.Ising.ising_gap` is an
unconditional gap theorem for `IsingLongitudinal.model (1/4) (1/20)` (no `sorry`, standard
axioms only). The spin-1 trial state is checked outright (`trial_check_72`, `trial_check_120`:
the integer matrix-product state of bond dimension 26 on 72 and 120 sites), so the spin-1 gap
follows from the trace table alone (`SpinOneHeisenberg.gap_of_table`). lean/README.md has the file
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
import SeedKernel.Heisenberg.Trial
import SeedKernel.Heisenberg.Seed36
import SeedKernel.Heisenberg.Seed60
import SeedKernel.Heisenberg.Reseed120
import SeedKernel.Heisenberg.OpenTrial
import SeedKernel.Heisenberg.BoundaryTraces
import SeedKernel.Heisenberg.ThermalTraces
import SeedKernel.Heisenberg.Integers
import SeedKernel.Ising.Generated
import SeedKernel.Ising.GeneratedChunks
import SeedKernel.Ising.ChunkCoverage
import SeedKernel.Ising.NegativeChecks
import SeedKernel.Ising.Hamiltonian
import SeedKernel.Ising.Bridge
import SeedKernel.Ising.Seed
