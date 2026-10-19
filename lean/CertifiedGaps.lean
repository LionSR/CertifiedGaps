-- The abstract criterion and its rational certificates
import CertifiedGaps.PurityCriterion
import CertifiedGaps.Chain
import CertifiedGaps.ChainCheck
-- Transfer-matrix tools
import CertifiedGaps.Transfer.PiTensor
import CertifiedGaps.Transfer.MomentForm
import CertifiedGaps.Transfer.ProductLimit
-- Continuous-time history prerequisites
import CertifiedGaps.Transfer.TimeOrdering
import CertifiedGaps.Transfer.PoissonHistory
import CertifiedGaps.Transfer.HistoryShuffle
import CertifiedGaps.Transfer.ProductPoissonHistory
-- Weighted histories for arbitrary atomless time laws
import CertifiedGaps.Transfer.WeightedProductPoissonHistory
-- Reflection-positive rings
import CertifiedGaps.Ring.Model
import CertifiedGaps.Ring.Transfer
import CertifiedGaps.Ring.Gap
import CertifiedGaps.Ring.Conjugation
-- The spin-one Heisenberg ring (Theorem 1 of the note)
import CertifiedGaps.Heisenberg.Hamiltonian
import CertifiedGaps.Heisenberg.Gap
import CertifiedGaps.Heisenberg.Main
-- Open chains with boundary fields (Appendices A–C of the note)
import CertifiedGaps.Open.MomentBounds
import CertifiedGaps.Open.Criterion
import CertifiedGaps.Open.Certificate
import CertifiedGaps.Open.Model
import CertifiedGaps.Open.Transfer
import CertifiedGaps.Open.Bounds
import CertifiedGaps.Open.Gap
import CertifiedGaps.Heisenberg.OpenChain
-- Local building blocks of the models
import CertifiedGaps.Local.Pauli
import CertifiedGaps.Local.SpinHalf
import CertifiedGaps.Local.SpinHalfRotation
import CertifiedGaps.Local.SpinOne
import CertifiedGaps.Local.CellSwap
-- Models (one file per directory `models/<slug>`; `*Physical` for rotated frames)
import CertifiedGaps.Models.SpinOneLargeD
import CertifiedGaps.Models.SpinOneDimerized
import CertifiedGaps.Models.HeisenbergLadder
import CertifiedGaps.Models.HeisenbergLadderPhysical
import CertifiedGaps.Models.StaggeredHeisenberg
import CertifiedGaps.Models.StaggeredHeisenbergPhysical
import CertifiedGaps.Models.DimerizedHeisenberg
import CertifiedGaps.Models.DimerizedHeisenbergPhysical
import CertifiedGaps.Models.TwoLegIsing
import CertifiedGaps.Models.IsingLongitudinal
