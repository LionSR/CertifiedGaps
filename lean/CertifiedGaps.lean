-- The abstract criterion and its rational certificates
import CertifiedGaps.PurityCriterion
import CertifiedGaps.Chain
import CertifiedGaps.ChainCheck
-- Transfer-matrix tools
import CertifiedGaps.Transfer.PiTensor
import CertifiedGaps.Transfer.MomentForm
import CertifiedGaps.Transfer.ProductLimit
import CertifiedGaps.Transfer.PolynomialFilter
import CertifiedGaps.Transfer.MatrixProduct
import CertifiedGaps.Transfer.MatrixLemmas
import CertifiedGaps.Transfer.SingleSite
import CertifiedGaps.Transfer.SingleSiteLimit
import CertifiedGaps.Transfer.SingleSiteSectors
import CertifiedGaps.Transfer.SingleSiteD2
import CertifiedGaps.Transfer.CyclicSectors
import CertifiedGaps.Transfer.SingleSiteA4
import CertifiedGaps.Transfer.SpectralForm
-- Spectral counting: eigenvalue counts, inertia, Schur bounds, trace bounds, operator bounds
import CertifiedGaps.Spectral.Basic
import CertifiedGaps.Spectral.OperatorBounds
import CertifiedGaps.Spectral.FunctionalCalculus
import CertifiedGaps.Spectral.Counting
import CertifiedGaps.Spectral.Inertia
import CertifiedGaps.Spectral.Schur
import CertifiedGaps.Spectral.TraceBounds
-- Seeds from short-ring trace tables and trial states
import CertifiedGaps.Seed.CappedMass
import CertifiedGaps.Seed.Scheme
import CertifiedGaps.Seed.Limits
import CertifiedGaps.Seed.Dominant
import CertifiedGaps.Seed.Certificate
import CertifiedGaps.Seed.Ring
import CertifiedGaps.Seed.Lift
import CertifiedGaps.Seed.LiftedSeed
-- Reflection-positive rings
import CertifiedGaps.Ring.Hamiltonian
import CertifiedGaps.Ring.Model
import CertifiedGaps.Ring.Transfer
import CertifiedGaps.Ring.Gap
import CertifiedGaps.Ring.Conjugation
import CertifiedGaps.Ring.Variational
import CertifiedGaps.Ring.EnergyBound
import CertifiedGaps.Ring.Symmetry
-- Rigorous Gibbs-trace enclosures: Poisson truncation, integer Horner columns, Frobenius form,
-- representative columns
import CertifiedGaps.Trace.Poisson
import CertifiedGaps.Trace.Horner
import CertifiedGaps.Trace.Frobenius
import CertifiedGaps.Trace.Orbit
-- The spin-one Heisenberg ring (Theorem 1 of the note)
import CertifiedGaps.Heisenberg.Hamiltonian
import CertifiedGaps.Heisenberg.Gap
import CertifiedGaps.Heisenberg.Main
import CertifiedGaps.Heisenberg.SzBridge
import CertifiedGaps.Heisenberg.TraceTable
import CertifiedGaps.Heisenberg.SectorTransfer
import CertifiedGaps.Heisenberg.SectorTransferC
import CertifiedGaps.Heisenberg.TrialBound
import CertifiedGaps.Heisenberg.Seed36
import CertifiedGaps.Heisenberg.Seed60
import CertifiedGaps.Heisenberg.EnergyBound
import CertifiedGaps.Heisenberg.ThermalTable
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
import CertifiedGaps.Local.SpinOneRotation
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
