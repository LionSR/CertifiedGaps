/-
Copyright (c) 2026 CertifiedGaps contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: CertifiedGaps contributors
-/
import Leaderboard.RPRing
import Leaderboard.Chains.IsingLongitudinal
import SeedKernel.Ising.Seed

/-! # Records: the longitudinal-field Ising chain -/

namespace Leaderboard.Records.IsingLongitudinal

open CertifiedGaps CertifiedGaps.IsingLongitudinal

/-- The longitudinal-field Ising chain, `g = 1/4`, `h = 1/20`, with its seed proved in Lean. -/
theorem g_1_4_h_1_20 :
    Record.Unconditional (model (1 / 4) (1 / 20)).spectrum 6 13531935e-3 (97 / 4) :=
  (RPRing.seeded (model (1 / 4) (1 / 20)) _ Chains.IsingLongitudinal.g_1_4_h_1_20_valid
    (n₀ := 6) (s := 9064e-9) (t := 3141e-8)).unconditional
    (RPRing.seedBounds _ (by push_cast; exact SeedKernel.Ising.seed_spatial)
      (by push_cast; exact SeedKernel.Ising.seed_thermal))

/-- The longitudinal-field Ising chain, `g = 1/4`, `h = 1/10`. -/
theorem g_1_4_h_1_10 :
    Record.Seeded (model (1 / 4) (1 / 10)).spectrum 6 23909642e-3 (55 / 4) 6
      46744447e-13 17985668e-12 :=
  RPRing.seeded (model (1 / 4) (1 / 10)) _ Chains.IsingLongitudinal.g_1_4_h_1_10_valid

/-- The longitudinal-field Ising chain, `g = 1/4`, `h = 1/4`. -/
theorem g_1_4_h_1_4 :
    Record.Seeded (model (1 / 4) (1 / 4)).spectrum 6 48720928e-3 (13 / 2) 6
      22252215e-13 88566881e-13 :=
  RPRing.seeded (model (1 / 4) (1 / 4)) _ Chains.IsingLongitudinal.g_1_4_h_1_4_valid

/-- The longitudinal-field Ising chain, `g = 1/4`, `h = 1/2`. -/
theorem g_1_4_h_1_2 :
    Record.Seeded (model (1 / 4) (1 / 2)).spectrum 6 52375218e-3 (11 / 2) 6
      20629271e-13 8241803e-12 :=
  RPRing.seeded (model (1 / 4) (1 / 2)) _ Chains.IsingLongitudinal.g_1_4_h_1_2_valid

/-- The longitudinal-field Ising chain, `g = 1/4`, `h = 1`. -/
theorem g_1_4_h_1 :
    Record.Seeded (model (1 / 4) 1).spectrum 6 43323337e-3 (21 / 4) 6 2493943e-12 99637962e-13 :=
  RPRing.seeded (model (1 / 4) 1) _ Chains.IsingLongitudinal.g_1_4_h_1_valid

/-- The longitudinal-field Ising chain, `g = 1/2`, `h = 1/20`. -/
theorem g_1_2_h_1_20 :
    Record.Seeded (model (1 / 2) (1 / 20)).spectrum 6 3771579e-4 13 6 13135704e-10 18719939e-12 :=
  RPRing.seeded (model (1 / 2) (1 / 20)) _ Chains.IsingLongitudinal.g_1_2_h_1_20_valid

/-- The longitudinal-field Ising chain, `g = 1/2`, `h = 1/10`. -/
theorem g_1_2_h_1_10 :
    Record.Seeded (model (1 / 2) (1 / 10)).spectrum 6 21500127e-4 8 6 22668546e-11 1128044e-11 :=
  RPRing.seeded (model (1 / 2) (1 / 10)) _ Chains.IsingLongitudinal.g_1_2_h_1_10_valid

/-- The longitudinal-field Ising chain, `g = 1/2`, `h = 1/4`. -/
theorem g_1_2_h_1_4 :
    Record.Seeded (model (1 / 2) (1 / 4)).spectrum 6 41239512e-3 (9 / 2) 6
      72364808e-13 71371608e-13 :=
  RPRing.seeded (model (1 / 2) (1 / 4)) _ Chains.IsingLongitudinal.g_1_2_h_1_4_valid

/-- The longitudinal-field Ising chain, `g = 1/2`, `h = 1/2`. -/
theorem g_1_2_h_1_2 :
    Record.Seeded (model (1 / 2) (1 / 2)).spectrum 6 66266676e-3 4 6 16249764e-13 65164761e-13 :=
  RPRing.seeded (model (1 / 2) (1 / 2)) _ Chains.IsingLongitudinal.g_1_2_h_1_2_valid

/-- The longitudinal-field Ising chain, `g = 1/2`, `h = 1`. -/
theorem g_1_2_h_1 :
    Record.Seeded (model (1 / 2) 1).spectrum 6 54746212e-3 4 6 19735499e-13 78848807e-13 :=
  RPRing.seeded (model (1 / 2) 1) _ Chains.IsingLongitudinal.g_1_2_h_1_valid

/-- The longitudinal-field Ising chain, `g = 3/4`, `h = 1/20`. -/
theorem g_3_4_h_1_20 :
    Record.Seeded (model (3 / 4) (1 / 20)).spectrum 6 12631258e-6 8 6 37097109e-9 36532524e-11 :=
  RPRing.seeded (model (3 / 4) (1 / 20)) _ Chains.IsingLongitudinal.g_3_4_h_1_20_valid

/-- The longitudinal-field Ising chain, `g = 3/4`, `h = 1/10`. -/
theorem g_3_4_h_1_10 :
    Record.Seeded (model (3 / 4) (1 / 10)).spectrum 6 57330474e-6 (21 / 4) 6
      85725139e-10 69305257e-12 :=
  RPRing.seeded (model (3 / 4) (1 / 10)) _ Chains.IsingLongitudinal.g_3_4_h_1_10_valid

/-- The longitudinal-field Ising chain, `g = 3/4`, `h = 1/4`. -/
theorem g_3_4_h_1_4 :
    Record.Seeded (model (3 / 4) (1 / 4)).spectrum 6 1566019e-3 (7 / 2) 6
      23714941e-11 13301067e-11 :=
  RPRing.seeded (model (3 / 4) (1 / 4)) _ Chains.IsingLongitudinal.g_3_4_h_1_4_valid

/-- The longitudinal-field Ising chain, `g = 3/4`, `h = 1/2`. -/
theorem g_3_4_h_1_2 :
    Record.Seeded (model (3 / 4) (1 / 2)).spectrum 6 37835372e-3 (7 / 2) 6
      43472772e-13 1064594e-11 :=
  RPRing.seeded (model (3 / 4) (1 / 2)) _ Chains.IsingLongitudinal.g_3_4_h_1_2_valid

/-- The longitudinal-field Ising chain, `g = 3/4`, `h = 1`. -/
theorem g_3_4_h_1 :
    Record.Seeded (model (3 / 4) 1).spectrum 6 58743294e-3 (7 / 2) 6 18555903e-13 7341292e-12 :=
  RPRing.seeded (model (3 / 4) 1) _ Chains.IsingLongitudinal.g_3_4_h_1_valid

/-- The longitudinal-field Ising chain, `g = 1`, `h = 1/20`. -/
theorem g_1_h_1_20 :
    Record.Seeded (model 1 (1 / 20)).spectrum 6 16422325e-7 (13 / 2) 6 18400219e-8 74240448e-10 :=
  RPRing.seeded (model 1 (1 / 20)) _ Chains.IsingLongitudinal.g_1_h_1_20_valid

/-- The longitudinal-field Ising chain, `g = 1`, `h = 1/10`. -/
theorem g_1_h_1_10 :
    Record.Seeded (model 1 (1 / 10)).spectrum 6 85257569e-7 (17 / 4) 6 51226486e-9 4666523e-9 :=
  RPRing.seeded (model 1 (1 / 10)) _ Chains.IsingLongitudinal.g_1_h_1_10_valid

/-- The longitudinal-field Ising chain, `g = 1`, `h = 1/4`. -/
theorem g_1_h_1_4 :
    Record.Seeded (model 1 (1 / 4)).spectrum 6 12491264e-5 3 6 24686783e-10 22440912e-10 :=
  RPRing.seeded (model 1 (1 / 4)) _ Chains.IsingLongitudinal.g_1_h_1_4_valid

/-- The longitudinal-field Ising chain, `g = 1`, `h = 1/2`. -/
theorem g_1_h_1_2 :
    Record.Seeded (model 1 (1 / 2)).spectrum 6 38538963e-4 3 6 50787566e-12 9928824e-11 :=
  RPRing.seeded (model 1 (1 / 2)) _ Chains.IsingLongitudinal.g_1_h_1_2_valid

/-- The longitudinal-field Ising chain, `g = 1`, `h = 1`. -/
theorem g_1_h_1 :
    Record.Seeded (model 1 1).spectrum 6 55318466e-3 3 6 2206019e-12 76891336e-13 :=
  RPRing.seeded (model 1 1) _ Chains.IsingLongitudinal.g_1_h_1_valid

/-- The longitudinal-field Ising chain, `g = 3/2`, `h = 1/20`. -/
theorem g_3_2_h_1_20 :
    Record.Seeded (model (3 / 2) (1 / 20)).spectrum 6 18076225e-7 (7 / 2) 6
      14773446e-8 54233031e-9 :=
  RPRing.seeded (model (3 / 2) (1 / 20)) _ Chains.IsingLongitudinal.g_3_2_h_1_20_valid

/-- The longitudinal-field Ising chain, `g = 3/2`, `h = 1/10`. -/
theorem g_3_2_h_1_10 :
    Record.Seeded (model (3 / 2) (1 / 10)).spectrum 6 32227995e-7 (13 / 4) 6
      98470933e-9 40840418e-9 :=
  RPRing.seeded (model (3 / 2) (1 / 10)) _ Chains.IsingLongitudinal.g_3_2_h_1_10_valid

/-- The longitudinal-field Ising chain, `g = 3/2`, `h = 1/4`. -/
theorem g_3_2_h_1_4 :
    Record.Seeded (model (3 / 2) (1 / 4)).spectrum 6 18227075e-6 (11 / 4) 6
      18177525e-9 12773079e-9 :=
  RPRing.seeded (model (3 / 2) (1 / 4)) _ Chains.IsingLongitudinal.g_3_2_h_1_4_valid

/-- The longitudinal-field Ising chain, `g = 3/2`, `h = 1/2`. -/
theorem g_3_2_h_1_2 :
    Record.Seeded (model (3 / 2) (1 / 2)).spectrum 6 20502425e-5 (5 / 2) 6
      11947573e-10 16723226e-10 :=
  RPRing.seeded (model (3 / 2) (1 / 2)) _ Chains.IsingLongitudinal.g_3_2_h_1_2_valid

/-- The longitudinal-field Ising chain, `g = 3/2`, `h = 1`. -/
theorem g_3_2_h_1 :
    Record.Seeded (model (3 / 2) 1).spectrum 6 10534123e-3 (5 / 2) 6 22402302e-12 33442686e-12 :=
  RPRing.seeded (model (3 / 2) 1) _ Chains.IsingLongitudinal.g_3_2_h_1_valid

/-- The longitudinal-field Ising chain, `g = 2`, `h = 1/20`. -/
theorem g_2_h_1_20 :
    Record.Seeded (model 2 (1 / 20)).spectrum 6 77844656e-7 (9 / 4) 6 32513643e-9 3672724e-8 :=
  RPRing.seeded (model 2 (1 / 20)) _ Chains.IsingLongitudinal.g_2_h_1_20_valid

/-- The longitudinal-field Ising chain, `g = 2`, `h = 1/10`. -/
theorem g_2_h_1_10 :
    Record.Seeded (model 2 (1 / 10)).spectrum 6 87540205e-7 (9 / 4) 6 29403014e-9 3284532e-8 :=
  RPRing.seeded (model 2 (1 / 10)) _ Chains.IsingLongitudinal.g_2_h_1_10_valid

/-- The longitudinal-field Ising chain, `g = 2`, `h = 1/4`. -/
theorem g_2_h_1_4 :
    Record.Seeded (model 2 (1 / 4)).spectrum 6 17526371e-6 (9 / 4) 6 15774851e-9 16683003e-9 :=
  RPRing.seeded (model 2 (1 / 4)) _ Chains.IsingLongitudinal.g_2_h_1_4_valid

/-- The longitudinal-field Ising chain, `g = 2`, `h = 1/2`. -/
theorem g_2_h_1_2 :
    Record.Seeded (model 2 (1 / 2)).spectrum 6 90114498e-6 (9 / 4) 6 33287485e-10 31950636e-10 :=
  RPRing.seeded (model 2 (1 / 2)) _ Chains.IsingLongitudinal.g_2_h_1_2_valid

/-- The longitudinal-field Ising chain, `g = 2`, `h = 1`. -/
theorem g_2_h_1 :
    Record.Seeded (model 2 1).spectrum 6 10773755e-4 2 6 13015338e-11 38610967e-11 :=
  RPRing.seeded (model 2 1) _ Chains.IsingLongitudinal.g_2_h_1_valid

end Leaderboard.Records.IsingLongitudinal
