/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Leaderboard.Attr
import Leaderboard.RPRing
import CertifiedGaps.Models.IsingLongitudinal

/-!
# Records for the Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`

The hypothesis of each record is the seed inequality of `CertifiedGaps.IsingLongitudinal.seeds`
at that point, with the certified defects written out; the rate `r` and `β₀` are those of the
kernel-checked chain certificate `points[i].cert`. The seeds are certified outside Lean in
`models/tfim_longitudinal`.

The first grid point `g = 1/4`, `h = 1/20` has an unconditional record instead
(`Leaderboard.Records.Ising`).
-/

namespace Leaderboard.Records.IsingLongitudinal

open CertifiedGaps CertifiedGaps.IsingLongitudinal

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 1/4`, `h = 1/10`.
Gap `Δ_L > log 23810 / (55/4)` for every even `L ≥ 6`, given the seeds
`P_s(6, 55/4) ≥ 1 - 4675e-9` and `P_τ(12, 55/4) ≥ 1 - 1799e-8`. -/
@[gap_record] theorem g_1_4_h_1_10 : Record.Conditional
    (1 - (4675e-9 : ℝ) ≤ (model (1 / 4) (1 / 10)).spatialPurity 6 (55 / 4) ∧
      1 - (1799e-8 : ℝ) ≤ (model (1 / 4) (1 / 10)).thermalPurity 12 (55 / 4))
    (model (1 / 4) (1 / 10)).spectrum 6 23810 (55 / 4) :=
  RPRing.conditional _ (by norm_num) (by decide) points[1].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 1/4`, `h = 1/4`.
Gap `Δ_L > log 48540 / (13/2)` for every even `L ≥ 6`, given the seeds
`P_s(6, 13/2) ≥ 1 - 2226e-9` and `P_τ(12, 13/2) ≥ 1 - 8857e-9`. -/
@[gap_record] theorem g_1_4_h_1_4 : Record.Conditional
    (1 - (2226e-9 : ℝ) ≤ (model (1 / 4) (1 / 4)).spatialPurity 6 (13 / 2) ∧
      1 - (8857e-9 : ℝ) ≤ (model (1 / 4) (1 / 4)).thermalPurity 12 (13 / 2))
    (model (1 / 4) (1 / 4)).spectrum 6 48540 (13 / 2) :=
  RPRing.conditional _ (by norm_num) (by decide) points[2].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 1/4`, `h = 1/2`.
Gap `Δ_L > log 52170 / (11/2)` for every even `L ≥ 6`, given the seeds
`P_s(6, 11/2) ≥ 1 - 2063e-9` and `P_τ(12, 11/2) ≥ 1 - 8242e-9`. -/
@[gap_record] theorem g_1_4_h_1_2 : Record.Conditional
    (1 - (2063e-9 : ℝ) ≤ (model (1 / 4) (1 / 2)).spatialPurity 6 (11 / 2) ∧
      1 - (8242e-9 : ℝ) ≤ (model (1 / 4) (1 / 2)).thermalPurity 12 (11 / 2))
    (model (1 / 4) (1 / 2)).spectrum 6 52170 (11 / 2) :=
  RPRing.conditional _ (by norm_num) (by decide) points[3].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 1/4`, `h = 1`.
Gap `Δ_L > log 43160 / (21/4)` for every even `L ≥ 6`, given the seeds
`P_s(6, 21/4) ≥ 1 - 2494e-9` and `P_τ(12, 21/4) ≥ 1 - 9964e-9`. -/
@[gap_record] theorem g_1_4_h_1 : Record.Conditional
    (1 - (2494e-9 : ℝ) ≤ (model (1 / 4) 1).spatialPurity 6 (21 / 4) ∧
      1 - (9964e-9 : ℝ) ≤ (model (1 / 4) 1).thermalPurity 12 (21 / 4))
    (model (1 / 4) 1).spectrum 6 43160 (21 / 4) :=
  RPRing.conditional _ (by norm_num) (by decide) points[4].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 1/2`, `h = 1/20`.
Gap `Δ_L > log 3756e-1 / 13` for every even `L ≥ 6`, given the seeds
`P_s(6, 13) ≥ 1 - 1314e-6` and `P_τ(12, 13) ≥ 1 - 1872e-8`. -/
@[gap_record] theorem g_1_2_h_1_20 : Record.Conditional
    (1 - (1314e-6 : ℝ) ≤ (model (1 / 2) (1 / 20)).spatialPurity 6 13 ∧
      1 - (1872e-8 : ℝ) ≤ (model (1 / 2) (1 / 20)).thermalPurity 12 13)
    (model (1 / 2) (1 / 20)).spectrum 6 3756e-1 13 :=
  RPRing.conditional _ (by norm_num) (by decide) points[5].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 1/2`, `h = 1/10`.
Gap `Δ_L > log 2141 / 8` for every even `L ≥ 6`, given the seeds
`P_s(6, 8) ≥ 1 - 2267e-7` and `P_τ(12, 8) ≥ 1 - 1129e-8`. -/
@[gap_record] theorem g_1_2_h_1_10 : Record.Conditional
    (1 - (2267e-7 : ℝ) ≤ (model (1 / 2) (1 / 10)).spatialPurity 6 8 ∧
      1 - (1129e-8 : ℝ) ≤ (model (1 / 2) (1 / 10)).thermalPurity 12 8)
    (model (1 / 2) (1 / 10)).spectrum 6 2141 8 :=
  RPRing.conditional _ (by norm_num) (by decide) points[6].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 1/2`, `h = 1/4`.
Gap `Δ_L > log 41080 / (9/2)` for every even `L ≥ 6`, given the seeds
`P_s(6, 9/2) ≥ 1 - 7237e-9` and `P_τ(12, 9/2) ≥ 1 - 7138e-9`. -/
@[gap_record] theorem g_1_2_h_1_4 : Record.Conditional
    (1 - (7237e-9 : ℝ) ≤ (model (1 / 2) (1 / 4)).spatialPurity 6 (9 / 2) ∧
      1 - (7138e-9 : ℝ) ≤ (model (1 / 2) (1 / 4)).thermalPurity 12 (9 / 2))
    (model (1 / 2) (1 / 4)).spectrum 6 41080 (9 / 2) :=
  RPRing.conditional _ (by norm_num) (by decide) points[7].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 1/2`, `h = 1/2`.
Gap `Δ_L > log 66010 / 4` for every even `L ≥ 6`, given the seeds
`P_s(6, 4) ≥ 1 - 1625e-9` and `P_τ(12, 4) ≥ 1 - 6517e-9`. -/
@[gap_record] theorem g_1_2_h_1_2 : Record.Conditional
    (1 - (1625e-9 : ℝ) ≤ (model (1 / 2) (1 / 2)).spatialPurity 6 4 ∧
      1 - (6517e-9 : ℝ) ≤ (model (1 / 2) (1 / 2)).thermalPurity 12 4)
    (model (1 / 2) (1 / 2)).spectrum 6 66010 4 :=
  RPRing.conditional _ (by norm_num) (by decide) points[8].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 1/2`, `h = 1`.
Gap `Δ_L > log 54530 / 4` for every even `L ≥ 6`, given the seeds
`P_s(6, 4) ≥ 1 - 1974e-9` and `P_τ(12, 4) ≥ 1 - 7885e-9`. -/
@[gap_record] theorem g_1_2_h_1 : Record.Conditional
    (1 - (1974e-9 : ℝ) ≤ (model (1 / 2) 1).spatialPurity 6 4 ∧
      1 - (7885e-9 : ℝ) ≤ (model (1 / 2) 1).thermalPurity 12 4)
    (model (1 / 2) 1).spectrum 6 54530 4 :=
  RPRing.conditional _ (by norm_num) (by decide) points[9].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 3/4`, `h = 1/20`.
Gap `Δ_L > log 1258e-2 / 8` for every even `L ≥ 6`, given the seeds
`P_s(6, 8) ≥ 1 - 371e-4` and `P_τ(12, 8) ≥ 1 - 3654e-7`. -/
@[gap_record] theorem g_3_4_h_1_20 : Record.Conditional
    (1 - (371e-4 : ℝ) ≤ (model (3 / 4) (1 / 20)).spatialPurity 6 8 ∧
      1 - (3654e-7 : ℝ) ≤ (model (3 / 4) (1 / 20)).thermalPurity 12 8)
    (model (3 / 4) (1 / 20)).spectrum 6 1258e-2 8 :=
  RPRing.conditional _ (by norm_num) (by decide) points[10].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 3/4`, `h = 1/10`.
Gap `Δ_L > log 5711e-2 / (21/4)` for every even `L ≥ 6`, given the seeds
`P_s(6, 21/4) ≥ 1 - 8573e-6` and `P_τ(12, 21/4) ≥ 1 - 6931e-8`. -/
@[gap_record] theorem g_3_4_h_1_10 : Record.Conditional
    (1 - (8573e-6 : ℝ) ≤ (model (3 / 4) (1 / 10)).spatialPurity 6 (21 / 4) ∧
      1 - (6931e-8 : ℝ) ≤ (model (3 / 4) (1 / 10)).thermalPurity 12 (21 / 4))
    (model (3 / 4) (1 / 10)).spectrum 6 5711e-2 (21 / 4) :=
  RPRing.conditional _ (by norm_num) (by decide) points[11].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 3/4`, `h = 1/4`.
Gap `Δ_L > log 1559 / (7/2)` for every even `L ≥ 6`, given the seeds
`P_s(6, 7/2) ≥ 1 - 2372e-7` and `P_τ(12, 7/2) ≥ 1 - 1331e-7`. -/
@[gap_record] theorem g_3_4_h_1_4 : Record.Conditional
    (1 - (2372e-7 : ℝ) ≤ (model (3 / 4) (1 / 4)).spatialPurity 6 (7 / 2) ∧
      1 - (1331e-7 : ℝ) ≤ (model (3 / 4) (1 / 4)).thermalPurity 12 (7 / 2))
    (model (3 / 4) (1 / 4)).spectrum 6 1559 (7 / 2) :=
  RPRing.conditional _ (by norm_num) (by decide) points[12].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 3/4`, `h = 1/2`.
Gap `Δ_L > log 37670 / (7/2)` for every even `L ≥ 6`, given the seeds
`P_s(6, 7/2) ≥ 1 - 4348e-9` and `P_τ(12, 7/2) ≥ 1 - 1065e-8`. -/
@[gap_record] theorem g_3_4_h_1_2 : Record.Conditional
    (1 - (4348e-9 : ℝ) ≤ (model (3 / 4) (1 / 2)).spatialPurity 6 (7 / 2) ∧
      1 - (1065e-8 : ℝ) ≤ (model (3 / 4) (1 / 2)).thermalPurity 12 (7 / 2))
    (model (3 / 4) (1 / 2)).spectrum 6 37670 (7 / 2) :=
  RPRing.conditional _ (by norm_num) (by decide) points[13].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 3/4`, `h = 1`.
Gap `Δ_L > log 58500 / (7/2)` for every even `L ≥ 6`, given the seeds
`P_s(6, 7/2) ≥ 1 - 1856e-9` and `P_τ(12, 7/2) ≥ 1 - 7342e-9`. -/
@[gap_record] theorem g_3_4_h_1 : Record.Conditional
    (1 - (1856e-9 : ℝ) ≤ (model (3 / 4) 1).spatialPurity 6 (7 / 2) ∧
      1 - (7342e-9 : ℝ) ≤ (model (3 / 4) 1).thermalPurity 12 (7 / 2))
    (model (3 / 4) 1).spectrum 6 58500 (7 / 2) :=
  RPRing.conditional _ (by norm_num) (by decide) points[14].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 1`, `h = 1/20`.
Gap `Δ_L > log 1633e-3 / (13/2)` for every even `L ≥ 6`, given the seeds
`P_s(6, 13/2) ≥ 1 - 1841e-4` and `P_τ(12, 13/2) ≥ 1 - 7425e-6`. -/
@[gap_record] theorem g_1_h_1_20 : Record.Conditional
    (1 - (1841e-4 : ℝ) ≤ (model 1 (1 / 20)).spatialPurity 6 (13 / 2) ∧
      1 - (7425e-6 : ℝ) ≤ (model 1 (1 / 20)).thermalPurity 12 (13 / 2))
    (model 1 (1 / 20)).spectrum 6 1633e-3 (13 / 2) :=
  RPRing.conditional _ (by norm_num) (by decide) points[15].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 1`, `h = 1/10`.
Gap `Δ_L > log 8494e-3 / (17/4)` for every even `L ≥ 6`, given the seeds
`P_s(6, 17/4) ≥ 1 - 5123e-5` and `P_τ(12, 17/4) ≥ 1 - 4667e-6`. -/
@[gap_record] theorem g_1_h_1_10 : Record.Conditional
    (1 - (5123e-5 : ℝ) ≤ (model 1 (1 / 10)).spatialPurity 6 (17 / 4) ∧
      1 - (4667e-6 : ℝ) ≤ (model 1 (1 / 10)).thermalPurity 12 (17 / 4))
    (model 1 (1 / 10)).spectrum 6 8494e-3 (17 / 4) :=
  RPRing.conditional _ (by norm_num) (by decide) points[16].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 1`, `h = 1/4`.
Gap `Δ_L > log 1244e-1 / 3` for every even `L ≥ 6`, given the seeds
`P_s(6, 3) ≥ 1 - 2469e-6` and `P_τ(12, 3) ≥ 1 - 2245e-6`. -/
@[gap_record] theorem g_1_h_1_4 : Record.Conditional
    (1 - (2469e-6 : ℝ) ≤ (model 1 (1 / 4)).spatialPurity 6 3 ∧
      1 - (2245e-6 : ℝ) ≤ (model 1 (1 / 4)).thermalPurity 12 3)
    (model 1 (1 / 4)).spectrum 6 1244e-1 3 :=
  RPRing.conditional _ (by norm_num) (by decide) points[17].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 1`, `h = 1/2`.
Gap `Δ_L > log 3839 / 3` for every even `L ≥ 6`, given the seeds
`P_s(6, 3) ≥ 1 - 5079e-8` and `P_τ(12, 3) ≥ 1 - 9929e-8`. -/
@[gap_record] theorem g_1_h_1_2 : Record.Conditional
    (1 - (5079e-8 : ℝ) ≤ (model 1 (1 / 2)).spatialPurity 6 3 ∧
      1 - (9929e-8 : ℝ) ≤ (model 1 (1 / 2)).thermalPurity 12 3)
    (model 1 (1 / 2)).spectrum 6 3839 3 :=
  RPRing.conditional _ (by norm_num) (by decide) points[18].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 1`, `h = 1`.
Gap `Δ_L > log 55100 / 3` for every even `L ≥ 6`, given the seeds
`P_s(6, 3) ≥ 1 - 2207e-9` and `P_τ(12, 3) ≥ 1 - 769e-8`. -/
@[gap_record] theorem g_1_h_1 : Record.Conditional
    (1 - (2207e-9 : ℝ) ≤ (model 1 1).spatialPurity 6 3 ∧
      1 - (769e-8 : ℝ) ≤ (model 1 1).thermalPurity 12 3)
    (model 1 1).spectrum 6 55100 3 :=
  RPRing.conditional _ (by norm_num) (by decide) points[19].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 3/2`, `h = 1/20`.
Gap `Δ_L > log 1798e-3 / (7/2)` for every even `L ≥ 6`, given the seeds
`P_s(6, 7/2) ≥ 1 - 1478e-4` and `P_τ(12, 7/2) ≥ 1 - 5424e-5`. -/
@[gap_record] theorem g_3_2_h_1_20 : Record.Conditional
    (1 - (1478e-4 : ℝ) ≤ (model (3 / 2) (1 / 20)).spatialPurity 6 (7 / 2) ∧
      1 - (5424e-5 : ℝ) ≤ (model (3 / 2) (1 / 20)).thermalPurity 12 (7 / 2))
    (model (3 / 2) (1 / 20)).spectrum 6 1798e-3 (7 / 2) :=
  RPRing.conditional _ (by norm_num) (by decide) points[20].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 3/2`, `h = 1/10`.
Gap `Δ_L > log 321e-2 / (13/4)` for every even `L ≥ 6`, given the seeds
`P_s(6, 13/4) ≥ 1 - 9848e-5` and `P_τ(12, 13/4) ≥ 1 - 4085e-5`. -/
@[gap_record] theorem g_3_2_h_1_10 : Record.Conditional
    (1 - (9848e-5 : ℝ) ≤ (model (3 / 2) (1 / 10)).spatialPurity 6 (13 / 4) ∧
      1 - (4085e-5 : ℝ) ≤ (model (3 / 2) (1 / 10)).thermalPurity 12 (13 / 4))
    (model (3 / 2) (1 / 10)).spectrum 6 321e-2 (13 / 4) :=
  RPRing.conditional _ (by norm_num) (by decide) points[21].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 3/2`, `h = 1/4`.
Gap `Δ_L > log 1815e-2 / (11/4)` for every even `L ≥ 6`, given the seeds
`P_s(6, 11/4) ≥ 1 - 1818e-5` and `P_τ(12, 11/4) ≥ 1 - 1278e-5`. -/
@[gap_record] theorem g_3_2_h_1_4 : Record.Conditional
    (1 - (1818e-5 : ℝ) ≤ (model (3 / 2) (1 / 4)).spatialPurity 6 (11 / 4) ∧
      1 - (1278e-5 : ℝ) ≤ (model (3 / 2) (1 / 4)).thermalPurity 12 (11 / 4))
    (model (3 / 2) (1 / 4)).spectrum 6 1815e-2 (11 / 4) :=
  RPRing.conditional _ (by norm_num) (by decide) points[22].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 3/2`, `h = 1/2`.
Gap `Δ_L > log 2042e-1 / (5/2)` for every even `L ≥ 6`, given the seeds
`P_s(6, 5/2) ≥ 1 - 1195e-6` and `P_τ(12, 5/2) ≥ 1 - 1673e-6`. -/
@[gap_record] theorem g_3_2_h_1_2 : Record.Conditional
    (1 - (1195e-6 : ℝ) ≤ (model (3 / 2) (1 / 2)).spatialPurity 6 (5 / 2) ∧
      1 - (1673e-6 : ℝ) ≤ (model (3 / 2) (1 / 2)).thermalPurity 12 (5 / 2))
    (model (3 / 2) (1 / 2)).spectrum 6 2042e-1 (5 / 2) :=
  RPRing.conditional _ (by norm_num) (by decide) points[23].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 3/2`, `h = 1`.
Gap `Δ_L > log 10490 / (5/2)` for every even `L ≥ 6`, given the seeds
`P_s(6, 5/2) ≥ 1 - 2241e-8` and `P_τ(12, 5/2) ≥ 1 - 3345e-8`. -/
@[gap_record] theorem g_3_2_h_1 : Record.Conditional
    (1 - (2241e-8 : ℝ) ≤ (model (3 / 2) 1).spatialPurity 6 (5 / 2) ∧
      1 - (3345e-8 : ℝ) ≤ (model (3 / 2) 1).thermalPurity 12 (5 / 2))
    (model (3 / 2) 1).spectrum 6 10490 (5 / 2) :=
  RPRing.conditional _ (by norm_num) (by decide) points[24].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 2`, `h = 1/20`.
Gap `Δ_L > log 7755e-3 / (9/4)` for every even `L ≥ 6`, given the seeds
`P_s(6, 9/4) ≥ 1 - 3252e-5` and `P_τ(12, 9/4) ≥ 1 - 3673e-5`. -/
@[gap_record] theorem g_2_h_1_20 : Record.Conditional
    (1 - (3252e-5 : ℝ) ≤ (model 2 (1 / 20)).spatialPurity 6 (9 / 4) ∧
      1 - (3673e-5 : ℝ) ≤ (model 2 (1 / 20)).thermalPurity 12 (9 / 4))
    (model 2 (1 / 20)).spectrum 6 7755e-3 (9 / 4) :=
  RPRing.conditional _ (by norm_num) (by decide) points[25].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 2`, `h = 1/10`.
Gap `Δ_L > log 872e-2 / (9/4)` for every even `L ≥ 6`, given the seeds
`P_s(6, 9/4) ≥ 1 - 2941e-5` and `P_τ(12, 9/4) ≥ 1 - 3285e-5`. -/
@[gap_record] theorem g_2_h_1_10 : Record.Conditional
    (1 - (2941e-5 : ℝ) ≤ (model 2 (1 / 10)).spatialPurity 6 (9 / 4) ∧
      1 - (3285e-5 : ℝ) ≤ (model 2 (1 / 10)).thermalPurity 12 (9 / 4))
    (model 2 (1 / 10)).spectrum 6 872e-2 (9 / 4) :=
  RPRing.conditional _ (by norm_num) (by decide) points[26].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 2`, `h = 1/4`.
Gap `Δ_L > log 1745e-2 / (9/4)` for every even `L ≥ 6`, given the seeds
`P_s(6, 9/4) ≥ 1 - 1578e-5` and `P_τ(12, 9/4) ≥ 1 - 1669e-5`. -/
@[gap_record] theorem g_2_h_1_4 : Record.Conditional
    (1 - (1578e-5 : ℝ) ≤ (model 2 (1 / 4)).spatialPurity 6 (9 / 4) ∧
      1 - (1669e-5 : ℝ) ≤ (model 2 (1 / 4)).thermalPurity 12 (9 / 4))
    (model 2 (1 / 4)).spectrum 6 1745e-2 (9 / 4) :=
  RPRing.conditional _ (by norm_num) (by decide) points[27].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 2`, `h = 1/2`.
Gap `Δ_L > log 8976e-2 / (9/4)` for every even `L ≥ 6`, given the seeds
`P_s(6, 9/4) ≥ 1 - 3329e-6` and `P_τ(12, 9/4) ≥ 1 - 3196e-6`. -/
@[gap_record] theorem g_2_h_1_2 : Record.Conditional
    (1 - (3329e-6 : ℝ) ≤ (model 2 (1 / 2)).spatialPurity 6 (9 / 4) ∧
      1 - (3196e-6 : ℝ) ≤ (model 2 (1 / 2)).thermalPurity 12 (9 / 4))
    (model 2 (1 / 2)).spectrum 6 8976e-2 (9 / 4) :=
  RPRing.conditional _ (by norm_num) (by decide) points[28].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

/-- The Ising chain `-∑ ZⱼZⱼ₊₁ - g ∑ Xⱼ - h ∑ Zⱼ`, `g = 2`, `h = 1`.
Gap `Δ_L > log 1073 / 2` for every even `L ≥ 6`, given the seeds
`P_s(6, 2) ≥ 1 - 1302e-7` and `P_τ(12, 2) ≥ 1 - 3862e-7`. -/
@[gap_record] theorem g_2_h_1 : Record.Conditional
    (1 - (1302e-7 : ℝ) ≤ (model 2 1).spatialPurity 6 2 ∧
      1 - (3862e-7 : ℝ) ≤ (model 2 1).thermalPurity 12 2)
    (model 2 1).spectrum 6 1073 2 :=
  RPRing.conditional _ (by norm_num) (by decide) points[29].cert
    (points_valid _ (List.getElem_mem _)) (by rfl) (by norm_num) (by norm_num)
    (by norm_num [points, ChainCertificate.sAt]) (by norm_num [points, ChainCertificate.tAt])
    (by norm_num)

end Leaderboard.Records.IsingLongitudinal
