/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import Leaderboard.RPRing
import SeedKernel.Ising.Seed

/-! # Records: the longitudinal-field Ising chain -/

namespace Leaderboard.Records.IsingLongitudinal

open CertifiedGaps CertifiedGaps.IsingLongitudinal

/-- The longitudinal-field Ising chain, `g = 1/4`, `h = 1/20`, with its seed proved in Lean. -/
theorem g_1_4_h_1_20 :
    Record.Unconditional (model (1 / 4) (1 / 20)).spectrum 6 13480 (97 / 4) :=
  ⟨⟨by norm_num, by norm_num⟩, fun _ hL hLe ↦ by
    push_cast
    exact SeedKernel.Ising.ising_gap hL hLe⟩

/-- The longitudinal-field Ising chain, `g = 1/4`, `h = 1/10`. -/
theorem g_1_4_h_1_10 :
    Record.Seeded (model (1 / 4) (1 / 10)).spectrum 6 23810 (55 / 4) 6 4675e-9 1799e-8 :=
  RPRing.seeded (model (1 / 4) (1 / 10)) points[1].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 1/4`, `h = 1/4`. -/
theorem g_1_4_h_1_4 :
    Record.Seeded (model (1 / 4) (1 / 4)).spectrum 6 48540 (13 / 2) 6 2226e-9 8857e-9 :=
  RPRing.seeded (model (1 / 4) (1 / 4)) points[2].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 1/4`, `h = 1/2`. -/
theorem g_1_4_h_1_2 :
    Record.Seeded (model (1 / 4) (1 / 2)).spectrum 6 52170 (11 / 2) 6 2063e-9 8242e-9 :=
  RPRing.seeded (model (1 / 4) (1 / 2)) points[3].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 1/4`, `h = 1`. -/
theorem g_1_4_h_1 :
    Record.Seeded (model (1 / 4) 1).spectrum 6 43160 (21 / 4) 6 2494e-9 9964e-9 :=
  RPRing.seeded (model (1 / 4) 1) points[4].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 1/2`, `h = 1/20`. -/
theorem g_1_2_h_1_20 :
    Record.Seeded (model (1 / 2) (1 / 20)).spectrum 6 3756e-1 13 6 1314e-6 1872e-8 :=
  RPRing.seeded (model (1 / 2) (1 / 20)) points[5].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 1/2`, `h = 1/10`. -/
theorem g_1_2_h_1_10 :
    Record.Seeded (model (1 / 2) (1 / 10)).spectrum 6 2141 8 6 2267e-7 1129e-8 :=
  RPRing.seeded (model (1 / 2) (1 / 10)) points[6].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 1/2`, `h = 1/4`. -/
theorem g_1_2_h_1_4 :
    Record.Seeded (model (1 / 2) (1 / 4)).spectrum 6 41080 (9 / 2) 6 7237e-9 7138e-9 :=
  RPRing.seeded (model (1 / 2) (1 / 4)) points[7].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 1/2`, `h = 1/2`. -/
theorem g_1_2_h_1_2 :
    Record.Seeded (model (1 / 2) (1 / 2)).spectrum 6 66010 4 6 1625e-9 6517e-9 :=
  RPRing.seeded (model (1 / 2) (1 / 2)) points[8].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 1/2`, `h = 1`. -/
theorem g_1_2_h_1 :
    Record.Seeded (model (1 / 2) 1).spectrum 6 54530 4 6 1974e-9 7885e-9 :=
  RPRing.seeded (model (1 / 2) 1) points[9].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 3/4`, `h = 1/20`. -/
theorem g_3_4_h_1_20 :
    Record.Seeded (model (3 / 4) (1 / 20)).spectrum 6 1258e-2 8 6 371e-4 3654e-7 :=
  RPRing.seeded (model (3 / 4) (1 / 20)) points[10].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 3/4`, `h = 1/10`. -/
theorem g_3_4_h_1_10 :
    Record.Seeded (model (3 / 4) (1 / 10)).spectrum 6 5711e-2 (21 / 4) 6 8573e-6 6931e-8 :=
  RPRing.seeded (model (3 / 4) (1 / 10)) points[11].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 3/4`, `h = 1/4`. -/
theorem g_3_4_h_1_4 :
    Record.Seeded (model (3 / 4) (1 / 4)).spectrum 6 1559 (7 / 2) 6 2372e-7 1331e-7 :=
  RPRing.seeded (model (3 / 4) (1 / 4)) points[12].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 3/4`, `h = 1/2`. -/
theorem g_3_4_h_1_2 :
    Record.Seeded (model (3 / 4) (1 / 2)).spectrum 6 37670 (7 / 2) 6 4348e-9 1065e-8 :=
  RPRing.seeded (model (3 / 4) (1 / 2)) points[13].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 3/4`, `h = 1`. -/
theorem g_3_4_h_1 :
    Record.Seeded (model (3 / 4) 1).spectrum 6 58500 (7 / 2) 6 1856e-9 7342e-9 :=
  RPRing.seeded (model (3 / 4) 1) points[14].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 1`, `h = 1/20`. -/
theorem g_1_h_1_20 :
    Record.Seeded (model 1 (1 / 20)).spectrum 6 1633e-3 (13 / 2) 6 1841e-4 7425e-6 :=
  RPRing.seeded (model 1 (1 / 20)) points[15].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 1`, `h = 1/10`. -/
theorem g_1_h_1_10 :
    Record.Seeded (model 1 (1 / 10)).spectrum 6 8494e-3 (17 / 4) 6 5123e-5 4667e-6 :=
  RPRing.seeded (model 1 (1 / 10)) points[16].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 1`, `h = 1/4`. -/
theorem g_1_h_1_4 :
    Record.Seeded (model 1 (1 / 4)).spectrum 6 1244e-1 3 6 2469e-6 2245e-6 :=
  RPRing.seeded (model 1 (1 / 4)) points[17].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 1`, `h = 1/2`. -/
theorem g_1_h_1_2 :
    Record.Seeded (model 1 (1 / 2)).spectrum 6 3839 3 6 5079e-8 9929e-8 :=
  RPRing.seeded (model 1 (1 / 2)) points[18].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 1`, `h = 1`. -/
theorem g_1_h_1 :
    Record.Seeded (model 1 1).spectrum 6 55100 3 6 2207e-9 769e-8 :=
  RPRing.seeded (model 1 1) points[19].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 3/2`, `h = 1/20`. -/
theorem g_3_2_h_1_20 :
    Record.Seeded (model (3 / 2) (1 / 20)).spectrum 6 1798e-3 (7 / 2) 6 1478e-4 5424e-5 :=
  RPRing.seeded (model (3 / 2) (1 / 20)) points[20].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 3/2`, `h = 1/10`. -/
theorem g_3_2_h_1_10 :
    Record.Seeded (model (3 / 2) (1 / 10)).spectrum 6 321e-2 (13 / 4) 6 9848e-5 4085e-5 :=
  RPRing.seeded (model (3 / 2) (1 / 10)) points[21].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 3/2`, `h = 1/4`. -/
theorem g_3_2_h_1_4 :
    Record.Seeded (model (3 / 2) (1 / 4)).spectrum 6 1815e-2 (11 / 4) 6 1818e-5 1278e-5 :=
  RPRing.seeded (model (3 / 2) (1 / 4)) points[22].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 3/2`, `h = 1/2`. -/
theorem g_3_2_h_1_2 :
    Record.Seeded (model (3 / 2) (1 / 2)).spectrum 6 2042e-1 (5 / 2) 6 1195e-6 1673e-6 :=
  RPRing.seeded (model (3 / 2) (1 / 2)) points[23].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 3/2`, `h = 1`. -/
theorem g_3_2_h_1 :
    Record.Seeded (model (3 / 2) 1).spectrum 6 10490 (5 / 2) 6 2241e-8 3345e-8 :=
  RPRing.seeded (model (3 / 2) 1) points[24].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 2`, `h = 1/20`. -/
theorem g_2_h_1_20 :
    Record.Seeded (model 2 (1 / 20)).spectrum 6 7755e-3 (9 / 4) 6 3252e-5 3673e-5 :=
  RPRing.seeded (model 2 (1 / 20)) points[25].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 2`, `h = 1/10`. -/
theorem g_2_h_1_10 :
    Record.Seeded (model 2 (1 / 10)).spectrum 6 872e-2 (9 / 4) 6 2941e-5 3285e-5 :=
  RPRing.seeded (model 2 (1 / 10)) points[26].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 2`, `h = 1/4`. -/
theorem g_2_h_1_4 :
    Record.Seeded (model 2 (1 / 4)).spectrum 6 1745e-2 (9 / 4) 6 1578e-5 1669e-5 :=
  RPRing.seeded (model 2 (1 / 4)) points[27].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 2`, `h = 1/2`. -/
theorem g_2_h_1_2 :
    Record.Seeded (model 2 (1 / 2)).spectrum 6 8976e-2 (9 / 4) 6 3329e-6 3196e-6 :=
  RPRing.seeded (model 2 (1 / 2)) points[28].cert (points_valid _ (List.getElem_mem _))

/-- The longitudinal-field Ising chain, `g = 2`, `h = 1`. -/
theorem g_2_h_1 :
    Record.Seeded (model 2 1).spectrum 6 1073 2 6 1302e-7 3862e-7 :=
  RPRing.seeded (model 2 1) points[29].cert (points_valid _ (List.getElem_mem _))

end Leaderboard.Records.IsingLongitudinal
