import SeedPrototype.Packed
namespace SeedPrototype

/-- Reference radix encoding, used to state carry conditions without bitwise notation. -/
def radixPack (b : Nat) : List Nat → Nat
  | [] => 0
  | a::xs => a + b*radixPack b xs

/-- Reference extraction by repeated Euclidean division. -/
def radixLane (b : Nat) : Nat → Nat → Nat
  | 0, x => x%b
  | i+1, x => radixLane b i (x/b)

theorem radixPack_head (b a : Nat) (xs : List Nat) (ha : a < b) :
    radixPack b (a::xs) % b = a := by
  simp [radixPack, Nat.add_mod, Nat.mod_eq_of_lt ha]

theorem radixPack_tail (b a : Nat) (xs : List Nat) (ha : a < b) :
    radixPack b (a::xs) / b = radixPack b xs := by
  rw [radixPack, Nat.add_mul_div_left _ _ (by omega), Nat.div_eq_of_lt ha]
  simp

/-- Lane extraction is exact whenever all input digits are below the radix. -/
theorem radixPack_extract (b : Nat) (xs : List Nat)
    (hx : ∀ a ∈ xs, a < b) (i : Nat) (hi : i < xs.length) :
    radixLane b i (radixPack b xs) = xs[i] := by
  induction xs generalizing i with
  | nil => simp at hi
  | cons a xs ih =>
    cases i with
    | zero => simpa [radixLane] using radixPack_head b a xs (hx a (by simp))
    | succ i =>
      rw [radixLane, radixPack_tail b a xs (hx a (by simp))]
      exact ih (fun a ha => hx a (by simp [ha])) i (by simpa using hi)

/-- Local no-carry addition, including contributions from arbitrary higher lanes. -/
theorem add_noCarry (b x y : Nat) (h : x%b + y%b < b) :
    (x+y)%b = x%b+y%b := by
  rw [Nat.add_mod, Nat.mod_eq_of_lt h]

/-- The executable left-shift encoder equals the mathematical radix encoder. -/
theorem pack_eq_radixPack (w : Nat) (xs : List Nat) :
    pack w xs = radixPack (2^w) xs := by
  induction xs with
  | nil => rfl
  | cons a xs ih =>
    change a + (pack w xs <<< w) = a + 2^w * radixPack (2^w) xs
    rw [Nat.shiftLeft_eq, ih, Nat.mul_comm]

end SeedPrototype
