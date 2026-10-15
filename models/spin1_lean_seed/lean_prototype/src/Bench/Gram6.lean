import SeedPrototype.Packed
set_option maxRecDepth 100000
set_option maxHeartbeats 0
open SeedPrototype
-- Four orthogonal product-state trials, a nontrivial full 729-lane Gram computation.
example : gram 128 (3^6) [1, 1 <<< (128*neel 6),
    1 <<< (128*(3^6-1-neel 6)), 1 <<< (128*(3^6-1))] =
    [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]] := by decide +kernel
