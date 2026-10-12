# Routes considered

The requested seed search is finite-volume ED. The ladder is interacting after a Jordan–Wigner transformation; a free-fermion partition function cannot be assumed. A large-field expansion could reach larger rings but requires a new rigorous remainder theorem. Certified MPOs could also reach larger sizes, at substantially higher implementation cost. Neither is silently substituted for the requested ED experiment.

Chosen route: dense complete symmetry blocks at affordable sizes, numerical screening over inverse temperature, and verification of the winning even seed using a polar-decomposition/Weyl residual enclosure. The ladder has leg exchange, rung translations/reflections and global spin-flip symmetry. Commuting involutions permit real symmetry blocks with entries that are integer multiples of square roots of rational numbers. Block dimension sums must equal 4^L. Dense diagonalization costs the sum of cubes of block sizes; residual validation uses BLAS error bounds and directed-rounding scalar arithmetic. Odd n are diagnostic only: the repository theorem explicitly requires even n0.

A ten-rung ring has 1,048,576 states; full dense diagonalization may exceed the per-run time budget even after symmetries. If that happens, report the reached size and the uncompleted search explicitly, without claiming absence of a certificate at ten rungs.
