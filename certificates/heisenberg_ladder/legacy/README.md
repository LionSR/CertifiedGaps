# Original validated snapshot

These files are byte-for-byte copies of the validated ThermalMPOGap run at
repository base ef0a3e6, before upstream renamed the core to CertifiedGaps.
The original hashes are preserved, not rewritten. The current certificate and
independent replay audit live in ../results and were generated on current main.

To reproduce this historical snapshot, check out ef0a3e6 in a separate clone,
copy source/ into certificates/heisenberg_ladder/, copy these results/ into that
certificate directory, and run its test/runtests.jl with --project=julia.
Core source files and manifests are supplied by that Git commit. The archived
source/ files are not intended to run in place from this nested directory.

Historical numerical scans retain their original provenance. This archive
covers the rigorously validated J_perp=4 snapshot, not every intermediate
version used during exploratory numerical optimization.

The archived producer used 10^-12 seed rounding. Small variations in fresh
residual/orthogonality allowances can exceed that tight margin, even though
they do not invalidate the original computed enclosure. Current-head
certification uses outward 10^-8 rounding and passes strict fresh replay.
The historical files themselves are preserved without rewriting their hashes.
