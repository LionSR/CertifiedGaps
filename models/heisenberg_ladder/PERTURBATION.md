# Arbitrary-length strong-rung norm sandwich: a rigorous negative result

This is a negative result for a particular bounding route, not for the ladder gap.
It applies to every even seed length n >= 2, every beta > 0, and 0 < J <= 4.
It explains why increasing n cannot rescue a naive independent-rung perturbation.

Let H0 contain only the rungs, V contain the 2L leg bonds, and
z(beta) = exp(3 beta J/4) + 3 exp(-beta J/4). Each leg bond h has eigenvalues
-3/4 and 1/4. Center it: ||h + I/4|| = 1/2. Therefore

    ||V + L I/2|| <= L,
    |log Z_L(beta) - [L log z(beta) + beta L/2]| <= beta L.

Weyl's eigenvalue inequality proves the second line; it is valid for all L,
without ED, expansions, or unverified spectral input. Bulk shifts cancel from
both purities. Applying this sandwich independently to the three traces gives

    eps_s = 1 - exp(-4 beta n),
    eps_t = 1 - exp(-8 beta n) p(beta)^(2n),
    p(beta) = (1 + 3 q^2)/(1 + 3 q)^2, q = exp(-beta J).

These are rigorous upper bounds on the seed defects, but cannot enter the
admissible recurrence region for any n or beta in the stated coupling range:

* If beta n >= 1/8, then (1-eps_s)^2 = exp(-8 beta n) <= exp(-1) < 1/2.
  The strict last inequality follows already from exp(1) > 1+1+1/2 > 2.
* If beta n < 1/8, then beta J <= 4 beta < 1/(2n) <= 1/4. Hence
  3/4 < 1-beta J <= q <= 1, and p <= 4/(1+9/4)^2 = 64/169 < 1/2.
  In particular 1-eps_t < 1/2.

In both cases A=(1-eps_s)^2(1-eps_t)<1/2. Consequently z=1-A>1/2;
f(z)=z^2/[2(1-z)^2]>1/2. The second component of the next doubling is
also at least 1/2, and the order-preserving recurrence cannot converge to zero.
This proves failure of this sandwich even beyond ED and even after optimizing
beta continuously. It does not exclude a linked-cluster wrapping bound, a
transfer-contraction proof, a better perturbative estimate, or the actual seed.

## Useful local continuation instead

At a certified coupling J*, set H_L(J)-H_L(J*) = (J-J*) W_L.
The centered rung sum satisfies ||W_L + L I/4|| <= L/2. For |J-J*| <= delta,
subtracting the scalar shift gives |change in log Z_L| <= beta delta L/2.
Thus from a seed (s,t) at n,beta, valid enlarged upper defects are

    s_delta = 1 - (1-s) exp(-2 beta n delta),
    t_delta = 1 - (1-t) exp(-4 beta n delta).

Evaluating these expressions outward and closing the recurrence certifies an
entire coupling interval, without assuming monotonicity of the ladder gap or
its seed defects. `check.jl` recomputes this continuation for δ = 1/10000.

## Smaller-seed identity checked

Directly from partition functions,

    P_tau(2n,beta) = P_s(n,2beta) P_tau(n,beta)^2 / P_s(n,beta)^2.

Since RP gives P_s(n,beta)<=1 for even n, the product in the numerator is a
valid lower bound. This does not remove Z_2n(2beta) unless P_s(n,2beta) is
bounded by additional structure. We therefore use the direct complete trace,
and do not claim a smaller largest ring from an unproved substitute bound.
