# Congruence subgroups of Bianchi groups

Two SageMath scripts for congruence subgroups `Γ₀(a)` of the five Euclidean Bianchi
groups, at prime-power level `a = πᵉ`:

| script | computes |
|---|---|
| `bianchi_abelianisation.sage` | the abelianisation of `Γ₀(a)` in `SL₂(O_d)`, `PSL`, `GL₂(O_d)` or `PGL` |
| `bianchi_hecke.sage` | the Hecke action on `H₁`, over the free part or mod `p` |

Here `d ∈ {-1, -2, -3, -7, -11}`, `O_d` is the ring of integers of `Q(√d)`, and `O_d = Z[x]`
with `x = √d` for `d = -1, -2` and `x = (1 + √d)/2` for `d = -3, -7, -11`.

Both scripts are interactive: run one and answer its questions.

These scripts extend the work done in my master's thesis "Torsion in the homology of Bianchi groups and Galois representations", for which I wrote SageMath programs to compute abelianisations of Gamma0 type congruence subgroups of the (SL_2) euclidian Bianchi groups, and additionally to compute the action of the Hecke algebra on the first homology groups. 

Since completing my master's degree, in the few moments of time I have to think about math, I was able to come up with ways to extend my algorithms to include computations at prime power levels. I don't believe the algorithm I have in place here can easily extend to all levels. Regardless, this still expands the computations of abelianisations from just SL_2 and for prime level, to now include PSL_2, GL_2 and PGL_2, in addition to including prime power levels.

The computation of the Hecke algebra also extends to these new cases. 

I hope there is someone that finds this code useful! I thoroughly enjoyed revisiting my work and expanding on it (and making it significantly more efficient!) please use these scripts, tweak the code (I'm sure lots of it is still suboptimal), use it to do interesting math!

-Emre
