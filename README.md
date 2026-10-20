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

bianchi_abelianisation.sage asks the following questions:

1) Which Bianchi group [possible inputs are -1,-2,-3,-7,-11] (the only supported Bianchi groups are the euclidian ones)
2) Which rational prime divides the norm of the level [input a rational prime] (if the rational prime splits/ramifies in the number field, the program selects an appropriate prime in O)
3) What is the exponent of the level [input a positive integer] (This determines the level: the prime from question (2) raised to the exponent entered here.)
4) Which group? [input SL, PSL, GL, PGL] (not case sensitive)

It then computes the abelianisation

Bianchi_hecke.sage starts with the same 4 questions, and then continues:

5) Which coefficients for the Hecke action? [input: "free" for complex coefficients, prime number "p" for coefficients in algebraic closure of F_p]
6) How many Hecke operators should be computed [input: a positive integer] (the program always includes hecke operators at the prime dividing the level, and at the even prime(s), For SL, PSL, it also includes a "Hecke operator at a unit")
7) Any particular rational primes you want Hecke operators for? [input: rational primes] (If you're interested in specific eigenvalues you can therefore input them here)

- The program also computes the action of the Fricke involution, and for complex-conjugation invariant levels, the action of complex conjugation.

- It returns:
  1) the dimension of the space (over QQ or F_p)
  2) Characteristic polynomials of the Hecke operators (factored into irreducibles over the base field)
  3) The Hecke eigenspaces, with dimension, generalized dimension, the eigenvalues at each operator and whether it is one of the two trivial classes described in my thesis.
  4) Determines if Fricke involution preserves the eigenclass, and if so, gives the eigenvalue.
  5) Determines if complex conjugation preserves the eigenclass, if so, gives the eigenvalue, if not, identifies which eigenclass it pairs with.

These scripts extend the work done in my master's thesis "Torsion in the homology of Bianchi groups and Galois representations", for which I wrote SageMath programs to compute abelianisations of Gamma0 type congruence subgroups of the (SL_2) euclidian Bianchi groups, and additionally to compute the action of the Hecke algebra on the first homology groups. 

Since completing my master's degree, in the few moments of time I have to think about math, I was able to come up with ways to extend my algorithms to include computations at prime power levels. I don't believe the algorithm I have in place here can easily extend to all levels. Regardless, this still expands the computations of abelianisations from just SL_2 and for prime level, to now include PSL_2, GL_2 and PGL_2, in addition to including prime power levels.

The computation of the Hecke algebra also extends to these new cases. 

I hope there is someone that finds this code useful! I thoroughly enjoyed revisiting my work and expanding on it (and making it significantly more efficient!) please use these scripts, tweak the code (I'm sure lots of it is still suboptimal), use it to do interesting math!

For contact please reach out to my email:
emrek-at-hotmail.nl

-Emre
