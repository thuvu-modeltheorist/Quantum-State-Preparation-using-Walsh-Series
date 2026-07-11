# Quantum State Preparation for Max-Cut Using Walsh Series

This README explains how to encode the Max-Cut problem as a quantum
Hamiltonian using Pauli $Z$ operators, and how the Walsh Series Loader (WSL)
can prepare a quantum state whose measurement probabilities are biased
toward good Max-Cut solutions.

---

## Section 1. The Max-Cut Problem and Its Hamiltonian

### 1.1 The Max-Cut problem

Let $G=(V,E)$ be an undirected graph with $V=\{1,\dots,n\}$. Max-Cut asks
for a partition $S,\bar S$ of $V$ maximizing the number of edges crossing
the partition.

Encode a partition as a bitstring $x=(x_1,\dots,x_n)\in\{0,1\}^n$, with
$x_i=1$ iff $i\in S$. Edge $(i,j)$ is cut iff $x_i\neq x_j$, with indicator
$x_i(1-x_j)+x_j(1-x_i)$. So

$$
C(x)=\sum_{(i,j)\in E}\bigl[x_i(1-x_j)+x_j(1-x_i)\bigr]
     =\sum_{(i,j)\in E}\bigl[x_i+x_j-2x_ix_j\bigr],
$$

and the goal is $\max_x C(x)$.

### 1.2 Rewriting in spin variables

Let $z_i=(-1)^{x_i}\in\{+1,-1\}$. Then $z_iz_j=(-1)^{x_i+x_j}$: equal to
$+1$ if $x_i=x_j$, $-1$ if $x_i\neq x_j$. So the cut indicator is
$\frac{1-z_iz_j}{2}$, giving

$$
C(x)=\sum_{(i,j)\in E}\frac{1-z_iz_j}{2}
     =\frac{|E|}{2}-\frac12\sum_{(i,j)\in E}(-1)^{x_i+x_j}.
$$

This **is** already a Walsh expansion of $C$ — see Section 2.

### 1.3 Quantum Hamiltonian

Assign one qubit per vertex, so $|x\rangle=|x_1\cdots x_n\rangle$ represents
one cut. Since $Z_i|x\rangle=(-1)^{x_i}|x\rangle$ and
$Z_iZ_j|x\rangle=(-1)^{x_i+x_j}|x\rangle$, the operator

$$
H_C=\sum_{(i,j)\in E}\frac{I-Z_iZ_j}{2}
   =\frac{|E|}{2}I-\frac12\sum_{(i,j)\in E}Z_iZ_j
$$

satisfies $H_C|x\rangle=C(x)|x\rangle$: it's diagonal in the computational
basis, with the Max-Cut value as its eigenvalue on each $|x\rangle$.

---

## Section 2. Walsh Series and State Preparation

### 2.1 Walsh functions on bitstrings

For $j\in\{0,\dots,2^n-1\}$ with binary expansion $j=\sum_{i=1}^n j_i2^{i-1}$,
and a bitstring $x=(x_1,\dots,x_n)$, define

$$
w_j(x)=(-1)^{j\cdot x}, \qquad j\cdot x = \sum_{i=1}^n j_ix_i.
$$

$w_j(x)$ is the parity of the bits of $x$ selected by $j$'s 1-bits.

> **Relation to Appendix A, Eq. (A1).** The paper's original definition,
> $w_j(x)=(-1)^{\sum_i j_ix_{i-1}}$, is written for a *continuous* point
> $x\in[0,1]$ via its dyadic expansion $x=\sum_i x_i/2^{i+1}$, with an
> index shift ($j_i \leftrightarrow x_{i-1}$) baked in. Once you fix
> $x$ to one of the $N=2^n$ discretization points $\mathcal X_n$, its
> dyadic bits *are* a bitstring, and the two definitions describe the
> same object — the version above is just the cleaner, shift-free
> form appropriate once you're working directly with $\{0,1\}^n$
> (which is the natural domain for Max-Cut) rather than $[0,1]$.

**Example:** $j$ with bits $(j_1,j_2,j_3,j_4)=(1,0,1,0)$ gives
$w_j(x)=(-1)^{x_1+x_3}$.

### 2.2 Walsh operators

$$
\widehat W_j = Z_1^{j_1}\otimes\cdots\otimes Z_n^{j_n}, \qquad
\widehat W_j|x\rangle = w_j(x)|x\rangle.
$$

### 2.3 Walsh expansion of a real-valued function

Any $f:\{0,1\}^n\to\mathbb R$ expands exactly as

$$
f(x)=\sum_{j=0}^{2^n-1}a_j^f w_j(x), \qquad
a_j^f=\frac{1}{2^n}\sum_{x} f(x)w_j(x)
$$

(the discrete Walsh–Hadamard transform of $f$). For Max-Cut, the $2^n$
discretization points of Appendix A correspond exactly to the $2^n$
possible cuts.

### 2.4 A tunable Max-Cut bias function

Let $m=|E|$ and $c(x)=C(x)/m\in[0,1]$ be the normalized cut value. Define

$$
f(x)=1+\lambda\Bigl(c(x)-\tfrac12\Bigr), \qquad 0\le\lambda<2.
$$

A random cut cuts each edge with probability $1/2$, so $c(x)-\tfrac12$
measures how much better or worse a cut is than random; $\lambda$ tunes
the strength of the bias, and $f(x)>0$ for all $x$ is guaranteed exactly
when $\lambda<2$ (the worst case $c=0$ gives $f=1-\lambda/2$).

### 2.5 The bias stays Walsh-sparse

Substituting Section 1.2's expansion of $c(x)-\tfrac12$:

$$
c(x)-\tfrac12 = -\frac{1}{2m}\sum_{(i,j)\in E}(-1)^{x_i+x_j}
\quad\Longrightarrow\quad
f(x)=1-\frac{\lambda}{2m}\sum_{(i,j)\in E}(-1)^{x_i+x_j}.
$$

Each $(-1)^{x_i+x_j}$ is the Walsh function of the index with 1-bits at
positions $i,j$ only. So $f$'s Walsh expansion has exactly $1+m$ nonzero
terms: one constant term, one per edge — sparse and *exact*, no
truncation needed, for any graph.

### 2.6 State preparation target — and how the bias actually arises

The goal is to prepare

$$
|\psi_f\rangle=\frac{1}{\sqrt Z}\sum_x f(x)|x\rangle, \qquad Z=\sum_x|f(x)|^2,
$$

so that measuring $|\psi_f\rangle$ favors good cuts. **This cannot come from
the diagonal unitary $\widehat U_{f,\epsilon_0}=e^{-i\widehat f\epsilon_0}$
alone.** $\widehat U_{f,\epsilon_0}$ is unitary and diagonal in the
computational basis, so applying it to the uniform state
$|s\rangle=H^{\otimes n}|0\rangle$ only changes phases:

$$
\widehat U_{f,\epsilon_0}|s\rangle
=\frac{1}{\sqrt{2^n}}\sum_x e^{-if(x)\epsilon_0}|x\rangle,
$$

every amplitude still has magnitude $1/\sqrt{2^n}$. Measuring this gives
the **uniform** distribution over cuts, not a biased one.

The bias comes from the full **ancilla interference scheme** of the
paper's Fig. 1:
<img width="804" height="384" alt="image" src="https://github.com/user-attachments/assets/1806f97f-06c3-4b1e-95cc-d4e84f439a08" />

1. Prepare $|s\rangle\otimes|+\rangle_{\rm anc}$ ($H^{\otimes n}$ on the
   register, $H$ on one ancilla qubit).
2. Apply $\widehat U_{f,\epsilon_0}$ **controlled by the ancilla**:
   $\;\frac{1}{\sqrt2}\bigl(|s\rangle|0\rangle+e^{-if(x)\epsilon_0}|s\rangle|1\rangle\bigr)$
   (phase applied only in the $|1\rangle_{\rm anc}$ branch).
   <img width="1172" height="776" alt="image" src="https://github.com/user-attachments/assets/697b1215-f9f3-4a21-abd0-d58f9688ecbe" />

4. Apply $H$, then $P=\mathrm{diag}(1,-i)$, to the ancilla. This turns the
   *relative* phase between branches into an *amplitude* difference:
   $$
   \tfrac12(I+e^{-if\epsilon_0})|s\rangle|0\rangle
   -\tfrac{i}{2}(I-e^{-if\epsilon_0})|s\rangle|1\rangle.
   $$
5. **Measure the ancilla; postselect on $|1\rangle$** (repeat-until-success
   if it comes out $|0\rangle$).

The surviving register state is $\propto (I-e^{-if\epsilon_0})|s\rangle$,
whose amplitude on $|x\rangle$ has magnitude $\propto \sin\bigl(f(x)\epsilon_0/2\bigr)$
— monotonically increasing in $f(x)$ for $f\epsilon_0/2\in[0,\pi/2]$, and
$\approx f(x)\epsilon_0/2$ (proportional to $f(x)$, as advertised) only in
the small-$\epsilon_0$ limit. This is the step Section 2.7 of the original
README skipped from "build the diagonal unitary" straight to "prepare a
state whose amplitudes approximate $f(x)$" — the ancilla and postselection
are what actually do that.

### 2.7 Appendix A algorithm, in Max-Cut language

1. One qubit per vertex; $|x\rangle$ = one cut.
2. Choose the bias $f(x)=1+\lambda(c(x)-\tfrac12)$.
3. Compute the (already-exact) Walsh coefficients: $a_0^f=1$, and
   $a_{\{i,j\}}^f=-\lambda/(2m)$ for each edge $(i,j)$, else $0$.
4. Build $\widehat f=\sum_j a_j^f\widehat W_j$ and $\widehat U_{f,\epsilon_0}=e^{-i\widehat f\epsilon_0}$.
5. Since all $\widehat W_j$ commute (tensor products of $Z$'s),
   $\widehat U_{f,\epsilon_0}=\prod_j\exp(-ia_j^f\epsilon_0\widehat W_j)$ —
   one factor per Walsh term.
6. Wrap the **controlled** version of this product in the ancilla scheme
   of Section 2.6 to get real amplitude bias, not just phase.
7. Measure, or hand the resulting state to QAOA as its initial state.

### 2.8 Decomposition into Walsh operators

For the linear Max-Cut bias,
$\widehat f = I-\frac{\lambda}{2m}\sum_{(i,j)\in E}Z_iZ_j$, so

$$
\widehat U_{f,\epsilon_0}=\exp\!\Bigl(-i\epsilon_0 I+i\tfrac{\lambda\epsilon_0}{2m}\sum_{(i,j)\in E}Z_iZ_j\Bigr).
$$

Dropping the global phase $e^{-i\epsilon_0 I}$ (note: this phase must
still be applied *inside* the ancilla-controlled block per Section 2.6 —
it's only "global" with respect to the register, not irrelevant overall,
since it's what makes the $\sin(\cdot)$ argument include the constant
term $a_0$):

$$
\widehat U_{f,\epsilon_0}\sim\prod_{(i,j)\in E}\exp\!\Bigl(i\tfrac{\lambda\epsilon_0}{2m}Z_iZ_j\Bigr).
$$

One $ZZ$-rotation per edge.

### 2.9 Circuit for one edge term

For edge $(i,j)$, implement $\exp(i\theta Z_iZ_j)$ with $\theta=\lambda\epsilon_0/(2m)$:

```text
CNOT(i -> j)
RZ(-2*theta) on qubit j      # convention: RZ(phi) = diag(e^{-i phi/2}, e^{i phi/2})
CNOT(i -> j)
```

Sign conventions vary by software package — the structural point is:
compute parity via CNOTs, rotate, uncompute the parity.

### 2.10 Complete pipeline

```text
Input: graph G=(V,E), n vertices, m edges

1.  One qubit per vertex; bitstring x in {0,1}^n encodes a cut.
2.  C(x) = sum_{(i,j) in E} [x_i + x_j - 2 x_i x_j];  c(x) = C(x)/m.
3.  f(x) = 1 + lambda*(c(x) - 1/2),  0 <= lambda < 2.
4.  f(x) = 1 - (lambda/(2m)) * sum_{(i,j) in E} (-1)^{x_i+x_j}   [exact, sparse]
5.  Register: H^{\otimes n}.  Ancilla: H.
6.  Controlled-U_{f,eps0}: constant phase (a_0 term) + one CNOT-RZ-CNOT
    block per edge, all controlled by the ancilla.
7.  H, then P = diag(1,-i), on the ancilla.
8.  Measure the ancilla; postselect on |1>  (repeat-until-success on |0>).
9.  Surviving register state has amplitude ~ sin(f(x)*eps0/2) on |x> --
    biased toward high-cut x.
10. Measure directly, or use as QAOA's initial state.
```

---

## Summary

$$
C(x)=\frac{|E|}{2}-\frac12\sum_{(i,j)\in E}(-1)^{x_i+x_j}
\quad\longleftrightarrow\quad
H_C=\frac{|E|}{2}I-\frac12\sum_{(i,j)\in E}Z_iZ_j.
$$

Every edge $\leftrightarrow$ one Walsh operator $Z_iZ_j$. The tunable bias
$f(x)=1+\lambda(c(x)-\tfrac12)$ stays exactly Walsh-sparse ($1+|E|$ terms).
The Walsh Series Loader's **ancilla interference + postselection** — not
the diagonal unitary by itself — is what converts this sparse phase
structure into a real amplitude bias toward good cuts, at rate
$\sin(f(x)\epsilon_0/2)$.
