# Open Problems & Next Steps

**Project:** Hematite nonrelativistic magnon band structure (Phase 1)  
**Reference:** Hoyer et al., arXiv:2503.11623 (2025)  
**Date:** 2026-05-06

---

## What works

| Check | Result |
|-------|--------|
| Crystal structure (R3̄c, 4 sublattices A–D) | ✓ |
| All 6 13th-neighbour bond distances equal 6.4228 Å | ✓ |
| 8×8 BdG kernel Hermiticity: max \|H − H†\| = 0 | ✓ |
| Block-diagonal structure: off-block max = 0 | ✓ |
| Acoustic Goldstone mode at Γ approaches 0 | ✓ |
| Band dispersion topology matches Fig. 5(a) qualitatively | ✓ |

---

## Problem 1 — Energy scale off by factor ~10

### Symptom
The optical magnon branch sits at ~97 meV in the current output. Hoyer et al. Fig. 5(a) shows it at ~9–10 meV. The ratio is ≈ √(97) ≈ 9.85 — suspiciously close to a square-root factor.

### Root cause hypothesis
The physical magnon energies of a bosonic BdG system are obtained from the **eigenvalues of** $GH_{\mathbf{k}}$, where $G = \mathrm{diag}(I_n, -I_n)$ is the BdG metric:

$$GH_{\mathbf{k}}\, \mathbf{u} = \varepsilon\, \mathbf{u}$$

The previous Colpa implementation computed instead $\sqrt{\lambda}$ where $\lambda$ are the eigenvalues of $K^\dagger G K$ (with $H = KK^\dagger$ the Cholesky factorisation). Algebraically these are the **same** quantity ($\lambda = \varepsilon$), so taking an extra square root gives $\sqrt{\varepsilon}$ — numerically landing near 10 meV only because $\sqrt{97\,\mathrm{meV}} \approx 9.85\,\mathrm{meV}$.

### What to verify next
1. **Read Eq. (8)–(9) and surrounding text of Hoyer et al. carefully**: confirm whether the on-site coefficient $A_{\mathrm{ISO}} = S(J_1 - 3J_2 + \ldots) \approx 97\,\mathrm{meV}$ is the correct value or whether it should already equal ~10 meV at the optical mode.
2. **Cross-check with a minimal model**: for a 2-sublattice Heisenberg AFM with a single exchange $J > 0$, the magnon dispersion is $\varepsilon_k = 2JSz\sqrt{1-\gamma_k^2}$. Plug in $J_3 = 5.452\,\mathrm{meV}$, $S = 5/2$, $z = 3$ (3rd-neighbour coordination) and evaluate at the zone boundary to compare with the paper.
3. **Check normalization convention**: some spin-wave papers factor out $S$ from the BdG kernel so that $H \to H/S$ or include a factor of $2S$ in the Holstein-Primakoff expansion. If Hoyer et al. normalize differently, the on-site coefficient changes by $\mathcal{O}(S)$.
4. **Use the original Colpa Cholesky method with sqrt** as a pragmatic workaround, with a clear note that the convention differs from the GH diagonalization by $\varepsilon_{\mathrm{Colpa}} = \sqrt{\varepsilon_{\mathrm{GH}}}$.

---

## Problem 2 — Altermagnetic spin splitting is identically zero

### Symptom
The altermagnetic splitting $\Delta\varepsilon_{\mathrm{ac}} = |\varepsilon_1 - \varepsilon_3|$ (lowest spin-$+\tfrac{1}{2}$ acoustic minus lowest spin-$-\tfrac{1}{2}$ acoustic) is exactly zero at **every** k-point on every path tested (Γ–B, Γ–X, and generic off-symmetry points). The splitting remains at floating-point noise level (~10⁻¹² meV).

### Root cause (confirmed analytically)
The two spin-sector 4×4 BdG blocks are related by a **unitary permutation**:

$$H_2 = P\, H_1\, P^\dagger, \qquad P = \begin{pmatrix} 0&1&0&0\\1&0&0&0\\0&0&0&1\\0&0&1&0 \end{pmatrix}$$

where $P$ swaps $(a^\dagger \leftrightarrow d^\dagger)$ within the particle sector and $(b_{-\mathbf{k}} \leftrightarrow c_{-\mathbf{k}})$ within the hole sector. This permutation **commutes with the BdG metric** $G_4 = \mathrm{diag}(+1,+1,-1,-1)$:

$$P G_4 = G_4 P$$

Therefore $G_4 H_2 = P\,(G_4 H_1)\,P^\dagger$, which is a similarity transformation. $G_4 H_1$ and $G_4 H_2$ share identical eigenvalues, and the two spin blocks yield the same magnon energies at every $\mathbf{k}$.

This is an **exact mathematical identity**, independent of the values of $\Delta_{\mathbf{k}}$, $D_{\mathbf{k}}$, etc. No choice of k-path or diagonalization method can produce non-zero splitting while this symmetry holds.

### Physical interpretation
The $P$-symmetry survives because the H_AM matrix adds $\Delta_{\mathbf{k}}$ to the $(a^\dagger, d^\dagger)$ particle-particle coupling **and** the $(b_{-\mathbf{k}}, c_{-\mathbf{k}})$ hole-hole coupling in block 1, and symmetrically to the corresponding positions in block 2. This means the altermagnetic anisotropy enters both spin sectors in an identical way, leaving the spectra degenerate.

### What to verify next
1. **Check the exact positions of $\Delta_{\mathbf{k}}$ in Eq. (8) of Hoyer et al.** — specifically whether $\Delta_{\mathbf{k}}$ appears at positions (2,3)/(3,2) and (6,7)/(7,6) in the 8×8 matrix (i.e., in the $b^\dagger$–$c^\dagger$ and $b_{-\mathbf{k}}$–$c_{-\mathbf{k}}$ channels), or whether those entries should be zero (since the altermagnetic exchange is only on A–D bonds, not B–C bonds).
2. **If positions (2,3)/(3,2)/(6,7)/(7,6) should be zero**: removing them breaks the $P$-symmetry. Specifically, $H_1[3,4]$ would become $D^*$ (without $\Delta^*$) while $H_2[3,4]$ remains $D + \Delta$, making the two blocks non-isospectral and producing genuine spin splitting.
3. **Check whether nonrelativistic splitting is physical at all**: in a fully collinear antiferromagnet without spin-orbit coupling, time-reversal symmetry protects spin degeneracy everywhere in the BZ. If this is the case, the splitting in Hoyer et al. Fig. 5(c) may only appear in the relativistic (EAP or WFP) phase with DMI/SIA, and the nonrelativistic spin coloring in Fig. 5(a) is a spin-sector label rather than an energy splitting.

---

## Suggested resolution order

| Priority | Action |
|----------|--------|
| 1 (quick) | Open Hoyer et al. Eq. (8) and confirm exact positions of $\Delta_{\mathbf{k}}$ in the 8×8 H_AM matrix |
| 2 (quick) | Confirm whether Fig. 5(c) shows nonrelativistic or EAP-phase splitting |
| 3 (medium) | Cross-check energy scale with a 2-sublattice minimal model |
| 4 (medium) | If $\Delta_{\mathbf{k}}$ at (2,3)/(6,7) should be zero, remove and rerun |
| 5 (medium) | If energy scale is confirmed wrong, add normalization factor and rerun |
| 6 (long) | Move to Phase 2: add SIA + DMI to `hamiltonian_eap.jl` |
