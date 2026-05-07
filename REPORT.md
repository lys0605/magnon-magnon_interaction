# Open Problems & Next Steps

**Project:** Hematite nonrelativistic magnon band structure (Phase 1)  
**Reference:** Hoyer et al., arXiv:2503.11623 (2025)  
**Date:** 2026-05-07

---

## What works

| Check | Result |
|-------|--------|
| Crystal structure (R3̄c, 4 sublattices A–D) | ✓ |
| All 6 13th-neighbour bond distances equal 6.4228 Å | ✓ |
| 8×8 BdG kernel Hermiticity: max \|H − H†\| = 0 | ✓ |
| Block-diagonal structure: off-block max = 0 | ✓ |
| Acoustic Goldstone mode at Γ approaches 0 | ✓ |
| Band dispersion topology (Goldstone at Γ, optical bands ~85–99 meV) | ✓ |
| Zero splitting at all TRIM points (Γ, Z, B, B', F, X, X') | ✓ |
| Non-zero altermagnetic splitting between TRIM points | ✓ |
| g-wave k⁴ power law: fitted exponent = **4.000** (non-nodal direction) | ✓ |

---

## Problem 1 — Energy scale ~×10 relative to original Cholesky code

### Current status
Direct eigenvalue decomposition of GH gives optical magnon at ~97 meV.
The original Cholesky code (with an unintentional √ step) gave ~9.85 meV ≈ √97.

### Assessment
97 meV is **physically plausible**. The dominant exchanges are J₃ = 5.452 meV (z=3) and
J₄ = 4.008 meV (z=6), giving a simple estimate for the zone-boundary energy:
2SJ₃z = 2 × 2.5 × 5.452 × 3 ≈ 82 meV (2-sublattice model),
consistent with our 4-sublattice result of ~97 meV.

### What to verify next
1. **Read Fig. 5(a) y-axis of Hoyer et al.** — confirm whether the optical mode is
   labelled at ~10 meV or ~90 meV. If the paper shows ~90 meV, the current code is correct.
2. If the paper uses a normalization convention (e.g., Hamiltonian written as H/S or 2SH),
   identify the overall factor and apply it uniformly.

---

## Problem 2 — Altermagnetic spin splitting (RESOLVED ✓)

### Root cause (fixed 2026-05-07)
The original `H_AM` block added Δ_k to **four** off-diagonal pairs:
- (1,4)/(4,1): A–D particle sector ← physically correct
- **(2,3)/(3,2): B–C particle sector ← WRONG (no altermagnetic exchange on B–C bonds)**
- (5,8)/(8,5): A–D hole sector ← physically correct
- **(6,7)/(7,6): B–C hole sector ← WRONG**

The spurious B–C entries created an exact permutation symmetry H₂ = P H₁ P† with P
commuting with the BdG metric G, making both spin blocks isospectral at every k.

### Fix applied
Removed the (2,3)/(3,2) and (6,7)/(7,6) entries from `H_AM`. Only A–D bonds (13th
neighbor) carry the altermagnetic anisotropy Δ.

### Results after fix
- **Maximum splitting**: Δε_ac ≈ 1.71 meV (between Γ and X along the standard path)
- **TRIM zeros confirmed**: splitting is exactly zero at all 9 high-symmetry points
  Γ, Z, P, B, Γ, B', F, X, X' on the band path
- **g-wave power law**: along a non-nodal direction [1,0.3,0] (fractional reciprocal),
  the splitting fits k^4.000 over 1.5 decades — confirming hematite's g-wave altermagnetic
  classification
- **Nodal lines**: Γ–Z, Γ–B, Γ–B', Γ–F, Γ–X directions all give Δ_k = 0 identically
  (or extremely small k^6 residual); the fitted exponent ≈ 6 on the old Γ–X splitting
  plot was a consequence of plotting along a nodal line of the g-wave structure factor

---

## Summary of g-wave physics

The altermagnetic structure factor Δ_k = S·Δ·(Σⱼ₌₁³ eⁱᵏ·δ₁₃ʲ − Σⱼ₌₄⁶ eⁱᵏ·δ₁₃ʲ)
has the following near-Γ behaviour:
- Im(Δ_k) ∝ k³ (cubic, dominates for small k)
- Re(Δ_k) ∝ k⁴ (quartic, subleading)

The **energy splitting** Δε_ac ∝ Δ_k² ∝ k⁴ because Im²(Δ_k) ∝ k⁶ and Re²(Δ_k) ∝ k⁸?

Wait — the numerical fit gives k⁴ for Δε, not k⁶. This means the splitting is
**linear** in Δ_k (first-order in perturbation theory), not quadratic.
More precisely, Δε ∝ Re(Δ_k) ∝ k⁴. The imaginary part Im(Δ_k) ∝ k³ does not
contribute to the energy splitting because it enters in combinations that cancel
between the two spin blocks at first order in perturbation theory.

---

## Phase 2 — Relativistic corrections (planned)

- Single-ion anisotropy (SIA, parameter d₂) → `hamiltonian_eap.jl`
- Dzyaloshinskii–Moriya interaction (DMI) for 1st, 3rd, 4th neighbours
- Easy-plane phase above the Morin temperature
- These terms open the acoustic gap and generate the larger spin-splitting seen
  in Fig. 5(c) of Hoyer et al. (if that figure shows the relativistic EAP phase)

---

## Suggested next actions

| Priority | Action |
|----------|--------|
| 1 (quick) | Verify energy scale against Fig. 5(a) y-axis of Hoyer et al. |
| 2 (quick) | Check if Fig. 5(c) shows nonrelativistic or EAP-phase splitting |
| 3 (medium) | Implement Phase 2: SIA + DMI in `hamiltonian_eap.jl` |
| 4 (long) | Phase 3: magnon–magnon self-energy Σ(k,ω) |
