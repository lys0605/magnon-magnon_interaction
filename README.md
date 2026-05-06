# Magnon–Magnon Interaction: Renormalization & Topology

[![Julia v1.10](https://img.shields.io/badge/Julia-v1.10-9558B2?logo=julia&logoColor=white)](https://julialang.org)
[![CairoMakie](https://img.shields.io/badge/Plotting-CairoMakie-389826?logo=julia&logoColor=white)](https://docs.makie.org/stable/)
[![JLD2](https://img.shields.io/badge/Data-JLD2%2FHDF5-blue)](https://github.com/JuliaIO/JLD2.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Status: Active](https://img.shields.io/badge/Status-Active-brightgreen)]()

Research code for studying **magnon–magnon interaction-induced renormalization** of magnon bands and the stability of topological and quantum-geometric properties in altermagnetic spin-wave systems.

The primary material system is **hematite (α-Fe₂O₃)** — a four-sublattice g-wave altermagnet modelled via linear spin-wave theory in the bosonic BdG framework, following Hoyer *et al.* (2025). The second focus (planned) is magnon–magnon interactions in **topological magnon insulators**, following Habel *et al.* (2023).

---

## Background

| Topic | Material | Reference |
|-------|----------|-----------|
| Altermagnetic magnon splitting | α-Fe₂O₃ (hematite) | Hoyer *et al.*, [arXiv:2503.11623](https://arxiv.org/abs/2503.11623) (2025) |
| Breakdown of chiral edge modes | Topological magnon insulator | Habel *et al.*, [arXiv:2308.03168](https://arxiv.org/abs/2308.03168) (2023) |

---

## Project Structure

```
.
├── README.md
├── REPORT.md                          # Open problems and next steps
│
└── magnon/                            # Main Julia project
    ├── Project.toml                   # Julia package manifest
    ├── Manifest.toml
    ├── CLAUDE.md                      # AI-assistant workflow notes
    │
    ├── src/                           # Core library (include-based)
    │   ├── math_utils.jl              # BdG solver, k-path builder, meshgrid
    │   ├── plot_utils.jl              # CairoMakie theme + band-plot helpers
    │   ├── hematite_lattice.jl        # R3̄c crystal structure, bond vectors, BZ
    │   └── hematite_hamiltonian.jl    # 8×8 BdG kernel, spin polarization
    │
    ├── compute_hematite_bands.jl      # Compute bands → data/
    ├── test_sanity.jl                 # Quick physics checks
    │
    ├── data/                          # Computed JLD2 files (gitignored)
    └── figures/
        ├── fig_hematite_bands.jl      # Dispersion colored by ⟨s^N̂⟩
        ├── fig_hematite_splitting.jl  # Acoustic splitting vs k
        └── output/                    # Generated PNGs (gitignored)
```

---

## Physics Overview

### Phase 1 — Hematite nonrelativistic magnon bands (in progress)

Hematite (α-Fe₂O₃, space group R$\bar{3}$c, $S = 5/2$) is a g-wave altermagnet with a four-sublattice Néel order. Linear spin-wave theory via Holstein–Primakoff gives an **8×8 bosonic BdG kernel**:

$$H_{\mathbf{k}} = H_{\mathbf{k}}^{\mathrm{ISO}} + H_{\mathbf{k}}^{\mathrm{AM}}$$

- $H^{\mathrm{ISO}}$: isotropic Heisenberg exchange through 13th neighbours
- $H^{\mathrm{AM}}$: altermagnetic term from two inequivalent 13th-neighbour A–D bond types (exchange anisotropy $\Delta$)

Physical energies are extracted from the positive eigenvalues of $GH_{\mathbf{k}}$ where $G = \mathrm{diag}(+1,+1,+1,+1,-1,-1,-1,-1)$ is the bosonic BdG metric. In the nonrelativistic collinear limit $H_{\mathbf{k}}$ is block-diagonal in the two spin sectors $\langle s^{\hat{N}}\rangle = \pm\tfrac{1}{2}$, which should yield the characteristic altermagnetic **g-wave ($k^4$) splitting** near $\Gamma$.

> **Current status:** Lattice, bond vectors, BdG kernel, and k-path are implemented and pass structural sanity checks (Hermiticity, block structure, bond distances). Two open problems are being resolved — see [`REPORT.md`](REPORT.md).

### Phase 2 — Relativistic corrections (planned)

- Single-ion anisotropy (SIA, parameter $d_2$) → `hamiltonian_eap.jl`
- Dzyaloshinskii–Moriya interaction (DMI) for 1st, 3rd, 4th neighbours
- Easy-plane phase above the Morin temperature

### Phase 3 — Magnon–magnon interactions (planned)

- $1/S$ correction beyond linear spin-wave theory
- Self-energy $\Sigma(\mathbf{k},\omega)$ from four-magnon vertices
- Renormalization of topological edge modes

---

## Getting Started

### Requirements

- Julia ≥ 1.10

### Setup

```julia
# From the magnon/ directory
using Pkg; Pkg.activate("."); Pkg.instantiate()
```

### Workflow

```julia
include("test_sanity.jl")                    # verify structure first
include("compute_hematite_bands.jl")         # writes data/*.jld2
include("figures/fig_hematite_bands.jl")     # dispersion plot
include("figures/fig_hematite_splitting.jl") # splitting plot
```

---

## Key Parameters (Tab. 1, Hoyer et al. 2025)

| Symbol | Value | Unit | Source |
|--------|-------|------|--------|
| $S$ | 5/2 | — | — |
| $J_1$ | −0.982 | meV | Neutron scattering |
| $J_2$ | −0.169 | meV | Neutron scattering |
| $J_3$ | +5.452 | meV | Neutron scattering |
| $J_4$ | +4.008 | meV | Neutron scattering |
| $J_5$ | +0.337 | meV | Neutron scattering |
| $J_{13}$ | +0.143 | meV | DFT |
| $\Delta$ | +0.148 | meV | DFT (altermagnetic exchange anisotropy) |

---

## License

MIT — see [LICENSE](LICENSE).
