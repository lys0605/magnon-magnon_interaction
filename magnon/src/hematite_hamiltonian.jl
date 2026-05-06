isdefined(Main, :_HEMATITE_HAMILTONIAN_LOADED) && return
const _HEMATITE_HAMILTONIAN_LOADED = true

include(joinpath(@__DIR__, "math_utils.jl"))
include(joinpath(@__DIR__, "hematite_lattice.jl"))

using LinearAlgebra

# ---------------------------------------------------------------------------
# Hematite magnon Hamiltonian — nonrelativistic limit
#
# Four-sublattice linear spin-wave theory for α-Fe₂O₃.
# Reference: Hoyer et al., arXiv:2503.11623 (2025), Sec. III.
#
# Hamiltonian: H = H_ISO + H_AM
#   H_ISO = isotropic Heisenberg exchange (J₁–J₅, J₁₃)
#   H_AM  = altermagnetic term (±Δ on 13th-neighbor bonds)
#
# Basis: Ψ†_k = (a†_k, b†_k, c†_k, d†_k, a_{-k}, b_{-k}, c_{-k}, d_{-k})
# where a/b/c/d are Holstein-Primakoff bosons on sublattices A/B/C/D.
#
# Parameters (Tab. 1 of Hoyer et al.):
#   S   = 5/2
#   J₁  = −0.982 meV   (nearest-neighbor,  A-B type)
#   J₂  = −0.169 meV   (2nd neighbor,      A-D type)
#   J₃  =  5.452 meV   (3rd neighbor,      A-B type)
#   J₄  =  4.008 meV   (4th neighbor,      A-C type)
#   J₅  =  0.337 meV   (5th neighbor,      B-C type)
#   J₁₃ =  0.143 meV   (13th neighbor,     A-D type, isotropic part)
#   Δ   =  0.148 meV   (altermagnetic exchange, anisotropic part)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Default parameters (all in meV, spin dimensionless)
# ---------------------------------------------------------------------------
struct HematiteParams
    S   :: Float64   # spin quantum number
    J₁  :: Float64   # meV
    J₂  :: Float64
    J₃  :: Float64
    J₄  :: Float64
    J₅  :: Float64
    J₁₃ :: Float64
    Δ   :: Float64   # altermagnetic exchange
end

"""
    default_params() -> HematiteParams

Return the nonrelativistic parameters from Tab. 1 of Hoyer et al. (2025).
J values from neutron scattering (Samuelsen 1970); J₁₃ and Δ from DFT.
"""
function default_params()
    HematiteParams(
        5/2,      # S
       -0.982,    # J₁
       -0.169,    # J₂
        5.452,    # J₃
        4.008,    # J₄
        0.337,    # J₅
        0.143,    # J₁₃
        0.148,    # Δ
    )
end

# ---------------------------------------------------------------------------
# Fourier-space structure factors  (Eqs. 7 of Hoyer et al.)
# ---------------------------------------------------------------------------

"""
    structure_factors(k, p) -> (A_ISO, B_k, C_k, D_k, Δ_k)

Compute the momentum-dependent coefficients of the harmonic Hamiltonian.

# Arguments
- `k`: 3-vector of crystal momentum (Cartesian, Å⁻¹)
- `p`: HematiteParams

# Returns (all Complex{Float64} unless noted)
- `A_ISO` (Real): on-site energy (k-independent)
- `B_k`: anomalous coupling between (A,B) and (C,D) sublattice pairs
- `C_k`: anomalous coupling between (A,C) and (B,D) sublattice pairs
- `D_k`: normal hopping between parallel sublattice pairs (A↔D, C↔B)
- `Δ_k`: altermagnetic off-diagonal coupling
"""
function structure_factors(k::AbstractVector, p::HematiteParams)
    S, J₁, J₂, J₃, J₄, J₅, J₁₃, Δ = p.S, p.J₁, p.J₂, p.J₃, p.J₄, p.J₅, p.J₁₃, p.Δ

    # On-site energy (k-independent part of diagonal, Eq. 7a)
    # Number of neighbors: 1st(1), 2nd(3), 3rd(3), 4th(6), 5th(1), 13th(6)
    A_ISO = S * (J₁ - 3*J₂ + 3*J₃ + 6*J₄ - J₅ - 6*J₁₃)

    # D_k: normal hopping A↔D and C↔B (Eq. 7b)
    # sum over 3 second-neighbor bonds + δ₅ + 6 thirteenth-neighbor bonds
    D_k = S * (
          J₂  * sum(exp(im * dot(k, δ₂[j]))  for j in 1:3)
        + J₅  * exp(im * dot(k, δ₅))
        + J₁₃ * sum(exp(im * dot(k, δ₁₃[j])) for j in 1:6)
    )

    # B_k: anomalous pairing between (A,-B) and (C,-D) (Eq. 7c)
    B_k = S * (
          J₁ * exp(im * dot(k, δ₁))
        + J₃ * sum(exp(im * dot(k, δ₃[j])) for j in 1:3)
    )

    # C_k: anomalous pairing between (A,-C) and (B,-D) (Eq. 7d)
    C_k = 2*S * J₄ * sum(cos(dot(k, δ₄[j])) for j in 1:3)

    # Δ_k: altermagnetic off-diagonal term (Eq. 7e)
    Δ_k = S * Δ * (
          sum(exp(im * dot(k, δ₁₃[j])) for j in 1:3)   # type-1 (red)
        - sum(exp(im * dot(k, δ₁₃[j])) for j in 4:6)   # type-2 (orange)
    )

    return A_ISO, B_k, C_k, D_k, Δ_k
end

# ---------------------------------------------------------------------------
# 8×8 BdG Hamiltonian kernel  H_k  (Eq. 8 of Hoyer et al.)
# ---------------------------------------------------------------------------

"""
    hamiltonian_nonrel(k, p) -> H (8×8 Complex matrix)

Construct the nonrelativistic BdG kernel H_k in the Nambu basis
Ψ†_k = (a†, b†, c†, d†, a_{-k}, b_{-k}, c_{-k}, d_{-k}).

H_k = H_ISO(k) + H_AM(k)   [Eq. 8 of Hoyer et al.]

The matrix is Hermitian. For physical energies diagonalize with colpa_energies().
"""
function hamiltonian_nonrel(k::AbstractVector, p::HematiteParams)
    A, B, C, D, Δ = structure_factors(k, p)

    # H_ISO block (Eq. 8 of Hoyer et al., first matrix)
    # Basis: (a†, b†, c†, d†, a_{-k}, b_{-k}, c_{-k}, d_{-k})
    # C_k is real (sum of cosines), so conj(C) = C everywhere.
    H = ComplexF64[
        A        0        0        D       0        B        C        0;
        0        A        conj(D)  0       conj(B)  0        0        C;
        0        D        A        0       C        0        0        B;
        conj(D)  0        0        A       0        C        conj(B)  0;
        0        B        C        0       A        0        0        D;
        conj(B)  0        0        C       0        A        conj(D)  0;
        C        0        0        B       0        D        A        0;
        0        C        conj(B)  0       conj(D)  0        0        A
    ]

    # H_AM block (Eq. 8, second matrix)
    H[1,4] += Δ;       H[4,1] += conj(Δ)
    H[2,3] += conj(Δ); H[3,2] += Δ
    H[5,8] += Δ;       H[8,5] += conj(Δ)
    H[6,7] += conj(Δ); H[7,6] += Δ

    return H
end

# ---------------------------------------------------------------------------
# Magnon energies via Colpa diagonalization
# ---------------------------------------------------------------------------

"""
    magnon_energies(k, p) -> energies (length-4 vector, meV)

Return the four physical magnon energies at crystal momentum k,
sorted in ascending order.
"""
function magnon_energies(k::AbstractVector, p::HematiteParams=default_params())
    H = hamiltonian_nonrel(k, p)
    return colpa_energies(H)
end

# ---------------------------------------------------------------------------
# Spin polarization  ⟨s^N̂⟩
#
# In the nonrelativistic collinear case H_k is block-diagonal when
# re-expressed in the basis:
#   Block 1:  (a†, d†, b_{-k}, c_{-k})  → spin quantum number  +1/2
#   Block 2:  (b†, c†, a_{-k}, d_{-k})  → spin quantum number  −1/2
#
# After Colpa diagonalization of the full 8×8 matrix, we identify which
# eigenvector belongs to which block by computing the overlap with block
# projectors.  In practice, for the nonrelativistic Hamiltonian we can
# simply diagonalize the two 4×4 blocks separately and assign spin.
# ---------------------------------------------------------------------------

# Block structure of H_k (both ISO and AM terms are block-diagonal):
# Block 1: rows/cols {1,4,6,7} = (a†, d†, b_{-k}, c_{-k}) → spin +½
# Block 2: rows/cols {2,3,5,8} = (b†, c†, a_{-k}, d_{-k}) → spin −½
const _IDX1 = [1, 4, 6, 7]
const _IDX2 = [2, 3, 5, 8]

"""
    magnon_energies_spin(k, p) -> (energies, spin_vals)

Return all 4 magnon energies and their spin expectation values ⟨s^N̂⟩.

In the nonrelativistic collinear case H_k is block-diagonal in the bases:
  Block 1:  (a†, d†, b_{-k}, c_{-k})  → spin +½
  Block 2:  (b†, c†, a_{-k}, d_{-k})  → spin −½

# Returns
- `energies`: 4-vector (meV); first 2 from block 1 (spin +½), last 2 from block 2 (spin −½)
- `spin_vals`: 4-vector of ⟨s^N̂⟩ values (+0.5 or −0.5)
"""
function magnon_energies_spin(k::AbstractVector, p::HematiteParams=default_params())
    H = hamiltonian_nonrel(k, p)

    H1 = H[_IDX1, _IDX1]
    H2 = H[_IDX2, _IDX2]

    e1 = colpa_energies(Hermitian(H1))   # 2 energies, spin +½
    e2 = colpa_energies(Hermitian(H2))   # 2 energies, spin −½

    energies  = vcat(e1, e2)
    spin_vals = vcat(fill(+0.5, length(e1)), fill(-0.5, length(e2)))

    return energies, spin_vals
end

# ---------------------------------------------------------------------------
# Band structure along a k-path
# ---------------------------------------------------------------------------

"""
    compute_bands_nonrel(kpts, p) -> (energies, spin_vals)

Compute magnon band structure along a k-path.

# Arguments
- `kpts`: (N × 3) matrix of k-vectors (Cartesian, Å⁻¹)
- `p`: HematiteParams

# Returns
- `energies`: (N × 4) matrix of energies (meV)
- `spin_vals`: (N × 4) matrix of ⟨s^N̂⟩ values
"""
function compute_bands_nonrel(kpts::AbstractMatrix, p::HematiteParams=default_params())
    N = size(kpts, 1)
    energies  = zeros(Float64, N, 4)
    spin_vals = zeros(Float64, N, 4)

    for i in 1:N
        k = kpts[i, :]
        e, s = magnon_energies_spin(k, p)
        energies[i, :]  = e
        spin_vals[i, :] = s
    end

    return energies, spin_vals
end

# ---------------------------------------------------------------------------
# Altermagnetic splitting of acoustic bands
# Δε_ac = ε₃(k) − ε₄(k)   [Eq. 9 of Hoyer et al.]
# (band 3 = spin +½ acoustic top, band 4 = spin −½ acoustic = lowest)
# ---------------------------------------------------------------------------

"""
    acoustic_splitting(kpts, p) -> Δε_ac (N-vector, meV)

Compute the altermagnetic splitting of the two acoustic magnon bands
along the k-path `kpts`.
"""
function acoustic_splitting(kpts::AbstractMatrix, p::HematiteParams=default_params())
    N = size(kpts, 1)
    splitting = zeros(Float64, N)
    for i in 1:N
        k = kpts[i, :]
        e, s = magnon_energies_spin(k, p)
        # acoustic bands are the two lowest-energy bands (indices 1 and 2
        # in the sorted e1 block, and indices 1 and 2 in the e2 block)
        # The two acoustic bands have one from each spin block
        e_ac_plus  = e[1]   # lowest spin +½
        e_ac_minus = e[3]   # lowest spin −½
        splitting[i] = abs(e_ac_plus - e_ac_minus)
    end
    return splitting
end
