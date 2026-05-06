isdefined(Main, :_HEMATITE_LATTICE_LOADED) && return
const _HEMATITE_LATTICE_LOADED = true

using LinearAlgebra

# ---------------------------------------------------------------------------
# Hematite (α-Fe₂O₃) crystal structure
#
# Space group:  R3̄c  (No. 167)
# Bravais lattice: rhombohedral (trigonal)
# Reference: Hoyer et al., arXiv:2503.11623 (2025)
#            Pailhé et al., J. Solid State Chem. 181, 2933 (2008)
#
# Conventions
# -----------
# All lengths in Ångström (Å).
# Lattice vectors a₁, a₂, a₃ are columns of the matrix A.
# Iron positions r_A,...,r_D are given in fractional coordinates of the
# primitive (rhombohedral) unit cell.
# Bond vectors δ are in Cartesian coordinates (Å).
# Reciprocal lattice vectors b₁, b₂, b₃ satisfy aᵢ · bⱼ = 2π δᵢⱼ.
# ---------------------------------------------------------------------------

# --- Conventional hexagonal lattice constants (Å) ---
const a_c = 5.0342   # in-plane
const c_c = 13.7519  # out-of-plane

# --- Conventional hexagonal lattice vectors (as columns) ---
# a₁ᶜ, a₂ᶜ, a₃ᶜ  (Eq. 1 of Hoyer et al.)
const A_conv = [
    sqrt(2 + sqrt(3))/2 * a_c    -sqrt(2 - sqrt(3))/2 * a_c    0.0;
   -sqrt(2 - sqrt(3))/2 * a_c     sqrt(2 + sqrt(3))/2 * a_c    0.0;
    0.0                            0.0                           c_c
]

# --- Transformation to primitive rhombohedral vectors (Eq. 2 of Hoyer et al.) ---
# (a₁ a₂ a₃) = (a₁ᶜ a₂ᶜ a₃ᶜ) M
const M_conv2prim = [
   -1/3   -1/3    2/3;
    1/3   -2/3    1/3;
    1/3    1/3    1/3
]

# Primitive lattice vectors as columns: A[:,i] = aᵢ  (Å, Cartesian)
const A = A_conv * M_conv2prim   # 3×3 matrix, columns are a₁, a₂, a₃

const a₁ = A[:, 1]
const a₂ = A[:, 2]
const a₃ = A[:, 3]

# --- Reciprocal lattice vectors (Å⁻¹) ---
# bᵢ satisfy aᵢ · bⱼ = 2π δᵢⱼ
const B = 2π * inv(A)'   # columns are b₁, b₂, b₃

const b₁ = B[:, 1]
const b₂ = B[:, 2]
const b₃ = B[:, 3]

# ---------------------------------------------------------------------------
# Iron atom positions in fractional coordinates of the primitive cell
# δ_Fe = 0.35498 − 1/3  (Pailhé et al. 2008, Table 1)
# ---------------------------------------------------------------------------
const δ_Fe = 0.35498 - 1/3

const r_A_frac = [1/6 - δ_Fe,  1/6 - δ_Fe,  1/6 - δ_Fe]
const r_B_frac = [1/3 + δ_Fe,  1/3 + δ_Fe,  1/3 + δ_Fe]
const r_C_frac = [2/3 - δ_Fe,  2/3 - δ_Fe,  2/3 - δ_Fe]
const r_D_frac = [5/6 + δ_Fe,  5/6 + δ_Fe,  5/6 + δ_Fe]

# Cartesian positions (Å)
frac_to_cart(f) = A * f

const r_A = frac_to_cart(r_A_frac)
const r_B = frac_to_cart(r_B_frac)
const r_C = frac_to_cart(r_C_frac)
const r_D = frac_to_cart(r_D_frac)

# ---------------------------------------------------------------------------
# C₃ rotation matrix (120° about c-axis, i.e. the [111] direction of rhomb.)
# The c-axis in Cartesian is ẑ = (0, 0, 1).
# For the rhombohedral system, the 3-fold axis is along the body diagonal.
# We use the explicit 3×3 matrix from the space-group symmetry.
# ---------------------------------------------------------------------------
const θ_C3 = 2π / 3
const C3 = [
    cos(θ_C3)   -sin(θ_C3)   0.0;
    sin(θ_C3)    cos(θ_C3)   0.0;
    0.0          0.0          1.0
]
const C3_inv = C3'   # C3⁻¹ = C3†

# ---------------------------------------------------------------------------
# Bond vectors (Cartesian, Å)
# Following the notation of Hoyer et al., Eqs. (5a)–(5f)
# ---------------------------------------------------------------------------

# 1st neighbor: A → B
const δ₁ = r_B - r_A

# 2nd neighbor: A → D  (three bonds related by C₃)
const δ₂₁ = (r_D - r_A) - a₂ - a₃
const δ₂₂ = C3  * δ₂₁
const δ₂₃ = C3_inv * δ₂₁
const δ₂ = [δ₂₁, δ₂₂, δ₂₃]

# 3rd neighbor: A → B  (three bonds related by C₃)
const δ₃₁ = (r_B - r_A) - a₃
const δ₃₂ = C3  * δ₃₁
const δ₃₃ = C3_inv * δ₃₁
const δ₃ = [δ₃₁, δ₃₂, δ₃₃]

# 4th neighbor: A → C  (three bonds related by C₃)
const δ₄₁ = (r_C - r_A) - a₂ - a₃
const δ₄₂ = C3  * δ₄₁
const δ₄₃ = C3_inv * δ₄₁
const δ₄ = [δ₄₁, δ₄₂, δ₄₃]

# 5th neighbor: B → C  (one bond along c)
const δ₅ = r_B - r_C

# ---------------------------------------------------------------------------
# Mirror matrices (Householder reflections in Cartesian space)
# M_n̂ = I − 2 n̂ n̂ᵀ  where n̂ = (aᵢ − aⱼ)/|aᵢ − aⱼ|
# M_{-110}: perpendicular to (a₁ − a₂)
# M_{-101}: perpendicular to (a₁ − a₃)
# M_{0-11}: perpendicular to (a₂ − a₃)
# ---------------------------------------------------------------------------
function _mirror_mat(n::AbstractVector)
    n̂ = n / norm(n)
    return I - 2 * n̂ * n̂'
end

const M_m110 = _mirror_mat(a₁ - a₂)   # M_{-110}
const M_m101 = _mirror_mat(a₁ - a₃)   # M_{-101}
const M_0m11 = _mirror_mat(a₂ - a₃)   # M_{0-11}

# 13th neighbor: A → D  (six bonds, two inequivalent sets of 3)
# Type 1 "red"   bonds: J₁₃ + Δ  — three related by C₃
# Type 2 "orange" bonds: J₁₃ − Δ  — three obtained by mirrors of δ₁₃₁
#
# Note: the paper's Eq.(5f) writes δ₁₃₃ = C₃⁻¹δ₂₁, which is a typo;
# the physically correct expression is C₃⁻¹δ₁₃₁ (same 13th-shell distance).
const δ₁₃₁ = (r_D - r_A) - a₂ - 2*a₃
const δ₁₃₂ = C3     * δ₁₃₁
const δ₁₃₃ = C3_inv * δ₁₃₁
const δ₁₃₄ = M_m110 * δ₁₃₁
const δ₁₃₅ = M_m101 * δ₁₃₁
const δ₁₃₆ = M_0m11 * δ₁₃₁
# Collect: first 3 = "type-1" (red), last 3 = "type-2" (orange)
const δ₁₃ = [δ₁₃₁, δ₁₃₂, δ₁₃₃, δ₁₃₄, δ₁₃₅, δ₁₃₆]   # 6 bonds total

# ---------------------------------------------------------------------------
# Verify bond distances (for sanity)
# ---------------------------------------------------------------------------
function print_bond_distances()
    println("Bond distances (Å):")
    println("  1st: ", round(norm(δ₁),    digits=4))
    println("  2nd: ", round(norm(δ₂₁),   digits=4))
    println("  3rd: ", round(norm(δ₃₁),   digits=4))
    println("  4th: ", round(norm(δ₄₁),   digits=4))
    println("  5th: ", round(norm(δ₅),    digits=4))
    println("  13th:", round(norm(δ₁₃₁),  digits=4), " (paper: 6.42278 Å)")
end

# ---------------------------------------------------------------------------
# High-symmetry points in the rhombohedral BZ (Cartesian, Å⁻¹)
# Reference: Hoyer et al. Fig. 5(b); standard R3̄c BZ conventions.
#
# The rhombohedral BZ has the shape of a trigonally-distorted cube.
# Key points in fractional reciprocal coordinates (q₁,q₂,q₃) of b₁,b₂,b₃:
#   Γ  = (0, 0, 0)
#   Z  = (1/2, 1/2, 1/2)        — zone boundary along the 3-fold axis
#   F  = (1/2, 1/2, 0)          — face center
#   B  = (1/2, 0, 0)            — (or symmetry-equivalent point)
#   L  = (1/2, 0, 1/2)
#   X  = (1/4, 1/2, -1/4)       — approximate; adjust as needed
#   P  = (3/8, 3/4, 3/8)
# ---------------------------------------------------------------------------
_to_cart(q) = B * q   # fractional → Cartesian reciprocal

const Γ_k  = _to_cart([0.0,   0.0,   0.0])
const Z_k  = _to_cart([0.5,   0.5,   0.5])
const F_k  = _to_cart([0.5,   0.5,   0.0])
const L_k  = _to_cart([0.5,   0.0,   0.5])
const B_k  = _to_cart([0.5,   0.0,   0.0])
const Bp_k = _to_cart([0.0,   0.5,   0.0])   # B' = C₃-related partner of B
const P_k  = _to_cart([0.375, 0.375, 0.75])   # between Z and F
const X_k  = _to_cart([0.5,   0.0,  -0.5])
const Xp_k = _to_cart([0.5,   0.5,  -0.5])   # X'

"""
    hematite_kpath() -> waypoints

Return the standard k-path for hematite used in Hoyer et al. Fig. 5(a):
Γ – Z – P – B – Γ – B' – F – Γ – X – X'
"""
function hematite_kpath()
    return [
        ("Γ",  Γ_k),
        ("Z",  Z_k),
        ("P",  P_k),
        ("B",  B_k),
        ("Γ",  Γ_k),
        ("B'", Bp_k),
        ("F",  F_k),
        ("Γ",  Γ_k),
        ("X",  X_k),
        ("X'", Xp_k),
    ]
end
