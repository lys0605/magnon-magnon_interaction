# test_sanity.jl
#
# Quick sanity checks — run from magnon/ directory before compute scripts.
#
#   include("test_sanity.jl")

include("src/math_utils.jl")
include("src/hematite_lattice.jl")
include("src/hematite_hamiltonian.jl")

using LinearAlgebra

println("=== Hematite sanity checks ===\n")

# 1. Bond distances
println("Bond distances (Å):")
println("  1st  : ", round(norm(δ₁),   digits=4), "  (nearest A-B)")
println("  2nd  : ", round(norm(δ₂₁),  digits=4))
println("  3rd  : ", round(norm(δ₃₁),  digits=4))
println("  4th  : ", round(norm(δ₄₁),  digits=4))
println("  5th  : ", round(norm(δ₅),   digits=4))
println("  13th : ", round(norm(δ₁₃₁), digits=4), "  (paper: 6.42278 Å)")

# All 13th-neighbor bonds should have the same length
d13 = [norm(δ₁₃[i]) for i in 1:6]
println("\nAll 6 13th-neighbor distances: ", round.(d13, digits=4))
@assert all(abs.(d13 .- d13[1]) .< 1e-6) "13th-neighbor bonds not all equal length!"
println("  ✓ All equal")

# 2. Hamiltonian is Hermitian at a test k-point
println("\nHamiltonian Hermiticity check:")
p = default_params()
k_test = [0.1, 0.2, 0.05]
H = hamiltonian_nonrel(k_test, p)
err = maximum(abs.(H - H'))
println("  max |H - H†| = $(round(err, sigdigits=3))")
@assert err < 1e-12 "Hamiltonian is not Hermitian!"
println("  ✓ Hermitian")

# 3. Block-diagonal structure: off-block elements should be zero
println("\nBlock-diagonal structure check:")
off_block = maximum(abs.(H[_IDX1, _IDX2]))
println("  max |H[block1, block2]| = $(round(off_block, sigdigits=3))")
@assert off_block < 1e-12 "H is not block-diagonal!"
println("  ✓ Block-diagonal")

# 4. Energies at Γ: acoustic modes should be ~0, optical ~finite
println("\nEnergies at Γ:")
e_G, s_G = magnon_energies_spin(Γ_k, p)
println("  Block 1 (spin +½): ", round.(e_G[1:2], digits=4), " meV")
println("  Block 2 (spin -½): ", round.(e_G[3:4], digits=4), " meV")
@assert e_G[1] < 1.0 "Acoustic mode at Γ should be ~0 meV"
println("  ✓ Acoustic modes near zero at Γ")

# 5. Energies at B point
# B = (1/2,0,0) is a TRIM: -k ≡ k (mod b₁), so antiunitary symmetry forces
# both spin sectors to be degenerate there.  Splitting must be zero by symmetry.
println("\nEnergies at B (TRIM — splitting must vanish by symmetry):")
e_B, s_B = magnon_energies_spin(B_k, p)
println("  Block 1 (spin +½): ", round.(e_B[1:2], digits=4), " meV")
println("  Block 2 (spin -½): ", round.(e_B[3:4], digits=4), " meV")
Δε_B = abs(e_B[1] - e_B[3])
println("  Acoustic splitting at B: $(round(Δε_B, sigdigits=3)) meV (expect ~0 at TRIM)")
@assert Δε_B < 1e-6 "Splitting at TRIM B should be zero by symmetry"
println("  ✓ Zero splitting at TRIM B")

# 5b. Generic off-TRIM k-point: splitting should be finite
println("\nEnergies at generic k = (0.3, 0.1, 0.05) (should show splitting):")
k_gen = _to_cart([0.3, 0.1, 0.05])
e_gen, _ = magnon_energies_spin(k_gen, p)
Δε_gen = abs(e_gen[1] - e_gen[3])
println("  Acoustic splitting: $(round(Δε_gen, digits=4)) meV")
@assert Δε_gen > 0.001 "Expected finite altermagnetic splitting at generic k"
println("  ✓ Finite splitting at generic off-TRIM k")

# 6. Energies at Z (nodal path Γ-Z: should be degenerate)
println("\nEnergies at Z (Γ-Z is a nodal path, should be degenerate):")
e_Z, s_Z = magnon_energies_spin(Z_k, p)
println("  Block 1 (spin +½): ", round.(e_Z[1:2], digits=4), " meV")
println("  Block 2 (spin -½): ", round.(e_Z[3:4], digits=4), " meV")
Δε_Z = abs(e_Z[1] - e_Z[3])
println("  Acoustic splitting at Z: $(round(Δε_Z, digits=6)) meV (should be ~0)")

println("\n=== All checks passed ===")
