# compute_hematite_bands.jl
#
# Compute nonrelativistic hematite magnon band structure and save to data/.
#
# Run from the magnon/ directory:
#   include("compute_hematite_bands.jl")
#
# Output files (saved to data/):
#   hematite_bands_nonrel.jld2   — bands + spin along standard k-path
#   hematite_splitting.jld2      — acoustic splitting along Γ–B

include("src/math_utils.jl")
include("src/hematite_lattice.jl")
include("src/hematite_hamiltonian.jl")

using JLD2

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
p = default_params()

println("Hematite nonrelativistic magnon band structure")
println("  S = $(p.S),  J₁ = $(p.J₁) meV,  J₃ = $(p.J₃) meV")
println("  J₁₃ = $(p.J₁₃) meV,  Δ = $(p.Δ) meV")

# ---------------------------------------------------------------------------
# Standard k-path (Γ–Z–P–B–Γ–B'–F–Γ–X–X')
# ---------------------------------------------------------------------------
waypoints = hematite_kpath()
kpts, dists, tick_dists, tick_labels = build_kpath(waypoints; N_total=600)

println("\nComputing bands along k-path ($(size(kpts,1)) points)...")
energies, spin_vals = compute_bands_nonrel(kpts, p)
println("  Done. Energy range: $(round(minimum(energies), digits=2)) – $(round(maximum(energies), digits=2)) meV")

mkpath("data")
@save "data/hematite_bands_nonrel.jld2" kpts dists tick_dists tick_labels energies spin_vals p
println("  Saved: data/hematite_bands_nonrel.jld2")

# ---------------------------------------------------------------------------
# Altermagnetic splitting along Γ–B
# ---------------------------------------------------------------------------
println("\nComputing acoustic splitting along Γ–B...")
splitting_path = [("Γ", Γ_k), ("X", X_k)]
kpts_GB, dists_GB, _, _ = build_kpath(splitting_path; N_total=300)
splitting = acoustic_splitting(kpts_GB, p)

@save "data/hematite_splitting.jld2" kpts_GB dists_GB splitting p
println("  Saved: data/hematite_splitting.jld2")

println("\nAll done.")
