# fig_hematite_splitting.jl
#
# Reproduce Fig. 5(c) of Hoyer et al. (2025):
#   Altermagnetic acoustic splitting Δε_ac along Γ–B.
#   Shows the characteristic k⁴ scaling near Γ.
#
# Requires: run compute_hematite_bands.jl first.
#
# Run from the magnon/ directory:
#   include("figures/fig_hematite_splitting.jl")

include(joinpath(@__DIR__, "..", "src", "plot_utils.jl"))

using JLD2, CairoMakie, LaTeXStrings

set_paper_theme!(fontsize=12, width=460, height=340)

# --- load data ---
@load joinpath(@__DIR__, "..", "data", "hematite_splitting.jld2") dists_GB splitting

# --- figure ---
fig = Figure(size=(480, 340))
ax  = Axis(fig[1, 1];
           xlabel         = L"|k|\ (\mathrm{\AA}^{-1})",
           ylabel         = L"\Delta\varepsilon_\mathrm{ac}\ \mathrm{(meV)}",
           xticklabelsize = 12,
           yticklabelsize = 12,
           xlabelsize     = 13,
           ylabelsize     = 13)

lines!(ax, dists_GB, splitting; color=:royalblue, linewidth=2.2)

# Overlay k⁴ fit near Γ (first 30% of path for fitting)
N_fit  = round(Int, 0.3 * length(dists_GB))
k_fit  = dists_GB[2:N_fit]
s_fit  = splitting[2:N_fit]
# linear fit of log(Δε) vs log(|k|) to extract exponent
log_k  = log.(k_fit)
log_s  = log.(max.(s_fit, 1e-12))
A_mat  = hcat(ones(length(log_k)), log_k)
coeffs = A_mat \ log_s
exponent = coeffs[2]
prefactor = exp(coeffs[1])

k_dense  = range(dists_GB[2], dists_GB[N_fit], length=100)
k4_curve = prefactor .* k_dense .^ exponent
lines!(ax, k_dense, k4_curve;
       color=:orange, linewidth=1.8, linestyle=:dash,
       label=L"k^{%$(round(exponent, digits=1))}\ \mathrm{fit}")

axislegend(ax; position=:rb, framevisible=false, labelsize=11)

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)
save(joinpath(outdir, "fig_hematite_splitting.png"), fig; px_per_unit=3)
println("Saved: figures/output/fig_hematite_splitting.png")
println("  Fitted exponent: $(round(exponent, digits=3))  (expected ≈ 4.0)")

display(fig)
