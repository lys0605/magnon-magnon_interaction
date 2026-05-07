# fig_hematite_splitting.jl
#
# Altermagnetic acoustic splitting Δε_ac of hematite magnons.
#
# Panel 1 (left):  splitting along the full band path Γ–Z–P–B–Γ–B'–F–Γ–X–X'
#   showing the characteristic pattern: zero at all TRIM points, peaked between them.
# Panel 2 (right): log-log plot along a non-nodal direction [1,0.3,0] (frac. rec.)
#   showing the g-wave k⁴ power law near Γ.
#
# Requires: run compute_hematite_bands.jl first.

include(joinpath(@__DIR__, "..", "src", "plot_utils.jl"))
include(joinpath(@__DIR__, "..", "src", "math_utils.jl"))
include(joinpath(@__DIR__, "..", "src", "hematite_lattice.jl"))
include(joinpath(@__DIR__, "..", "src", "hematite_hamiltonian.jl"))

using JLD2, CairoMakie, LaTeXStrings, LinearAlgebra

set_paper_theme!(fontsize=12, width=860, height=350)

# --- load full band path data ---
@load joinpath(@__DIR__, "..", "data", "hematite_bands_nonrel.jld2") dists tick_dists tick_labels energies spin_vals

# Acoustic splitting from full path: col 1 = lowest spin+½, col 3 = lowest spin-½
splitting_full = abs.(energies[:, 1] .- energies[:, 3])

tex_labels = map(tick_labels) do l
    d = Dict("Γ"=>"Γ", "Z"=>"Z", "P"=>"P", "B"=>"B", "B'"=>"B'",
             "F"=>"F", "X"=>"X", "X'"=>"X'")
    get(d, l, l)
end

# --- panel 2: k^4 check along non-nodal direction [1,0.3,0] (frac. reciprocal) ---
p_params = default_params()
dir_frac  = [1.0, 0.3, 0.0]
dir_cart  = B * dir_frac; dir_cart /= norm(dir_cart)
ks_dense  = 10 .^ range(log10(0.003), log10(0.4), length=80)
split_nd  = Float64[]
for k in ks_dense
    kv = k * dir_cart
    H  = hamiltonian_nonrel(kv, p_params)
    e1 = colpa_energies(Hermitian(H[[1,4,6,7],[1,4,6,7]]))
    e2 = colpa_energies(Hermitian(H[[2,3,5,8],[2,3,5,8]]))
    push!(split_nd, abs(e1[1] - e2[1]))
end

# Fit log-log slope in range k = 0.005–0.05
mask = 0.005 .< ks_dense .< 0.05
lk   = log.(ks_dense[mask])
ls   = log.(split_nd[mask])
c    = hcat(ones(sum(mask)), lk) \ ls
exponent  = c[2]
prefactor = exp(c[1])

# --- figure: two panels side by side ---
fig = Figure(size=(900, 350))

# Panel 1: full-path splitting
ax1 = Axis(fig[1,1];
           ylabel         = L"\Delta\varepsilon_\mathrm{ac}\ \mathrm{(meV)}",
           xticks         = (tick_dists, tex_labels),
           xticklabelsize = 13,
           yticklabelsize = 11,
           ylabelsize     = 13,
           title          = "Splitting along Γ–Z–P–B–Γ–B'–F–Γ–X–X'",
           titlesize      = 11)

lines!(ax1, dists, splitting_full; color=:royalblue, linewidth=2.0)
for td in tick_dists
    vlines!(ax1, td; color=:black, linewidth=0.7, linestyle=:solid)
end
ylims!(ax1, 0, nothing)

# Panel 2: log-log k^4 verification
ax2 = Axis(fig[1,2];
           xlabel         = L"\log_{10}|k|\ (\mathrm{\AA}^{-1})",
           ylabel         = L"\log_{10}\,\Delta\varepsilon_\mathrm{ac}\ \mathrm{(meV)}",
           xticklabelsize = 11,
           yticklabelsize = 11,
           xlabelsize     = 12,
           ylabelsize     = 12,
           title          = "k⁴ scaling (non-nodal direction [1,0.3,0])",
           titlesize      = 11)

scatter!(ax2, log10.(ks_dense), log10.(max.(split_nd, 1e-20));
         color=:royalblue, markersize=4)

# Overlay fitted slope line
k_fit_range = [0.003, 0.08]
y_fit = @. log10(prefactor) + exponent * log10(k_fit_range)
lines!(ax2, log10.(k_fit_range), y_fit;
       color=:orange, linewidth=2.0, linestyle=:dash,
       label=L"k^{%$(round(exponent, digits=1))}\ \mathrm{fit}")

axislegend(ax2; position=:rb, framevisible=false, labelsize=11)

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)
save(joinpath(outdir, "fig_hematite_splitting.png"), fig; px_per_unit=3)
println("Saved: figures/output/fig_hematite_splitting.png")
println("  Fitted exponent (non-nodal direction): $(round(exponent, digits=3))  (expected 4.0)")
