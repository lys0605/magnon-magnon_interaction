# fig_hematite_bands.jl
#
# Reproduce Fig. 5(a) of Hoyer et al. (2025):
#   Magnon dispersion of hematite along Γ–Z–P–B–Γ–B'–F–Γ–X–X'
#   colored by spin expectation value ⟨s^N̂⟩.
#
# Requires: run compute_hematite_bands.jl first.
#
# Run from the magnon/ directory:
#   include("figures/fig_hematite_bands.jl")

include(joinpath(@__DIR__, "..", "src", "plot_utils.jl"))

using JLD2, CairoMakie, LaTeXStrings

set_paper_theme!(fontsize=12, width=620, height=420)

# --- load data ---
@load joinpath(@__DIR__, "..", "data", "hematite_bands_nonrel.jld2") dists tick_dists tick_labels energies spin_vals

# Build LaTeXString tick labels with standard notation
# (the stored labels are plain strings, convert them)
tex_labels = map(tick_labels) do l
    d = Dict("Γ"=>"Γ", "Z"=>"Z", "P"=>"P", "B"=>"B", "B'"=>"B'",
             "F"=>"F", "X"=>"X", "X'"=>"X'")
    get(d, l, l)
end

# --- figure ---
fig = Figure(size=(680, 440))
ax  = Axis(fig[1, 1];
           ylabel         = "Energy (meV)",
           xticks         = (tick_dists, tex_labels),
           xticklabelsize = 14,
           yticklabelsize = 12,
           ylabelsize     = 13,
           yminorticksvisible = true,
           yminorticks    = IntervalsBetween(5),
           limits         = (nothing, nothing, 0, nothing))

# Draw each band, colored by ⟨s^N̂⟩
# Use discrete coloring: +0.5 → magenta, −0.5 → cyan, 0 → gray
nbands = size(energies, 2)
cmap   = cgrad([:cyan, :gray80, :magenta], [0.0, 0.5, 1.0])

# Draw spin-down first so spin-up (magenta) appears on top at split points
for b in [3, 4, 1, 2]
    scatter!(ax, dists, energies[:, b];
             color      = spin_vals[:, b],
             colormap   = cmap,
             colorrange = (-0.5, 0.5),
             markersize = 2.5,
             strokewidth = 0)
end

# High-symmetry vertical lines
for td in tick_dists
    vlines!(ax, td; color=:black, linewidth=0.8, linestyle=:solid)
end

# Colorbar
cb = Colorbar(fig[1, 2];
              colormap   = cmap,
              limits     = (-0.5, 0.5),
              label      = L"\langle s^{\hat{N}} \rangle",
              labelsize  = 13,
              ticks      = [-0.5, 0.0, 0.5],
              width      = 14)

# Label: degenerate regions shown in gray
Label(fig[2, 1:2],
      "Degenerate bands shown in gray. Magenta: ⟨s⟩ = +1/2, Cyan: ⟨s⟩ = −1/2";
      fontsize=10, halign=:left, tellwidth=false)

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)
save(joinpath(outdir, "fig_hematite_bands_nonrel.png"), fig; px_per_unit=3)
println("Saved: figures/output/fig_hematite_bands_nonrel.png")

display(fig)
