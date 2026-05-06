isdefined(Main, :_PLOT_UTILS_LOADED) && return
const _PLOT_UTILS_LOADED = true

using CairoMakie
using LaTeXStrings

# ---------------------------------------------------------------------------
# Paper theme
# ---------------------------------------------------------------------------

"""
    set_paper_theme!(; fontsize=11, width=500, height=380)

Apply a clean publication-quality theme to all subsequent CairoMakie figures.
"""
function set_paper_theme!(; fontsize=11, width=500, height=380)
    set_theme!(
        fontsize = fontsize,
        size     = (width, height),
        Axis = (
            xgridvisible   = false,
            ygridvisible   = false,
            spinewidth     = 1.2,
            xtickwidth     = 1.2,
            ytickwidth     = 1.2,
            xticksize      = 5.0,
            yticksize      = 5.0,
            xticklabelsize = fontsize,
            yticklabelsize = fontsize,
            xlabelsize     = fontsize + 1,
            ylabelsize     = fontsize + 1,
            titlesize      = fontsize + 1,
        ),
        Colorbar = (
            ticksize       = 4.0,
            tickwidth      = 1.0,
            ticklabelsize  = fontsize,
            labelsize      = fontsize,
        ),
        Legend = (
            framevisible   = false,
            labelsize      = fontsize,
        ),
        Lines   = (linewidth = 1.8,),
        Scatter = (markersize = 5,),
    )
end

reset_theme!() = set_theme!()

# ---------------------------------------------------------------------------
# Magnon band plot with spin-polarization coloring
# ---------------------------------------------------------------------------

"""
    plot_bands_spin(dists, energies, spin_vals;
                   tick_dists, labels,
                   ylabel="Energy (meV)",
                   clims=(-0.5, 0.5),
                   colormap=:RdBu,
                   title="") -> fig

Plot magnon band structure colored by spin expectation value.

# Arguments
- `dists`: N-vector of cumulative k-path distances
- `energies`: N × nbands matrix of energies
- `spin_vals`: N × nbands matrix of ⟨s^N̂⟩ values (−0.5 to +0.5)
- `tick_dists`: positions of high-symmetry-point ticks on x-axis
- `labels`: LaTeXString labels at each tick
- `clims`: colormap range for spin polarization
"""
function plot_bands_spin(dists, energies, spin_vals;
                         tick_dists=Float64[],
                         labels=String[],
                         ylabel="Energy (meV)",
                         clims=(-0.5, 0.5),
                         colormap=:RdBu,
                         title="")
    fig = Figure()
    ax  = Axis(fig[1, 1];
               title    = title,
               ylabel   = ylabel,
               xticks   = (tick_dists, labels),
               xticklabelsize = 13,
               yticklabelsize = 11)

    nbands = size(energies, 2)
    for b in 1:nbands
        sc = scatter!(ax, dists, energies[:, b];
                      color       = spin_vals[:, b],
                      colormap    = colormap,
                      colorrange  = clims,
                      markersize  = 3,
                      strokewidth = 0)
    end

    # vertical lines at high-symmetry points
    for td in tick_dists
        vlines!(ax, td; color=:black, linewidth=0.8)
    end

    Colorbar(fig[1, 2];
             colormap  = colormap,
             limits    = clims,
             label     = L"\langle s^{\hat{N}} \rangle",
             labelsize = 12)

    return fig
end

# ---------------------------------------------------------------------------
# Simple band plot (no spin coloring)
# ---------------------------------------------------------------------------

"""
    plot_bands(dists, energies;
              tick_dists, labels, ylabel, title, colors) -> fig

Plain line-style band plot. `colors` can be a vector of colors (one per band)
or a single color.
"""
function plot_bands(dists, energies;
                    tick_dists=Float64[],
                    labels=String[],
                    ylabel="Energy (meV)",
                    title="",
                    colors=nothing)
    fig = Figure()
    ax  = Axis(fig[1, 1];
               title  = title,
               ylabel = ylabel,
               xticks = (tick_dists, labels),
               xticklabelsize = 13,
               yticklabelsize = 11)

    nbands = size(energies, 2)
    default_colors = Makie.wong_colors()
    for b in 1:nbands
        c = colors !== nothing ? (colors isa AbstractVector ? colors[b] : colors) :
                                  default_colors[mod1(b, length(default_colors))]
        lines!(ax, dists, energies[:, b]; color=c)
    end

    for td in tick_dists
        vlines!(ax, td; color=:black, linewidth=0.8)
    end

    return fig
end

# ---------------------------------------------------------------------------
# Line plot (generic)
# ---------------------------------------------------------------------------

"""
    plot_line(xs, ys; xlabel, ylabel, title, label, yscale) -> fig
"""
function plot_line(xs, ys;
                   xlabel="",
                   ylabel="",
                   title="",
                   label=nothing,
                   yscale=identity)
    fig = Figure()
    ax  = Axis(fig[1, 1]; title=title, xlabel=xlabel, ylabel=ylabel, yscale=yscale)
    label !== nothing ? lines!(ax, xs, ys; label=label) : lines!(ax, xs, ys)
    label !== nothing && axislegend(ax)
    return fig
end

# ---------------------------------------------------------------------------
# 2D heatmap
# ---------------------------------------------------------------------------

"""
    heatmap2d(xs, ys, Z; title, xlabel, ylabel, colormap, colorrange) -> fig
"""
function heatmap2d(xs, ys, Z;
                   title="", xlabel="", ylabel="",
                   colormap=:viridis,
                   colorrange=nothing)
    fig = Figure()
    ax  = Axis(fig[1, 1]; title=title, xlabel=xlabel, ylabel=ylabel, aspect=DataAspect())
    kw  = Dict{Symbol,Any}(:colormap => colormap)
    colorrange !== nothing && (kw[:colorrange] = colorrange)
    hm = heatmap!(ax, xs, ys, Z; kw...)
    Colorbar(fig[1, 2], hm)
    return fig
end

# ---------------------------------------------------------------------------
# Save helper
# ---------------------------------------------------------------------------

"""
    save_fig(fig, path; px_per_unit=3)

Save a CairoMakie figure. Format inferred from extension (png, pdf, svg).
"""
function save_fig(fig, path; px_per_unit=3)
    save(path, fig; px_per_unit=px_per_unit)
    println("Saved: $path")
end
