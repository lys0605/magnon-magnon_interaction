isdefined(Main, :_MATH_UTILS_LOADED) && return
const _MATH_UTILS_LOADED = true

using LinearAlgebra

# ---------------------------------------------------------------------------
# Grid helpers
# ---------------------------------------------------------------------------

"""
    meshgrid(xs, ys) -> (X, Y)

Return matrices X, Y with X[i,j] = xs[i] and Y[i,j] = ys[j].
"""
function meshgrid(xs, ys)
    X = [x for x in xs, _ in ys]
    Y = [y for _ in xs, y in ys]
    return X, Y
end

# ---------------------------------------------------------------------------
# k-path construction (works for any dimension: 2D or 3D)
# ---------------------------------------------------------------------------

"""
    build_kpath(waypoints; N_total=400) -> (kpts, dists, tick_dists, labels)

Build a k-path through high-symmetry points, distributing k-points
proportionally to segment length.

# Arguments
- `waypoints`: vector of `(label, kvec)` pairs, e.g.
               `[("Γ", [0,0,0]), ("Z", [0,0,0.5]), ...]`
- `N_total`: total number of k-points across all segments

# Returns
- `kpts`: (N × dim) matrix of k-vectors
- `dists`: N-vector of cumulative path distances (x-axis for band plots)
- `tick_dists`: positions of high-symmetry points on the x-axis
- `labels`: label strings at each high-symmetry point
"""
function build_kpath(waypoints; N_total=400)
    kvecs  = [w[2] for w in waypoints]
    labels = [w[1] for w in waypoints]
    nseg   = length(kvecs) - 1

    seg_len = [norm(kvecs[i+1] .- kvecs[i]) for i in 1:nseg]
    total   = sum(seg_len)

    seg_npts = [max(2, round(Int, N_total * l / total)) for l in seg_len]

    kpts       = Vector{Float64}[]
    dists      = Float64[]
    tick_dists = [0.0]
    cum        = 0.0

    for i in 1:nseg
        k0, k1 = kvecs[i], kvecs[i+1]
        n  = seg_npts[i]
        ts = i < nseg ? range(0, 1, length=n+1)[1:end-1] : range(0, 1, length=n)
        for t in ts
            push!(kpts,  k0 .+ t .* (k1 .- k0))
            push!(dists, cum + t * seg_len[i])
        end
        cum += seg_len[i]
        push!(tick_dists, cum)
    end

    return reduce(hcat, kpts)', dists, tick_dists, labels
end

# ---------------------------------------------------------------------------
# BdG diagonalization via direct eigendecomposition of G·H
#
# For a bosonic BdG kernel H (Hermitian), the physical magnon energies are
# the positive eigenvalues of the effective Hamiltonian M = G·H, where
# G = diag(+1,...,+1,−1,...,−1) is the bosonic metric.
# Eigenvalues of M come in ±ε pairs; only the positive ones are physical.
# This approach works at Goldstone points (no regularizer needed) and makes
# all spectral quantities directly accessible.
# ---------------------------------------------------------------------------

"""
    bdg_diagonalize(H, G) -> (energies, vecs)

Diagonalize the effective BdG Hamiltonian M = G·H directly.

# Arguments
- `H`: 2n×2n Hermitian BdG kernel
- `G`: 2n×2n bosonic metric, G = diag(I_n, -I_n)

# Returns
- `energies`: n-vector of positive magnon energies (meV), sorted ascending
- `vecs`:     2n×2n matrix whose columns are the eigenvectors of G·H,
              with positive-energy columns first (columns 1:n)
"""
function bdg_diagonalize(H::AbstractMatrix, G::AbstractMatrix)
    n2 = size(H, 1)
    n  = n2 ÷ 2

    M = G * Matrix(Hermitian(H))
    vals, vecs = eigen(M)   # eigenvalues of GH; come in ±ε pairs

    # Sort descending by real part: positive energies first
    idx  = sortperm(real.(vals), rev=true)
    vals = vals[idx]
    vecs = vecs[:, idx]

    energies = max.(real.(vals[1:n]), 0.0)   # physical (positive) energies
    return sort(energies), vecs
end

"""
    colpa_energies(H) -> energies

Return the n physical magnon energies (meV) from the 2n×2n BdG kernel H.
G = diag(I_n, -I_n) is built automatically.

Diagonalizes G·H directly; no Cholesky or regularizer needed.
"""
function colpa_energies(H::AbstractMatrix)
    n2 = size(H, 1)
    n  = n2 ÷ 2
    G  = Diagonal(vcat(ones(Float64, n), -ones(Float64, n)))
    energies, _ = bdg_diagonalize(H, G)
    return energies
end
