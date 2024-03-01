export coin_operator
export shift_timebins, beam_splitter_operator, mesh_evolution
export shift_timebins_sp, coin_operator_sp, mesh_evolution_sp
export phase_on_density_matrix

"Numerical cutoff value to determine wether a certain coefficient is ≈ 0"
global const weight_cutoff = 1e-16

"""
    shift_timebins_sp(state_vec::Vector)

Shift the time bins of a single photon state according to the current loop index.

The state after operation is returned as a `new_vec` Vector{ComplexF64} object.

See also shift_timebins.
"""
function shift_timebins_sp(state_vec::Vector)
    state_vec = convert(Vector{ComplexF64}, state_vec)::Vector{ComplexF64}
    new_vec = Vector{ComplexF64}(undef, length(state_vec)+n_loops)
    new_vec[2] = 0
    new_vec[end-1] = 0
    new_vec[1:2:end-3] = @view state_vec[1:2:end]
    new_vec[4:2:end] = @view state_vec[2:2:end]
    return new_vec
end

"""
    shift_timebins(state_vec::Vector)
    shift_timebins(state_vec::SparseVector)

Shift the time bins of a two-photon state according to the loop index.
Accepts both a dense or sparse `state_vec` argument. The return type of the new_vec matches the type of `state_vec`.

See also shift_timebins_sp, shift_j!.
"""
function shift_timebins end

function shift_timebins(state_vec::Vector)
    state_vec = convert(Vector{ComplexF64}, state_vec)::Vector{ComplexF64}
    N = Int64(sqrt(length(state_vec)/(n_loops2)))
    new_vec = zeros(ComplexF64, ((N+1)*n_loops)^2)
    for j in eachindex(state_vec)
        shift_j!(N, j, state_vec, new_vec)
    end
    return new_vec
end

function shift_timebins(state_vec::SparseVector)
    state_vec = convert(SparseVector{ComplexF64, Int64}, state_vec)::SparseVector{ComplexF64, Int64}
    N = Int64(sqrt(length(state_vec)/(n_loops2)))
    new_vec = spzeros(ComplexF64, ((N+1)*n_loops)^2)
    for j in state_vec.nzind
        shift_j!(N, j, state_vec, new_vec)
    end
    return new_vec
end

function shift_j!(N, j, state_vec, new_vec)
    l,c,m,k = j2lcmk(N,j)
    shifted_j = lcmk2j(N+1, l+c, c, m+k, k) # adapted system has one more time bin, so we need to put N+1
    new_vec[shifted_j] = state_vec[j]
end

"""
    beam_splitter_operator(θ)

Return the sparse 2x2 single-photon beam-splitter unitary parametrized by splitting angle `θ`.

See also coin_operator.
"""
function beam_splitter_operator(θ)
    θ = convert(Float64,θ)::Float64
    cs = cos(θ)
    sn = im*sin(θ)
    cols = [1,1,2,2]
    rows = [1,2,1,2]
    vals = [cs, sn, sn, cs]
   return sparse(cols,rows,vals)
end

"""
    coin_operator_sp(angles::Vector)

Compute the complete sparse single-photon beam-splitter unitary for all time bins within one round trip.

The `angles` argument contains the values of the parameter `θ` for each time bin of that round trip.

See also `beam_splitter_operator`, `coin_operator`.
"""
function coin_operator_sp(angles::Vector)
    real_angles = convert(Vector{Float64}, angles)::Vector{Float64}
    matrices = [beam_splitter_operator(θ) for θ in real_angles]
    coin_operator_single_photon = blockdiag(matrices...)
    return coin_operator_single_photon::SparseMatrixCSC{ComplexF64, Int64}
end


"""
    coin_operator(angles::Vector)

Compute the complete sparse two-photon beam-splitter unitary for all time bins within one round trip.

The `angles` argument contains the values of the parameter `θ` for each time bin of that round trip.

See also `beam_splitter_operator`, `coin_operator_sp`.
"""
function coin_operator(angles::Vector)
    coin_operator_single_photon = coin_operator_sp(angles)
    tensor_coin_operator = kron(coin_operator_single_photon,coin_operator_single_photon)
    return tensor_coin_operator::SparseMatrixCSC{ComplexF64, Int64}
end

"""
    mesh_evolution(ψ_init::Vector, angles)
    mesh_evolution(ψ_init::SparseVector, angles)

Numerically compute the unitary evolution of the initial state `Ψ_init` through the fiber mesh network parametrized through `angles`.

Accepts both a dense or sparse `ψ_init` argument. The internal numerics and the return type of the `state` object match the type of `ψ_init`.

# Arguments
- `ψ_init`: the pure intial state before the mesh network. Must be in the |lcmk> basis.
- `angles`: a vector or vectors of beam-splitter angles. The number of round trips matches `length(angles)`.

See also `explicit_state_evolution`.
"""
function mesh_evolution end

function mesh_evolution(ψ_init::Vector, angles)
    state = convert(Vector{ComplexF64}, ψ_init)::Vector{ComplexF64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = iterative_mesh_evolution(state, angles)
    return state
end

function mesh_evolution(ψ_init::SparseVector, angles)
    state = convert(SparseVector{ComplexF64, Int64}, ψ_init)::SparseVector{ComplexF64, Int64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = iterative_mesh_evolution(state, angles)
    return state
end

function iterative_mesh_evolution(state, angles)
    for i in eachindex(angles)
        coin_op = coin_operator(angles[i])
        state = coin_op * state # apply beam splitters
        state = shift_timebins(state) # shift time bins accordingly
    end
    return state
end

function mesh_evolution end

function mesh_evolution_sp(ψ_init::Vector, angles)
    state = convert(Vector{ComplexF64}, ψ_init)::Vector{ComplexF64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = iterative_mesh_evolution_sp(state, angles)
    return state
end

function mesh_evolution_sp(ψ_init::SparseVector, angles)
    state = convert(SparseVector{ComplexF64, Int64}, ψ_init)::SparseVector{ComplexF64, Int64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = iterative_mesh_evolution_sp(state, angles)
    return state
end

function iterative_mesh_evolution_sp(state, angles)
    for i in eachindex(angles)
        coin_op = coin_operator_sp(angles[i])
        state = coin_op * state # apply beam splitters
        state = shift_timebins_sp(state) # shift time bins accordingly
    end
    return state
end


"""
    phase_on_density_matrix(ρ, φ_arr)

Apply phases `φ_arr` to the correlated time bins of the density matrix `ρ`

Returns a new density matrix `ρ` after phase application. Each time bin |ii⟩, i ∈ {0,1,…} is subjected to phase `φ_arr[i+1]`.

See also `initial_state_phase_estimation`.
"""
function phase_on_density_matrix(ρ, φ_arr)
    ρ_rot = convert(Matrix{ComplexF64}, copy(ρ))::Matrix{ComplexF64}
    φ_arr = convert(Vector{Float64}, φ_arr)::Vector{Float64}
    N = Int64(sqrt(size(ρ)[1]/(n_loops2)))

    @argcheck length(φ_arr) == N
    
    tb_idxs = [lcmk2j(N,i,0,i,0) for i in 0:N-1]
    for (idx1, j1) in enumerate(tb_idxs), (idx2, j2) in enumerate(tb_idxs)
        ρ_rot[j1,j2] *= cis(φ_arr[idx1]-φ_arr[idx2])
    end
    return ρ_rot
end