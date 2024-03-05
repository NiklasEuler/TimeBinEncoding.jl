export coin_operator, beam_splitter_operator, mesh_evolution
export shift_timebins, shift_timebins_operator
export shift_timebins_sp, shift_timebins_operator_sp, coin_operator_sp, mesh_evolution_sp
export phase_on_density_matrix

const global N_LOOPS = 2
# number of fiber loops. Saved as const to avoid magic numbers.
const global N_LOOPS2 = 4
# number of fiber loops squared. Saved as const to avoid magic numbers.
global const WEIGHT_CUTOFF = 1e-16
#Numerical cutoff value to determine wether a certain coefficient is ≈ 0

"""
    shift_timebins_sp(state_vec::Vector)
    shift_timebins_sp(state_vec::SparseVector)

Shift the time bins of a single photon pure state according to the current loop index.

Accepts both a dense or sparse `state_vec` argument. The return type of the new_vec matches
the type of `state_vec`.

See also shift_timebins.
"""
function shift_timebins_sp end

function shift_timebins_sp(state_vec::Vector)
    state_vec = convert(Vector{ComplexF64}, state_vec)::Vector{ComplexF64}
    new_vec = Vector{ComplexF64}(undef, length(state_vec)+N_LOOPS)
    new_vec[2] = 0
    new_vec[end - 1] = 0
    new_vec[1:2:end - 3] = @view state_vec[1:2:end]
    new_vec[4:2:end] = @view state_vec[2:2:end]
    return new_vec
end

function shift_timebins_sp(state_vec::SparseVector)
    state_vec =
        convert(SparseVector{ComplexF64, Int64}, state_vec)::SparseVector{ComplexF64, Int64}
    new_vec = spzeros(ComplexF64, length(state_vec)+N_LOOPS)
    for j in state_vec.nzind
        shift_j_sp!(j, state_vec, new_vec)
   end

    return new_vec
end

"""
    shift_timebins(state_vec::Vector)
    shift_timebins(state_vec::SparseVector)

Shift the time bins of a two-photon pure state according to the loop index.

Accepts both a dense or sparse `state_vec` argument. The return type of the new_vec matches
the type of `state_vec`.

See also shift_timebins_sp, shift_j!.
"""
function shift_timebins end

function shift_timebins(state_vec::Vector)
    state_vec = convert(Vector{ComplexF64}, state_vec)::Vector{ComplexF64}
    N = Int64(sqrt(length(state_vec) / (N_LOOPS2)))
    new_vec = zeros(ComplexF64, ((N + 1) * N_LOOPS)^2)
    for j in eachindex(state_vec)
        shift_j!(N, j, state_vec, new_vec)
   end

    return new_vec
end

function shift_timebins(state_vec::SparseVector)
    state_vec =
        convert(SparseVector{ComplexF64, Int64}, state_vec)::SparseVector{ComplexF64, Int64}
    N = Int64(sqrt(length(state_vec) / (N_LOOPS2)))
    new_vec = spzeros(ComplexF64, ((N + 1) * N_LOOPS)^2)
    for j in state_vec.nzind
        shift_j!(N, j, state_vec, new_vec)
   end

    return new_vec
end

"""
    shift_j!(N, j, state_vec, new_vec)

Shift the two-photon coefficient of index `j` from the old `state_vec` with `N` time bins to
the new, longer `new_vec` state vector with `N + 1` timebins.

See also `shift_timebins`, `shift_j_sp!`.
"""
function shift_j!(N, j, state_vec, new_vec)
    l, c , m , k  = j2lcmk(N, j)
    shifted_j = lcmk2j(N + 1, l + c, c, m + k, k)
    # adapted system has one more time bin, so we need to put N + 1
    new_vec[shifted_j] = state_vec[j]
    return nothing
end

"""
    shift_j_sp!(j, state_vec, new_vec)

Shift the single-photon coefficient of index `j` from the old `state_vec` to the new, longer
`new_vec` state vector.

See also `shift_timebins`, `shift_j!`
"""
function shift_j_sp!(j, state_vec, new_vec)
    l, c  = j2lc(j)
    shifted_j = lc2j(l+c, c)
    new_vec[shifted_j] = state_vec[j]
    return nothing
end

"""
    shift_timebins_operator_sp(N)

Compute the matrix operator which shifts the time bins of a single-photon state with `N`
time bins in the `|lc⟩` basis.

See also `shift_timebins_operator`, `coin_operator_sp`, `mesh_evolution_sp`.
"""
function shift_timebins_operator_sp(N)
    N = convert(Int64, N)::Int64
    #s_diag = append!([[1, 0] for _ in 1:N]...)
    #l_diag = append!([[0, 1] for _ in 1:N]...)
    s_diag = spzeros(2 * N)
	s_diag[1:2:end - 1] .= 1
    # short bins remain unchanged, i.e. diagonal part of shift_operator_single_photon
    l_diag = spzeros(2 * N)
    l_diag[2:2:end] .= 1
    # long bins are shifted |n L⟩ → |N + 1 S⟩, so offdiagonal
    # in shift_operator_single_photon
    shift_operator_single_photon =
        dropzeros!(spdiagm(N_LOOPS * (N + 1), N_LOOPS * N, 0 => s_diag, -2 => l_diag))
    return shift_operator_single_photon::SparseMatrixCSC{Float64, Int64}
end

"""
    shift_timebins_operator(N)

Compute the matrix operator which shifts the time bins of a two-photon state with `N` time
bins in the `|lcmk⟩` basis.

See also `shift_timebins_operator_sp`, `coin_operator`, `mesh_evolution`.
"""
function shift_timebins_operator(N)
    shift_operator_single_photon = shift_timebins_operator_sp(N)
    tensor_coin_operator = kron(shift_operator_single_photon, shift_operator_single_photon)
    return tensor_coin_operator::SparseMatrixCSC{Float64, Int64}
end



"""
    beam_splitter_operator(θ)

Return the sparse 2x2 single-photon beam-splitter unitary parametrized by
splitting angle `θ`.

See also `coin_operator`, `coin_operator_sp`.
"""
function beam_splitter_operator(θ)
    θ = convert(Float64,θ)::Float64
    cs = cos(θ)
    sn = im * sin(θ)
    cols = [1, 1, 2, 2]
    rows = [1, 2, 1, 2]
    vals = [cs, sn, sn, cs]
   return sparse(cols, rows, vals)
end

"""
    coin_operator_sp(angles::Vector)

Compute the complete sparse single-photon beam-splitter unitary for all time bins within one
round trip.

The `angles` argument contains the values of the parameter `θ` for each time bin of that
round trip.

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

Compute the complete sparse two-photon beam-splitter unitary for all time bins within one
round trip.

The `angles` argument contains the values of the parameter `θ` for each time bin of that
round trip.

See also `beam_splitter_operator`, `coin_operator_sp`.
"""
function coin_operator(angles::Vector)
    coin_operator_single_photon = coin_operator_sp(angles)
    tensor_coin_operator = kron(coin_operator_single_photon, coin_operator_single_photon)
    return tensor_coin_operator::SparseMatrixCSC{ComplexF64, Int64}
end

"""
    mesh_evolution(initial_state::Vector, angles)
    mesh_evolution(initial_state::SparseVector, angles)
    mesh_evolution(initial_state, angles)

Numerically compute the unitary evolution of the two-photon initial state `initial_state`
through the fiber mesh network parametrized through `angles`.

Accepts both a dense or sparse Vector or a matrix `initial_state` argument. The internal
numerics and the return type of the `state` object after the evolution match the type of
`initial_state` in the case of a Vector.
For matrices, only a dense version is currently implemented.

# Arguments
- `initial_state`: the two-photon intial state before the mesh network. Can be in form of a
    Vector or SparseVector for a pure state or in form of a Matrix for a density matrix.
    Must be in the |lcmk⟩ basis.
- `angles`: a vector or vectors of beam-splitter angles. The number of round trips matches
    `length(angles)`.

See also `explicit_state_evolution`, `mesh_evolution_sp`.
"""
function mesh_evolution end

function mesh_evolution(initial_state::Vector, angles)
    state = convert(Vector{ComplexF64}, initial_state)::Vector{ComplexF64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = iterative_mesh_evolution(state, angles)
    return state::Vector{ComplexF64}
end

function mesh_evolution(initial_state::SparseVector, angles)
    state =
    convert(SparseVector{ComplexF64, Int64}, initial_state)::SparseVector{ComplexF64, Int64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = iterative_mesh_evolution(state, angles)
    return state::SparseVector{ComplexF64, Int64}
end

function mesh_evolution(initial_state, angles)
    state = convert(Matrix{ComplexF64}, initial_state)::Matrix{ComplexF64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = iterative_mesh_evolution_ρ(state, angles)
    return state
end

"""
    iterative_mesh_evolution(input_state, angles)

Iteratively apply the coin- and bin-shifting operators to the two-photon `input state`
Vector.

The return object `state` is either a SparseVector or a Vector type, depending on the type
of `input_state`. The internal numerics are matching that state.

See also `iterative_mesh_evolution_density_matrix`, `iterative_mesh_evolution_sp`,
`iterative_mesh_evolution_density_matrix_sp`.
"""
function iterative_mesh_evolution(input_state, angles)
    state = copy(input_state)
    for i in eachindex(angles)
        coin_op = coin_operator(angles[i])
        state = coin_op * state # apply beam splitters
        state = shift_timebins(state) # shift time bins accordingly
   end

    return state
end

"""
    iterative_mesh_evolution_ρ(input_state, angles)

Iteratively apply the coin- and bin-shifting operators to the two-photon `input state`
density matrix.

The return object `state` is a Matrix object, independent of the type of `input_state`.

See also `iterative_mesh_evolution`, `iterative_mesh_evolution_sp`,
`iterative_mesh_evolution_density_matrix_sp`.
"""
function iterative_mesh_evolution_ρ(input_state, angles)
    state = copy(input_state)
    for i in eachindex(angles)
        coin_op = coin_operator(angles[i])
        shift_op = shift_timebins_operator(length(angles[i]))
        state = coin_op * state * coin_op' # apply beam splitters
        state = shift_op * state * shift_op' # apply time-bin shift operator
   end

    return state
end

"""
    mesh_evolution_sp(initial_state::Vector, angles)
    mesh_evolution_sp(initial_state::SparseVector, angles)

Numerically compute the unitary evolution of the single-photon initial state `initial_state`
through the fiber mesh network parametrized through `angles`.

Accepts both a dense or sparse Vector or a matrix `initial_state` argument. The internal
numerics and the return type of the `state` object after the evolution match the type of
`initial_state` in the case of a Vector.
For matrices, only a dense version is currently implemented.

# Arguments
- `initial_state`: the pure single-photon intial state before the mesh network. Can be in
    form of a Vector or SparseVector for a pure state or in form of a Matrix for a density
    matrix. Must be in the |lc⟩ basis.
- `angles`: a vector or vectors of beam-splitter angles. The number of round trips matches
    `length(angles)`.

See also `explicit_state_evolution`, `mesh_evolution`.
"""
function mesh_evolution_sp end

function mesh_evolution_sp(initial_state::Vector, angles)
    state = convert(Vector{ComplexF64}, initial_state)::Vector{ComplexF64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = iterative_mesh_evolution_sp(state, angles)
    return state
end

function mesh_evolution_sp(initial_state::SparseVector, angles)
    state =
    convert(SparseVector{ComplexF64, Int64}, initial_state)::SparseVector{ComplexF64, Int64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = iterative_mesh_evolution_sp(state, angles)
    return state
end

function mesh_evolution_sp(initial_state, angles)
    state = convert(Matrix{ComplexF64}, initial_state)::Matrix{ComplexF64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = iterative_mesh_evolution_ρ_sp(state, angles)
    return state
end

"""
iterative_mesh_evolution_sp(input_state, angles)

Iteratively apply the coin- and bin-shifting operators to the single-photon `input state`
Vector.

The return object `state` is either a SparseVector or a Vector type, depending on the type
of `input_state`. The internal numerics are matching that state.

See also  `iterative_mesh_evolution_density_matrix_sp`, `iterative_mesh_evolution`,
`iterative_mesh_evolution_density_matrix`.
"""
function iterative_mesh_evolution_sp(input_state, angles)
    state = copy(input_state)
    for i in eachindex(angles)
        coin_op = coin_operator_sp(angles[i])
        state = coin_op * state # apply beam splitters
        state = shift_timebins_sp(state) # shift time bins accordingly
   end

    return state
end

"""
    iterative_mesh_evolution_ρ_sp(input_state, angles)

Iteratively apply the coin- and bin-shifting operators to the single-photon `input state`
density matrix.

The return object `state` is a Matrix object, independent of the type of `input_state`.

See also `iterative_mesh_evolution_sp`, `iterative_mesh_evolution`,
`iterative_mesh_evolution_density_matrix`.
"""
function iterative_mesh_evolution_ρ_sp(state, angles)
    for i in eachindex(angles)
        coin_op = coin_operator_sp(angles[i])
        shift_op = shift_timebins_operator_sp(length(angles[i]))
        state = coin_op * state * coin_op' # apply beam splitters
        state = shift_op * state * shift_op' # apply time-bin shift operator
   end

    return state
end


"""
    phase_on_density_matrix(ρ, φ_arr)

Apply phases `φ_arr` to the correlated time bins of the density matrix `ρ`

Returns a new density matrix `ρ` after phase application. Each time bin |ii⟩, i ∈ {0, 1,…}
is subjected to phase `φ_arr[i + 1]`.

See also `initial_state_phase_estimation`.
"""
function phase_on_density_matrix(ρ, φ_arr)
    ρ_rot = convert(Matrix{ComplexF64}, copy(ρ))::Matrix{ComplexF64}
    φ_arr = convert(Vector{Float64}, φ_arr)::Vector{Float64}
    N = Int64(sqrt(size(ρ)[1] / (N_LOOPS2)))

    @argcheck length(φ_arr) == N

    tb_idxs = [lcmk2j(N, i, 0, i, 0) for i in 0:N - 1]
    for (idx1, j1) in enumerate(tb_idxs), (idx2, j2) in enumerate(tb_idxs)
        ρ_rot[j1, j2] *= cis(φ_arr[idx1] - φ_arr[idx2])
   end

    return ρ_rot
end
