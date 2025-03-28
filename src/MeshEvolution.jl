export coin_operator, beam_splitter_operator
export shift_timebins, shift_timebins_operator
export shift_timebins_sp, shift_timebins_operator_sp, coin_operator_sp
export mesh_evolution, mesh_evolution_backwards
export mesh_evolution_sp
export phase_on_density_matrix

"""
    shift_timebins_sp(state_vec::Vector)
    shift_timebins_sp(state_vec::SparseVector)

Shift the time bins of a single photon pure state according to the current loop index.

Accepts both a dense or sparse `state_vec` argument. The return type of the new_vec matches
the type of `state_vec`.

See also [`shift_timebins`](@ref).
"""
function shift_timebins_sp end

function shift_timebins_sp(state_vec::Vector)
    state_vec = convert(Vector{ComplexF64}, state_vec)::Vector{ComplexF64}
    new_vec = Vector{ComplexF64}(undef, length(state_vec)+N_LOOPS)
    new_vec[2] = 0.0
    new_vec[end - 1] = 0.0
    new_vec[1:2:end - 3] = @view state_vec[1:2:end]
    new_vec[4:2:end] = @view state_vec[2:2:end]
    return new_vec
end

function shift_timebins_sp(state_vec::SparseVector)
    state_vec =
        convert(SparseVector{ComplexF64, Int64}, state_vec)::SparseVector{ComplexF64, Int64}
    new_vec = spzeros(ComplexF64, length(state_vec)+N_LOOPS)
    for j in state_vec.nzind
        _shift_j_sp!(j, state_vec, new_vec)
   end

    return new_vec
end

"""
    shift_timebins(state_vec::Vector)
    shift_timebins(state_vec::SparseVector)

Shift the time bins of a two-photon pure state according to the loop index.

Accepts both a dense or sparse `state_vec` argument. The return type of the new_vec matches
the type of `state_vec`.

See also [`shift_timebins_sp`](@ref), [`_shift_j!`](@ref).
"""
function shift_timebins end

function shift_timebins(state_vec::Vector)
    state_vec = convert(Vector{ComplexF64}, state_vec)::Vector{ComplexF64}
    N = Int64(sqrt(length(state_vec) / (N_LOOPS2)))
    new_vec = zeros(ComplexF64, ((N + 1) * N_LOOPS)^2)
    for j in eachindex(state_vec)
        _shift_j!(N, j, state_vec, new_vec)
   end

    return new_vec
end

function shift_timebins(state_vec::SparseVector)
    state_vec =
        convert(SparseVector{ComplexF64, Int64}, state_vec)::SparseVector{ComplexF64, Int64}
    N = Int64(sqrt(length(state_vec) / (N_LOOPS2)))
    new_vec = spzeros(ComplexF64, ((N + 1) * N_LOOPS)^2)
    for j in state_vec.nzind
        _shift_j!(N, j, state_vec, new_vec)
   end

    return new_vec
end

"""
    _shift_j!(N, j, state_vec, new_vec)

Shift the two-photon coefficient of index `j` from the old `state_vec` with `N` time bins to
the new, longer `new_vec` state vector with `N + 1` timebins.

See also [`shift_timebins`](@ref), [`_shift_j_sp!`](@ref).
"""
function _shift_j!(N, j, state_vec, new_vec)
    l, c , m , k  = j2lcmk(N, j)
    shifted_j = lcmk2j(N + 1, l + c, c, m + k, k)
    # adapted system has one more time bin, so we need to put N + 1
    new_vec[shifted_j] = state_vec[j]
    return nothing
end

"""
    _shift_j_sp!(j, state_vec, new_vec)

Shift the single-photon coefficient of index `j` from the old `state_vec` to the new, longer
`new_vec` state vector.

See also [`shift_timebins`](@ref), [`_shift_j!`](@ref).
"""
function _shift_j_sp!(j, state_vec, new_vec)
    l, c  = j2lc(j)
    shifted_j = lc2j(l + c, c)
    new_vec[shifted_j] = state_vec[j]
    return nothing
end

"""
    shift_timebins_operator_sp(N)

Compute the matrix operator which shifts the time bins of a single-photon state with `N`
time bins in the `|lc⟩` basis.

See also [`shift_timebins_operator`](@ref), [`coin_operator_sp`](@ref),
[`mesh_evolution_sp`](@ref).
"""
function shift_timebins_operator_sp(N)
    N = convert(Int64, N)::Int64

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

See also [`shift_timebins_operator_sp`](@ref), [`coin_operator`](@ref),
[`mesh_evolution`](@ref).
"""
function shift_timebins_operator(N)
    shift_operator_single_photon = shift_timebins_operator_sp(N)
    tensor_shift_operator = kron(shift_operator_single_photon, shift_operator_single_photon)
    return tensor_shift_operator::SparseMatrixCSC{Float64, Int64}
end



"""
    beam_splitter_operator(θ)

Return the sparse 2x2 single-photon beam-splitter unitary parametrized by
splitting angle `θ`.

See also [`coin_operator`](@ref), [`coin_operator_sp`](@ref).
"""
function beam_splitter_operator(θ)
    θ = convert(Float64,θ)::Float64
    cs = cos(θ)
    sn = im * sin(θ)
    #cols = [1, 1, 2, 2]
    #rows = [1, 2, 1, 2]
    vals = [cs, sn, sn, cs]
   return sparse(_cols, _rows, vals)
end

"""
    coin_operator_sp(angles::Vector)

Compute the complete sparse single-photon beam-splitter unitary for all time bins within one
round trip.

The `angles` argument contains the values of the parameter `θ` for each time bin of that
round trip.

See also [`beam_splitter_operator`](@ref), [`coin_operator`](@ref).
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

See also [`beam_splitter_operator`](@ref), [`coin_operator_sp`](@ref).
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

See also [`explicit_state_evolution`](@ref), [`mesh_evolution_sp`](@ref).
"""
function mesh_evolution end

function mesh_evolution(initial_state::Vector, angles)
    state = convert(Vector{ComplexF64}, initial_state)::Vector{ComplexF64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = _iterative_mesh_evolution(state, angles)
    return state::Vector{ComplexF64}
end

function mesh_evolution(initial_state::SparseVector, angles)
    state =
    convert(SparseVector{ComplexF64, Int64}, initial_state)::SparseVector{ComplexF64, Int64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = _iterative_mesh_evolution(state, angles)
    return state::SparseVector{ComplexF64, Int64}
end

function mesh_evolution(initial_state::AbstractMatrix, angles)
    state = convert(Matrix{ComplexF64}, initial_state)::Matrix{ComplexF64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = _iterative_mesh_evolution(state, angles)
    return state
end

"""
    mesh_evolution_backwards(final_state::Vector, angles)
    mesh_evolution_backwards(final_state::SparseVector, angles)
    mesh_evolution_backwards(final_state, angles)

Numerically compute the backwards unitary evolution of the two-photon final  `final_state`
through the fiber mesh network parametrized through `angles`, i.e. the inverse of the mesh
evolution given by `mesh_evolution`.

Accepts both a dense or sparse Vector or a matrix `final_state` argument. The internal
numerics and the return type of the `state` object after the backwards evolution match the
type of `final_state` in the case of a Vector.
For matrices, only a dense version is currently implemented.

# Arguments
- `final_state`: the two-photon final state after the mesh network. Can be in form of a
    Vector or SparseVector for a pure state or in form of a Matrix for a density matrix.
    Must be in the |lcmk⟩ basis.
- `angles`: a vector or vectors of beam-splitter angles. The number of round trips matches
    `length(angles)`.

See also [`explicit_state_evolution`](@ref), [`mesh_evolution_sp`](@ref).
"""
function mesh_evolution_backwards end

function mesh_evolution_backwards(final_state::Vector, angles)
    state = convert(Vector{ComplexF64}, final_state)::Vector{ComplexF64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = _iterative_mesh_evolution_backwards(state, angles)
    N = length(angles[1])
    j_arr_forbidden = j_arr_backwards_forbidden(N)
    state[j_arr_forbidden] .= 0
    return state::Vector{ComplexF64}
end

function mesh_evolution_backwards(final_state::SparseVector, angles)
    state =
    convert(SparseVector{ComplexF64, Int64}, final_state)::SparseVector{ComplexF64, Int64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = _iterative_mesh_evolution_backwards(state, angles)
    N = length(angles[1])
    j_arr_forbidden = j_arr_backwards_forbidden(N)
    state[j_arr_forbidden] .= 0
    return state::SparseVector{ComplexF64, Int64}
end

function mesh_evolution_backwards(final_state::AbstractMatrix, angles)
    state = convert(Matrix{ComplexF64}, final_state)::Matrix{ComplexF64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = _iterative_mesh_evolution_backwards(state, angles)
    N = length(angles[1])
    j_arr_forbidden = j_arr_backwards_forbidden(N)
    state[j_arr_forbidden, :] .= 0
    state[:, j_arr_forbidden] .= 0
    return state
end

"""
    j_arr_backwards_forbidden(N)

Generate a list of forbidden inintial state `j` indices for the backwards evolution of the
two-photon state. The forbidden indices are those where at least one photon is in the long
loop.

"""
function j_arr_backwards_forbidden(N)
    j_arr_forbidden = []
    for l in 0:N - 1
        for m in 0:N - 1
            j_sl = lcmk2j(N, l, 0, m, 1)
            j_ls = lcmk2j(N, l, 1, m, 0)
            j_ll = lcmk2j(N, l, 1, m, 1)
            push!(j_arr_forbidden, j_sl)
            push!(j_arr_forbidden, j_ls)
            push!(j_arr_forbidden, j_ll)
        end
    end

    return j_arr_forbidden
end

"""
    _iterative_mesh_evolution(input_state, angles)

Iteratively apply the coin- and bin-shifting operators to the two-photon `input state`
object.

The `input state` argument can be both a wave function in form of a Vector or SparseVector
or a density matrix (in form of a Matrix or Sparse Matrix). The type of the return object
`state` matches the type of `input_state`. The internal numerics are type-specialized too.

See also [`mesh_evolution`](@ref), [`_iterative_mesh_evolution`](@ref),
[`iterative_mesh_evolution_sp`](@ref), [`mesh_evolution_backwards`](@ref).
"""
function _iterative_mesh_evolution end

function _iterative_mesh_evolution(input_state::AbstractVector, angles)
    state = copy(input_state)
    for i in eachindex(angles)
        coin_op = coin_operator(angles[i])
        state = coin_op * state # apply beam splitters
        state = shift_timebins(state) # shift time bins accordingly
    end

    return state
end

function _iterative_mesh_evolution(input_state::AbstractMatrix, angles)
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
    _iterative_mesh_evolution_backwards(input_state, angles)

Iteratively apply the coin- and bin-shifting operators to the two-photon `input state`
object, going through the network backwards.

The `input state` argument can be both a wave function in form of a Vector or SparseVector
or a density matrix (in form of a Matrix or Sparse Matrix). The type of the return object
`state` matches the type of `input_state`. The internal numerics are type-specialized too.

See also [`mesh_evolution`](@ref), [`_iterative_mesh_evolution`](@ref),
[`iterative_mesh_evolution_sp`](@ref), [`mesh_evolution_backwards`](@ref).
"""
function _iterative_mesh_evolution_backwards end

function _iterative_mesh_evolution_backwards(input_state::AbstractVector, angles)
    state = copy(input_state)
    for i in lastindex(angles):-1:1
        shift_op = shift_timebins_operator(length(angles[i]))
        coin_op = coin_operator(angles[i])
        state = shift_op' * state # apply inverse time-bin shift operator
        state = coin_op' * state # apply inverse beam splitters

    end

    return state
end

function _iterative_mesh_evolution_backwards(input_state::AbstractMatrix, angles)
    state = copy(input_state)
    for i in lastindex(angles):-1:1
        coin_op = coin_operator(angles[i])
        shift_op = shift_timebins_operator(length(angles[i]))
        state = shift_op' * state * shift_op # apply inverse time-bin shift operator
        state = coin_op' * state * coin_op # apply inverse beam splitters

    end

    return state
end

"""
    mesh_evolution_sp(initial_state::Vector, angles)
    mesh_evolution_sp(initial_state::SparseVector, angles)
    mesh_evolution_sp(initial_state::AbstractMatrix, angles)

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

See also [`explicit_state_evolution`](@ref), [`mesh_evolution`](@ref).
"""
function mesh_evolution_sp end

function mesh_evolution_sp(initial_state::Vector, angles)
    state = convert(Vector{ComplexF64}, initial_state)::Vector{ComplexF64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = _iterative_mesh_evolution_sp(state, angles)
    return state
end

function mesh_evolution_sp(initial_state::SparseVector, angles)
    state =
    convert(SparseVector{ComplexF64, Int64}, initial_state)::SparseVector{ComplexF64, Int64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = _iterative_mesh_evolution_sp(state, angles)
    return state
end

function mesh_evolution_sp(initial_state::AbstractMatrix, angles)
    state = convert(Matrix{ComplexF64}, initial_state)::Matrix{ComplexF64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = _iterative_mesh_evolution_sp(state, angles)
    return state
end

"""
_iterative_mesh_evolution_sp(input_state, angles)

Iteratively apply the coin- and bin-shifting operators to the single-photon `input state`
Vector.

The `input state` argument can be both a wave function in form of a Vector or SparseVector
or a density matrix (in form of a Matrix or Sparse Matrix). The type of the return object
`state` matches the type of `input_state`. The internal numerics are type-specialized too.

See also  [`_iterative_mesh_evolution`](@ref), [`mesh_evolution_sp`](@ref).
"""
function _iterative_mesh_evolution_sp end

function _iterative_mesh_evolution_sp(input_state::AbstractVector, angles)
    state = copy(input_state)
    for i in eachindex(angles)
        coin_op = coin_operator_sp(angles[i])
        state = coin_op * state # apply beam splitters
        state = shift_timebins_sp(state) # shift time bins accordingly
   end

    return state
end

function _iterative_mesh_evolution_sp(state::AbstractMatrix, angles)
    for i in eachindex(angles)
        coin_op = coin_operator_sp(angles[i])
        shift_op = shift_timebins_operator_sp(length(angles[i]))

        state = coin_op * state * coin_op' # apply beam splitters
        state = shift_op * state * shift_op' # apply time-bin shift operator
   end

    return state
end
