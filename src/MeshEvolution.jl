export coin_operator, beam_splitter_operator
export shift_timebins, shift_timebins_operator
export shift_timebins_sp, shift_timebins_operator_sp, coin_operator_sp
export mesh_evolution
export mesh_evolution_sp
export phase_on_density_matrix
export coin_operator_sp_identical, coin_operator_identical
export shift_timebins_sp_identical, shift_timebins_operator_sp_identical
export mesh_evolution_sp_identical
export mesh_evolution_identical


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

function mesh_evolution(initial_state, angles)
    state = convert(Matrix{ComplexF64}, initial_state)::Matrix{ComplexF64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = _iterative_mesh_evolution(state, angles)
    return state
end

"""
    _iterative_mesh_evolution(input_state, angles)

Iteratively apply the coin- and bin-shifting operators to the two-photon `input state`
object.

The `input state` argument can be both a wave function in form of a Vector or SparseVector
or a density matrix (in form of a Matrix or Sparse Matrix). The type of the return object
`state` matches the type of `input_state`. The internal numerics are type-specialized too.

See also [`mesh_evolution`](@ref), [`iterative_mesh_evolution_sp`](@ref).
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

function coin_operator_sp_identical(angles::Vector)
    N = length(angles)
    d_hilbert_space = N * (2 * N + 1)
    real_angles = convert(Vector{Float64}, angles)::Vector{Float64}
    coin_op = spzeros(ComplexF64, d_hilbert_space, d_hilbert_space)

    for l in 0:N - 1
        for m in l:N - 1
            cc = cos(real_angles[l + 1]) * cos(real_angles[m + 1])
            ss = sin(real_angles[l + 1]) * sin(real_angles[m + 1])
            cs = sqrt(2) * im * cos(real_angles[l + 1]) * sin(real_angles[m + 1])

            j_ss = lcmk2j_identical(N, l, 0, m, 0)
            j_sl = lcmk2j_identical(N, l, 0, m, 1)
            j_ll = lcmk2j_identical(N, l, 1, m, 1)

            if l == m
                beam_splitter_op = [[cc, cs, -ss] [cs, cc - ss, cs] [-ss, cs, cc]]

                j_arr = [j_ss, j_sl, j_ll]
            else
                j_ls = lcmk2j_identical(N, l, 1, m, 0)

                beam_splitter_op_l = beam_splitter_operator(real_angles[l + 1])
                beam_splitter_op_m = beam_splitter_operator(real_angles[m + 1])
                beam_splitter_op = kron(beam_splitter_op_l, beam_splitter_op_m)


                j_arr = [j_ss, j_sl, j_ls, j_ll]

            end
            coin_op[j_arr, j_arr] = beam_splitter_op
        end
    end

    return coin_op::SparseMatrixCSC{ComplexF64, Int64}
end

function coin_operator_identical(angles::Vector)
    coin_operator_single_party = coin_operator_sp_identical(angles)
    tensor_coin_operator = kron(coin_operator_single_party, coin_operator_single_party)
    return tensor_coin_operator::SparseMatrixCSC{ComplexF64, Int64}
end

function shift_timebins_sp_identical end

function shift_timebins_sp_identical(state_vec::Vector)
    d_hilbert_space = length(state_vec)
    N = Int(-1 / 4 + sqrt(1 / 16 + d_hilbert_space / 2)) # p-q formular
    state_vec = convert(Vector{ComplexF64}, state_vec)::Vector{ComplexF64}
    new_vec = zeros(ComplexF64, (N + 1) * (2 * (N + 1) + 1))
    for j in 1:d_hilbert_space
        _shift_j_sp_identical!(N, j, state_vec, new_vec)
    end

    return new_vec
end

function shift_timebins_sp_identical(state_vec::SparseVector)
    d_hilbert_space = length(state_vec)
    N = Int(-1 / 4 + sqrt(1 / 16 + d_hilbert_space / 2)) # p-q formular
    state_vec =
        convert(SparseVector{ComplexF64, Int64}, state_vec)::SparseVector{ComplexF64, Int64}
    new_vec = spzeros(ComplexF64, (N + 1) * (2 * (N + 1) + 1))
    for j in state_vec.nzind
        _shift_j_sp_identical!(N, j, state_vec, new_vec)
   end

    return new_vec
end

function _shift_j_sp_identical!(N, j, state_vec, new_vec)
    l, c , m , k  = j2lcmk_identical(N, j)
    l_shift = l + c
    m_shift = m + k
    if l_shift  == m_shift && c > k
        j_shift = lcmk2j_identical(N + 1, l_shift, k, m_shift, c)
    else
        j_shift = lcmk2j_identical(N + 1, l_shift, c, m_shift, k)
    end

    # adapted system has one more time bin, so we need to put N + 1
    new_vec[j_shift] = state_vec[j]
    return nothing
end

function shift_timebins end

function shift_timebins_identical(state_vec::Vector)
    d_hilbert_space = Int(sqrt(length(state_vec)))
    N = Int(-1 / 4 + sqrt(1 / 16 + d_hilbert_space / 2)) # p-q formular
    state_vec = convert(Vector{ComplexF64}, state_vec)::Vector{ComplexF64}
    d_hilbert_space_shifted = ((N + 1) * (2 * (N + 1) + 1)) # square for two species
    new_vec = zeros(ComplexF64, d_hilbert_space_shifted^2)
    for j_super in eachindex(state_vec)
        _shift_j_identical!(N, j_super, state_vec, new_vec)
   end

    return new_vec
end

function shift_timebins_identical(state_vec::SparseVector)
    d_hilbert_space = Int(sqrt(length(state_vec)))
    N = Int(-1 / 4 + sqrt(1 / 16 + d_hilbert_space / 2)) # p-q formular
    state_vec =
        convert(SparseVector{ComplexF64, Int64}, state_vec)::SparseVector{ComplexF64, Int64}
    d_hilbert_space_shifted = ((N + 1) * (2 * (N + 1) + 1))
    new_vec = spzeros(ComplexF64, d_hilbert_space_shifted^2)
    for j_super in state_vec.nzind
        _shift_j_identical!(N, j_super, state_vec, new_vec)
    end

    return new_vec
end

function _shift_j_identical!(N, j_super, state_vec, new_vec)
    l1, c1, m1, k1, l2, c2, m2, k2 = j_super2lcmk_identical(N, j_super)
    l_shift1 = l1 + c1
    m_shift1 = m1 + k1
    if l_shift1 == m_shift1 && c1 > k1
        j_shift1 = lcmk2j_identical(N + 1, l_shift1, k1, m_shift1, c1)
    else
        j_shift1 = lcmk2j_identical(N + 1, l_shift1, c1, m_shift1, k1)
    end
    l_shift2 = l2 + c2
    m_shift2 = m2 + k2
    if l_shift2 == m_shift2 && c2 > k2
        j_shift2 = lcmk2j_identical(N + 1, l_shift2, k2, m_shift2, c2)
    else
        j_shift2 = lcmk2j_identical(N + 1, l_shift2, c2, m_shift2, k2)
    end
    d_hilbert_space_shift = (N + 1) * (2 * (N + 1) + 1)
    j_shift = lm2j(d_hilbert_space_shift, j_shift1 - 1, j_shift2 - 1)
    # adapted system has one more time bin, so we need to put N + 1
    new_vec[j_shift] = state_vec[j_super]

    return nothing
end

function shift_timebins_operator_sp_identical(N)
    N = convert(Int64, N)::Int64
    d_hilbert_space = N * (2 * N + 1)
    d_hilbert_space_shifted = (N + 1) * (2 * (N + 1) + 1)

    shift_op_single_party =
        spzeros(Int64, d_hilbert_space_shifted, d_hilbert_space_shifted)

	for j in 1:d_hilbert_space
        l, c, m, k = j2lcmk_identical(N, j)
        l_shift = l + c
        m_shift = m + k
        if l_shift  == m_shift && c > k
            j_shift = lcmk2j_identical(N + 1, l_shift, k, m_shift, c)
        else
            j_shift = lcmk2j_identical(N + 1, l_shift, c, m_shift, k)
        end
            shift_op_single_party[j_shift, j] = 1
        end

    return shift_op_single_party::SparseMatrixCSC{Float64, Int64}
end

function shift_timebins_operator_identical(N)
    shift_operator_single_party = shift_timebins_operator_sp_identical(N)
    tensor_shift_operator = kron(shift_operator_single_party, shift_operator_single_party)
    return tensor_shift_operator::SparseMatrixCSC{Float64, Int64}
end



function mesh_evolution_sp_identical end

function mesh_evolution_sp_identical(initial_state::Vector, angles)
    state = convert(Vector{ComplexF64}, initial_state)::Vector{ComplexF64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = _iterative_mesh_evolution_sp_identical(state, angles)
    return state
end

function mesh_evolution_sp_identical(initial_state::SparseVector, angles)
    state =
    convert(SparseVector{ComplexF64, Int64}, initial_state)::SparseVector{ComplexF64, Int64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = _iterative_mesh_evolution_sp_identical(state, angles)
    return state
end

function mesh_evolution_sp_identical(initial_state::AbstractMatrix, angles)
    state = convert(Matrix{ComplexF64}, initial_state)::Matrix{ComplexF64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = _iterative_mesh_evolution_sp_identical(state, angles)
    return state
end




function _iterative_mesh_evolution_sp_identical end

function _iterative_mesh_evolution_sp_identical(input_state::AbstractVector, angles)
    state = copy(input_state)
    for i in eachindex(angles)
        coin_op = coin_operator_sp_identical(angles[i])
        state = coin_op * state # apply beam splitters
        state = shift_timebins_sp_identical(state) # shift time bins accordingly
   end

    return state
end

function _iterative_mesh_evolution_sp_identical(state::AbstractMatrix, angles)
    for i in eachindex(angles)
        coin_op = coin_operator_sp_identical(angles[i])
        shift_op = shift_timebins_operator_sp_identical(length(angles[i]))

        state = coin_op * state * coin_op' # apply beam splitters
        state = shift_op * state * shift_op' # apply time-bin shift operator
   end

    return state
end

function mesh_evolution_identical end

function mesh_evolution_identical(initial_state::Vector, angles)
    state = convert(Vector{ComplexF64}, initial_state)::Vector{ComplexF64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = _iterative_mesh_evolution_identical(state, angles)
    return state::Vector{ComplexF64}
end

function mesh_evolution_identical(initial_state::SparseVector, angles)
    state =
    convert(SparseVector{ComplexF64, Int64}, initial_state)::SparseVector{ComplexF64, Int64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = _iterative_mesh_evolution_identical(state, angles)
    return state::SparseVector{ComplexF64, Int64}
end

function mesh_evolution_identical(initial_state, angles)
    state = convert(Matrix{ComplexF64}, initial_state)::Matrix{ComplexF64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = _iterative_mesh_evolution_identical(state, angles)
    return state
end


function _iterative_mesh_evolution end

function _iterative_mesh_evolution_identical(input_state::AbstractVector, angles)
    state = copy(input_state)
    for i in eachindex(angles)
        coin_op = coin_operator_identical(angles[i])
        state = coin_op * state # apply beam splitters
        state = shift_timebins_identical(state) # shift time bins accordingly
   end

    return state
end

function _iterative_mesh_evolution_identical(input_state::AbstractMatrix, angles)
    state = copy(input_state)
    for i in eachindex(angles)
        coin_op = coin_operator_identical(angles[i])
        shift_op = shift_timebins_operator_identical(length(angles[i]))
        state = coin_op * state * coin_op' # apply beam splitters
        state = shift_op * state * shift_op' # apply time-bin shift operator
   end

    return state
end
