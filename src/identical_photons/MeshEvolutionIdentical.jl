export coin_operator_sp_identical, coin_operator_identical
export shift_timebins_sp_identical, shift_timebins_operator_sp_identical
export mesh_evolution_sp_identical
export mesh_evolution_identical

"""
    coin_operator_sp_identical(angles::AbstractVector)

Compute the coin operator for a single party of two identical photons for a given set of
beam splitter `angles`. The coin operator is a unitary matrix that describes the action of
the central beam splitter on the two-photon state during one round trip.

See also [`coin_operator_identical`](@ref), [`shift_timebins_sp_identical`](@ref),
[`shift_timebins_operator_sp_identical`](@ref), [`mesh_evolution_sp_identical`](@ref),
[`mesh_evolution_identical`](@ref).

"""
function coin_operator_sp_identical(angles::AbstractVector)
    N = length(angles)
    d_hilbert_space = N * (2 * N + 1)
        # simplified expression for N_LOOPS = 2 canceling out geometric series
    real_angles = convert(Vector{Float64}, angles)::Vector{Float64}
    coin_op = spzeros(ComplexF64, d_hilbert_space, d_hilbert_space)

    cos_vals = cos.(real_angles)
    sin_vals = sin.(real_angles)

    for l in 0:N - 1
        for m in l:N - 1
            cc = cos_vals[l + 1] * cos_vals[m + 1]
            ss = sin_vals[l + 1] * sin_vals[m + 1]
            cs = sqrt(2) * im * cos_vals[l + 1] * sin_vals[m + 1]

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
            coin_op[j_arr, j_arr] .= beam_splitter_op
        end
    end

    return coin_op::SparseMatrixCSC{ComplexF64, Int64}
end


"""
    coin_operator_identical(angles::AbstractVector)
    coin_operator_identical(angles::AbstractVector, kron_mem)

Two-party version of `coin_operator_sp_identical`. Compute the tensor product of two single-
party coin operators for two identical photons for a given set of beam splitter `angles`.

Optionally, a preallocated memory for the tensor-product operators can be provided in the
`kron_mem` argument.

See also [`coin_operator_sp_identical`](@ref), [`shift_timebins_sp_identical`](@ref),
[`shift_timebins_operator_sp_identical`](@ref), [`mesh_evolution_sp_identical`](@ref),
[`mesh_evolution_identical`](@ref).

"""
function coin_operator_identical end

function coin_operator_identical(angles::Vector)
    coin_operator_single_party = coin_operator_sp_identical(angles)
    tensor_coin_operator = kron(coin_operator_single_party, coin_operator_single_party)
    return tensor_coin_operator::SparseMatrixCSC{ComplexF64, Int64}
end

function coin_operator_identical(angles::Vector, kron_mem)
    coin_operator_single_party = coin_operator_sp_identical(angles)
    if isdefined(Main, :kron!)
        kron!(kron_mem, coin_operator_single_party, coin_operator_single_party)
    else
        kron_mem .= kron(coin_operator_single_party, coin_operator_single_party)
    end

    return kron_mem::SparseMatrixCSC{ComplexF64, Int64}
end

"""
    shift_timebins_sp_identical(state_vec::Vector)
    shift_timebins_sp_identical(state_vec::SparseVector)

Shift the time bins of a single party of two identical photons for a given state vector
based on their loop indices. The state vector is assumed to be in the `|lcmk⟩` basis.

Photon states in the long loop are delayed by one time bin, while states in the short loop
remain in their current time bin. The state vector is reshaped accordingly.

See also [`shift_timebins_identical`](@ref), [`shift_timebins_operator_sp_identical`](@ref),
[`coin_operator_sp_identical`](@ref), [`coin_operator_identical`](@ref),
[`mesh_evolution_sp_identical`](@ref), [`mesh_evolution_identical`](@ref).

"""
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
        # if l == m, we need to swap the indices to maintain the correct order
        j_shift = lcmk2j_identical(N + 1, l_shift, k, m_shift, c)
    else
        j_shift = lcmk2j_identical(N + 1, l_shift, c, m_shift, k)
    end
    # adapted system has one more time bin, so we need to put N + 1
    new_vec[j_shift] = state_vec[j]
    return nothing
end


"""
    shift_timebins_identical(state_vec::Vector)
    shift_timebins_identical(state_vec::SparseVector)

Two-party version of `shift_timebins_sp_identical`. Shift the time bins of the full two-
party quantum state for a given state vector based on their loop indices. The state vector
is assumed to be in the `|l1 c1 m1 k1 l2 c2 m2 k2⟩` basis.

Photon states in the long loop are delayed by one time bin, while states in the short loop
remain in their current time bin. The state vector is reshaped accordingly.

See also [`shift_timebins_sp_identical`](@ref), [`shift_timebins_operator_identical`](@ref),
[`coin_operator_sp_identical`](@ref), [`coin_operator_identical`](@ref),
[`mesh_evolution_sp_identical`](@ref), [`mesh_evolution_identical`](@ref).

"""
function shift_timebins_identical end

function shift_timebins_identical(state_vec::Vector)
    d_hilbert_space = Int(sqrt(length(state_vec)))
    N = Int(-1 / 4 + sqrt(1 / 16 + d_hilbert_space / 2)) # p-q formular
    state_vec = convert(Vector{ComplexF64}, state_vec)::Vector{ComplexF64}
    d_hilbert_space_shifted = ((N + 1) * (2 * (N + 1) + 1)) # square for two species
    new_vec = zeros(ComplexF64, d_hilbert_space_shifted^2)
    if state_vec ≈ zero(state_vec)
        throw(ArgumentError("Input state is empty"))
    end
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
    nzinds_old = findnz(state_vec)[1]
    if isempty(nzinds_old)
        throw(ArgumentError("Input state is empty"))
    end
    for j_super in nzinds_old
        _shift_j_identical!(N, j_super, state_vec, new_vec)
    end

    return new_vec::SparseVector
end

function _shift_j_identical!(N, j_super, state_vec, new_vec)
    l1, c1, m1, k1, l2, c2, m2, k2 = j_super2lcmk_identical(N, j_super)
    l_shift1 = l1 + c1
    m_shift1 = m1 + k1
    if l_shift1 == m_shift1 && c1 > k1
        # if l_shift1 == m_shift1, we need to swap the indices to maintain the correct order
        j_shift1 = lcmk2j_identical(N + 1, l_shift1, k1, m_shift1, c1)
    else
        j_shift1 = lcmk2j_identical(N + 1, l_shift1, c1, m_shift1, k1)
    end
    l_shift2 = l2 + c2
    m_shift2 = m2 + k2
    if l_shift2 == m_shift2 && c2 > k2
        # if l_shift2 == m_shift2, we need to swap the indices to maintain the correct order
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


"""
    shift_timebins_operator_sp_identical(N)

Compute the time-bin shift operator for a single party of two identical photons for a given
number of time bins `N`. The shift operator is a unitary matrix that describes the action of
the time-bin shift on the two-photon state during one round trip.

See also [`shift_timebins_sp_identical`](@ref), [`shift_timebins_operator_identical`](@ref),
[`coin_operator_sp_identical`](@ref), [`coin_operator_identical`](@ref),
[`mesh_evolution_sp_identical`](@ref), [`mesh_evolution_identical`](@ref).

"""
function shift_timebins_operator_sp_identical(N)
    N = convert(Int64, N)::Int64
    d_hilbert_space = N * (2 * N + 1)
    # simplified expression for N_LOOPS = 2 canceling out geometric series
    d_hilbert_space_shifted = (N + 1) * (2 * (N + 1) + 1)

    shift_op_single_party =
        spzeros(Int64, d_hilbert_space_shifted, d_hilbert_space)

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

    return shift_op_single_party::SparseMatrixCSC{Int64, Int64}
end

"""
    shift_timebins_operator_identical(N)

Two-party version of `shift_timebins_operator_sp_identical`. Compute the tensor product of
two single-party time-bin shift operators for two identical photons for a given number of
time bins `N`.

See also [`shift_timebins_sp_identical`](@ref),
[`shift_timebins_operator_sp_identical`](@ref), [`coin_operator_sp_identical`](@ref),
[`coin_operator_identical`](@ref), [`mesh_evolution_sp_identical`](@ref),
[`mesh_evolution_identical`](@ref).

"""
function shift_timebins_operator_identical(N)
    shift_operator_single_party = shift_timebins_operator_sp_identical(N)
    tensor_shift_operator = kron(shift_operator_single_party, shift_operator_single_party)
    return tensor_shift_operator::SparseMatrixCSC{Int64, Int64}
end


"""
    mesh_evolution_sp_identical(initial_state, angles)

Compute the evolution of a single party of two identical photons for a given initial state
"""
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

function mesh_evolution_identical(initial_state::Vector, angles; kron_mem=nothing)
    state = convert(Vector{ComplexF64}, initial_state)::Vector{ComplexF64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = _iterative_mesh_evolution_identical(state, angles; kron_mem=kron_mem)

    return state::Vector{ComplexF64}
end

function mesh_evolution_identical(initial_state::SparseVector, angles; kron_mem=nothing)
    state =
    convert(SparseVector{ComplexF64, Int64}, initial_state)::SparseVector{ComplexF64, Int64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = _iterative_mesh_evolution_identical(state, angles; kron_mem=kron_mem)

    return state::SparseVector{ComplexF64, Int64}
end

function mesh_evolution_identical(initial_state::AbstractMatrix, angles; kron_mem=nothing)
    state = convert(Matrix{ComplexF64}, initial_state)::Matrix{ComplexF64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    state = _iterative_mesh_evolution_identical(state, angles; kron_mem=kron_mem)

    return state
end


function _iterative_mesh_evolution_identical end

function _iterative_mesh_evolution_identical(
    input_state::AbstractVector, angles; kron_mem=nothing
)
    state = copy(input_state)
    for i in eachindex(angles)
        if isnothing(kron_mem)
            coin_op = coin_operator_identical(angles[i])
        else
            coin_op = coin_operator_identical(angles[i], kron_mem[i])
        end
        state = coin_op * state # apply beam splitters
        state = shift_timebins_identical(state) # shift time bins accordingly
   end

    return state
end

function _iterative_mesh_evolution_identical(
    input_state::AbstractMatrix, angles; kron_mem=nothing
)
    state = copy(input_state)
    for i in eachindex(angles)
        if isnothing(kron_mem)
            coin_op = coin_operator_identical(angles[i])
        else
            coin_op = coin_operator_identical(angles[i], kron_mem[i])
        end
        shift_op = shift_timebins_operator_identical(length(angles[i]))
        state = coin_op * state * coin_op' # apply beam splitters
        state = shift_op * state * shift_op' # apply time-bin shift operator
   end

    return state
end
