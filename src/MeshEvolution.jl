export shift_timebins, beam_splitter_operator, coin_operator, mesh_evolution
export shift_timebins_single_photon
export phase_on_density_matrix

global const weight_cutoff = 1e-16

function shift_timebins_single_photon(state_vec::Vector)
    state_vec = convert(Vector{ComplexF64}, state_vec)::Vector{ComplexF64}
    new_vec = Vector{ComplexF64}(undef, length(state_vec)+n_loops)
    new_vec[2] = 0
    new_vec[end-1] = 0
    new_vec[1:2:end-3] = @view state_vec[1:2:end]
    new_vec[4:2:end] = @view state_vec[2:2:end]
    return new_vec
end

function shift_timebins(state_vec::Vector)
    state_vec = convert(Vector{ComplexF64}, state_vec)::Vector{ComplexF64}
    N = Int64(sqrt(length(state_vec)/(n_loops2)))
    new_vec = zeros(ComplexF64, ((N+1)*n_loops)^2)
    for j in eachindex(state_vec)
        l,c,m,k = j2lcmk(N,j)
        shifted_j = lcmk2j(N+1, l+c, c, m+k, k) # adapted system has one more time bin, so we need to put N+1
        new_vec[shifted_j] = state_vec[j]
    end
    return new_vec
end

function shift_timebins(state_vec::SparseVector)
    state_vec = convert(SparseVector{ComplexF64, Int64}, state_vec)::SparseVector{ComplexF64, Int64}
    N = Int64(sqrt(length(state_vec)/(n_loops2)))
    new_vec = spzeros(ComplexF64, ((N+1)*n_loops)^2)
    for j in state_vec.nzind
        l,c,m,k = j2lcmk(N,j)
        shifted_j = lcmk2j(N+1, l+c, c, m+k, k) # adapted system has one more time bin, so we need to put N+1
        new_vec[shifted_j] = state_vec[j]
    end
    return new_vec
end

function beam_splitter_operator(θ)
    θ = convert(Float64,θ)::Float64
    cs = cos(θ)
    sn = im*sin(θ)
    cols = [1,1,2,2]
    rows = [1,2,1,2]
    vals = [cs, sn, sn, cs]
   return sparse(cols,rows,vals)
end

function coin_operator(angles::Vector)
    real_angles = convert(Vector{Float64}, angles)::Vector{Float64}
    matrices = [beam_splitter_operator(θ) for θ in real_angles]
    single_photon_coin_operator = blockdiag(matrices...)
    tensor_coin_operator = kron(single_photon_coin_operator,single_photon_coin_operator)
    return tensor_coin_operator::SparseMatrixCSC{ComplexF64, Int64}
end

"""
    mesh_evolution(ψ_init, angles)

TBW
"""
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
        state = coin_op * state
        state = shift_timebins(state)
    end
    return state
end

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