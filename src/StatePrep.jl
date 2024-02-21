export correlated_timebin_state, insert_initial_state
export density_matrix, density_matrix_dephased, white_noise
export fidelity, purity



function correlated_timebin_state(wf_coeffs::Vector)#todo: normalize
    wf_coeffs = convert(Vector{ComplexF64}, wf_coeffs)::Vector{ComplexF64}
    N = length(wf_coeffs)
    coeffs = normalize(wf_coeffs)
    time_bin_state_vec = zeros(ComplexF64, N^2)
    for i in 0:N-1
        j = lm2j(N,i,i)
        time_bin_state_vec[j] = coeffs[i+1]
    end
    return time_bin_state_vec
end


function insert_initial_state(time_bin_state_vec::Vector)
    N = Int64(sqrt(length(time_bin_state_vec)))
    full_state_vec = zeros(ComplexF64, length(time_bin_state_vec)*n_loops2)
    for l in 0:N-1, m in 0:N-1
        j_time_bin = lm2j(N, l, m)
        j_full = lcmk2j(N, l, 0, m, 0) #insertion into the short-short ket
        full_state_vec[j_full] = time_bin_state_vec[j_time_bin]
    end
    return full_state_vec
end

function density_matrix(Ψ)
    ρ = kron(Ψ, Ψ')
    return ρ
end

function density_matrix_dephased(Ψ, ϵ)
    ϵ = convert(Float64, ϵ)::Float64
    N = Int64(sqrt(length(Ψ))/n_loops)::Int64
    @argcheck ϵ ≥ 0
    @argcheck ϵ ≤ 1
    
    ρ_pure = density_matrix(Ψ)
    ρ = (1-ϵ) * ρ_pure + ϵ * white_noise(N) 
    return ρ
end

function white_noise(N)
    ρ_noise = zeros(ComplexF64, n_loops2*N^2,n_loops2*N^2)
    weight = 1/N^2
    for l in 0:N-1, m in 0:N-1
        j = lcmk2j(N,l,0,m,0)
        ρ_noise[j,j] = weight
    end
    return ρ_noise
end

function fidelity(Ψ::Vector,ρ::Matrix)
    fidelity = Ψ' * ρ * Ψ
    return convert(Float64,real(fidelity))
end

function purity(ρ)
    return Float64(real(sum(diag(ρ*ρ))))
end