export coherence_extraction, initial_state_phase_estimation


function coherence_extraction(N, j_out, ρ, angles, extract_diagonal=true)
    N = convert(Int64, N)::Int64
    j_out = try 
		convert(Vector{Int64}, j_out)::Vector{Int64}
	catch 
		convert(Int64, j_out)::Int64
	end
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    pops = Float64.(diag(ρ))



    contr_j_idxs = correlated_short_bins_idxs(N)
    j1_arr, j2_arr, weights = explicit_final_state_coherence_map(j_out, angles)
    @argcheck weights ≠ []

	pop_j_out_extracted = explicit_final_state_projection_expval(ρ, j_out, angles)
	extracted_coherence = []
	# coherence_extracted = 0.0
	for idx in eachindex(j1_arr)
		j1 = j1_arr[idx]
		j2 = j2_arr[idx]
		if(j1 ∈ contr_j_idxs && j2 ∈ contr_j_idxs && (extract_diagonal || j1 ≠ j2))
			push!(extracted_coherence,(j1,j2))
		elseif j1 == j2
			pop_j_out_extracted -= pops[j1] * weights[idx]
		else			
			pop_j_out_extracted -= sqrt(pops[j1]*pops[j2]) * abs(weights[idx])
		end
	end
	pop_j_out_extracted /= N*weights[1]
	return convert(Float64, pop_j_out_extracted)#, extracted_coherence
end

function initial_state_phase_estimation(ρ_init)
    ρ = convert(Matrix{ComplexF64}, copy(ρ_init))::Matrix{ComplexF64}
    N = Int64(sqrt(size(ρ)[1]/(n_loops2)))
    
    ϕ_arr = (0:0.0001:2)*π
    #contr_j_idxs = correlated_short_bins_idxs(N)
    #contr_pops = Float64(sum([ρ[j,j] for j in contr_j_idxs]))
    #population_correction_term = contr_pops/N
    nn_phases = zeros(Float64, N)
    #tb_coefficients = [state[lcmk2j(N,i,0,i,0),lcmk2j(N,i,0,i,0)] for i in 0:N-1]
    
    angles_1_1 = zeros(Float64,N)
    angles_1_1[1:2:end] .= 0.5*π
    angles_1_2 = zeros(Float64,N+1)
    angles_1_2[2:2:end] .= 0.25*π
    angles_1 = [angles_1_1,angles_1_2]

    angles_2_1 = zeros(Float64,N)
    angles_2_1[2:2:end] .= 0.5*π
    angles_2_2 = zeros(Float64,N+1)
    angles_2_2[3:2:end] .= 0.25*π
    angles_2 = [angles_2_1,angles_2_2]

    angles_arr = [angles_1, angles_2]
    j_out_1_arr = [[lcmk2j(N+2,i,0,i,0),lcmk2j(N+2,i+1,1,i+1,1)] for i in 1:2:N-1]
    #if iseven(N)
        j_out_2_arr = [[lcmk2j(N+2,i,0,i,0),lcmk2j(N+2,i+1,1,i+1,1)] for i in 2:2:N-1]
    #else
    #    j_out_2_arr = [[lcmk2j(N+2,i,0,i,0),lcmk2j(N+2,i+1,1,i+1,1)]for i in 2:2:N-1]
    #end
    j_out_arr = [j_out_1_arr, j_out_2_arr]
    for l ∈ [1,2]
        j_out = j_out_arr[l]
        angles = angles_arr[l]
        for (k, j) in enumerate(j_out)
            φ_arr = zeros(Float64, N)
            bin_idx = 2*(k-1)+l+1
            φ_arr[bin_idx] = π/2
            ρ_rotated = phase_on_density_matrix(ρ, φ_arr)
			c_real = coherence_extraction(N, j, ρ, angles, false)
			c_imag = coherence_extraction(N, j, ρ_rotated, angles, false)
			c_contr = c_real .* cos.(-ϕ_arr) .+ c_imag .* sin.(-ϕ_arr)
			nn_phases[bin_idx] = ϕ_arr[argmax(c_contr)]
		end
    end
    relative_phases =  mod.(cumsum(nn_phases),2*π)
    ρ_corrected = phase_on_density_matrix(ρ, -1 * relative_phases)
    return ρ_corrected, relative_phases
end

