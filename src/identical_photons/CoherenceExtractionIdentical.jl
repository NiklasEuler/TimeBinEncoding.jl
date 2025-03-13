export coherence_extraction_identical, combined_measurement_coherence_extraction_identical
export j_out4bins, j_out_hom
export combined_weights_pops_4bins_all, combined_weights_pops_hom, combs_4bin
export combined_projector_weights_auto

"""
    coherence_extraction_identical(
        N,
        j_out,
        pops_init,
        pop_fs,
        angles,
        contr_j_tuples=correlated_short_bins_tuples_identical(N; extract_diagonal=false),
        projector_weights=ones(Float64, length(j_out))
    )

Extract the coherences from the final state population `pop_fs` corresponding to the index
(or vector of indices) `j_out` after a given beam-splitter configuration `angles`. Only
coherences indicated in `contr_j_tuples` are extracted, all others are bounded and then
subtracted.

"""
function coherence_extraction_identical(
    N,
    j_out,
    pops_init,
    pop_fs,
    angles,
    contr_j_tuples=correlated_short_bins_tuples_identical(N; extract_diagonal=false),
    projector_weights=ones(Float64, length(j_out))
)
    # original default was extract_diagonal=true. If problems arise, explicitly set to true.
    N = convert(Int64, N)::Int64
    j_out = try
		convert(Vector{Int64}, j_out)::Vector{Int64}
	catch
		convert(Int64, j_out)::Int64
	end

 #=    if length(j_out) == 1
        projector_weights = projector_weights[1]
    end =#

    projector_weights = try
		convert(Vector{Float64}, projector_weights)::Vector{Float64}
	catch
		convert(Float64, projector_weights)::Float64
	end
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}

    j1_arr, j2_arr, weights =
        explicit_fs_coherence_map_identical(j_out, angles, projector_weights)
        # extraction based on assumed ideal angles
    @argcheck weights ≠ []

	extracted_coherence = []
    extracted_weights = Float64[]
    for idx in eachindex(j1_arr)
		j1 = j1_arr[idx]
		j2 = j2_arr[idx]
		if (j1,j2) in contr_j_tuples
            # check whether both j1 and j2 are correlated time bins
			push!(extracted_coherence, (j1, j2))
            # relevant coherence, so indices saved to list of extracted coherences.
            push!(extracted_weights, weights[idx])
		elseif j1 == j2
			pop_fs -= pops_init[j1] * weights[idx]
            # non-contributing population. Can be removed exactly as exact value is known.
		else
			pop_fs -= sqrt(pops_init[j1] * pops_init[j2]) * abs(weights[idx])
            if sqrt(pops_init[j1] * pops_init[j2]) * abs(weights[idx]) != 0

                #println(
                    #"Subtracting: (",j1, " ", j2, ") ", sqrt(pops_init[j1]*pops_init[j2]) *
                    #abs(weights[idx])
                #)
            end
            # subtract non-contributing coherence bound.
		end
	end

    @argcheck extracted_weights ≠ []

    n_contr_sched = length(contr_j_tuples)
    n_extracted = length(extracted_weights)

    if (n_extracted != n_contr_sched )
        @warn "Some of the scheduled coherences have a vanishing weight in the given "*
        "final-state projectors. Please check again and consider adapting the scheduled "*
        "coherences in `contr_j_tuples`."
    end

    if !all(extracted_weights .≈ extracted_weights[1])
        @warn "The coherences scheduled for extraction have differing "*
        "weights in the chosen final-state projectors and can thus not straight-forwardly "*
        "be extracted. Please check input again."
        norm = 1 #minimum(extracted_weights) #maximum(extracted_weights)
        # solution for now, but not ideal.
        # Better solution would compute sensible normalization.
    else
        norm = extracted_weights[1]
    end


	pop_fs /= norm # normalization of the extracted coherences

	return convert(Float64, pop_fs)#, extracted_coherence
end

"""
    j_out4bins(N_post, bin_1, bin_2)

    Compute the `j_out` projector-state four-photon indicies for a four-bin interference
    setup. The photons can be detected in bins `bin_1` and `bin_2` out of `N_post` bins
    after the DTQW. The computed outputs are the indices for the full coincidence (SSSS,
    LLLL), species coincidence (SSLL, LLSS), and signal-idler coincidence (SLSL) configura-
    tions.

    See also [`j_out_hom`](@ref).
"""
function j_out4bins(N_post, bin_1, bin_2)
    j_SSSS = lcmk2j_super_identical(N_post, bin_1, 0, bin_1, 0, bin_1, 0, bin_1, 0)
    j_LLLL = lcmk2j_super_identical(N_post, bin_2, 1, bin_2, 1, bin_2, 1, bin_2, 1)
    j_SSLL = lcmk2j_super_identical(N_post, bin_1, 0, bin_1, 0, bin_2, 1, bin_2, 1)
    j_LLSS = lcmk2j_super_identical(N_post, bin_2, 1, bin_2, 1, bin_1, 0, bin_1, 0)
    j_SLSL = lcmk2j_super_identical(N_post, bin_1, 0, bin_2, 1, bin_1, 0, bin_2, 1)

    j_out = [j_SSSS, j_LLLL, j_SSLL, j_LLSS, j_SLSL]
    return j_out
end

"""
    j_out_hom(N_post, bin_1, bin_2)

    Compute the `j_out` projector-state four-photon indicies for a Hong-Ou-Mandel interfer-
    ence setup. The photons can be detected in bins `bin_1` and `bin_2` out of `N_post` bins
    after the DTQW. The computed outputs are the indices for the full coincidence (SSSS,
    LLLL) and signal-idler coincidence (SLSL) configurations.

    See also [`j_out4bins`](@ref).
"""
function j_out_hom(N_post, bin_1, bin_2)
    # using the Hong-Ou-Mandel interference pattern for coherence extraction
    @argcheck bin_1 < bin_2
    @argcheck bin_1 ≥ 0
    @argcheck bin_2 < N_post

    j_out_lm_SSSS = lcmk2j_super_identical(N_post, bin_1, 0, bin_1, 0, bin_1, 0, bin_1, 0)
    j_out_lm_LLLL = lcmk2j_super_identical(N_post, bin_2, 1, bin_2, 1, bin_2, 1, bin_2, 1)
    j_out_lm_SLSL = lcmk2j_super_identical(N_post, bin_1, 0, bin_2, 1, bin_1, 0, bin_2, 1,)

    j_out = [j_out_lm_SSSS, j_out_lm_LLLL, j_out_lm_SLSL]
    return j_out
end

"""
    combined_measurement_coherence_extraction_identical(
        N,
        combined_weights,
        pops_init,
        pop_fs_combined,
        contr_j_tuples=correlated_short_bins_tuples_identical(N);
    )

    Extract the coherences from a combined measurement of several different beam-splitter
    configurations or phase settings. This allowes for the cancelation of additional cohe-
    rences across different measurements and improves fidelity of the extracted coherences.

"""
function combined_measurement_coherence_extraction_identical(
    N,
    combined_weights,
    pops_init,
    pop_fs_combined,
    contr_j_tuples=correlated_short_bins_tuples_identical(N; extract_diagonal=false)
)
    d_full_hilbert_space = lcmk2j_super_identical(N, N-1, 1, N-1, 1, N-1, 1, N-1, 1)
    # full four-photon Hilbert space dimenion
	#extracted_coherence = []
    pop_fs = copy(pop_fs_combined)
    extracted_weights = Float64[]
    for j in combined_weights.nzind
		j1, j2 = j2lm(d_full_hilbert_space, j)
		j1 += 1
		j2 += 1
		if (j1,j2) in contr_j_tuples # check outsurced to contr_j_tuples creation
            push!(extracted_weights, combined_weights[j])
		elseif j1 == j2
			pop_fs -= pops_init[j1] * combined_weights[j]
            # non-contributing population. Can be removed exactly as exact value is known.
		else
			pop_fs -= sqrt(pops_init[j1] * pops_init[j2]) * abs(combined_weights[j])
            # subtract non-contributing coherence bound.
		end
	end

    @argcheck extracted_weights ≠ []

    n_contr_sched = length(contr_j_tuples)
    n_extracted = length(extracted_weights)

    if (n_extracted != n_contr_sched )
        @warn "Some of the scheduled coherences have a vanishing weight in the given "*
        "final-state projectors. Please check again and consider adapting the scheduled "*
        "coherences in `contr_j_tuples`."
    end

    if !all(extracted_weights .≈ extracted_weights[1])
        @warn "The coherences scheduled for extraction have differing "*
        "weights in the chosen final-state projectors and can thus not straight-forwardly "*
        "be extracted. Please check input again."
        norm = 1 #minimum(extracted_weights) #maximum(extracted_weights)
        # solution for now, but not ideal.
        # Better solution would compute sensible normalization.
    else
        norm = extracted_weights[1]
    end


	pop_fs /= norm # normalization of the extracted coherences

	return convert(Float64, pop_fs)#, extracted_coherence
end


"""
    combined_projector_weights_auto(N, extraction_weights=[1, 1, 1, 1, 1])

Compute the weights for the different projectors in the four-bin and Hong-Ou-Mandel
interference setup. The weights are computed automatically based on the number of time bins
`N` and the desired weights for the different coherences given in `extraction_weights`.

The `extraction weights` are given in the following order:
- dimer-no-dimer coherences with two different bins (e.g. ⟨iiii|ρ|ijij⟩)
- dimer-no-dimer coherences with three different bins, (e.g. ⟨iiii|ρ|jkjk⟩)
- dimer-dimer coherences with two different bins, (e.g. ⟨iiii|ρ|jjjj⟩)
- non-dimer-non-dimer coherences with four different bins, (e.g. ⟨ijij|ρ|klkl⟩)
- non-dimer-non-dimer coherences with three different bins. (e.g. ⟨ijij|ρ|jkjk⟩)

# Returns
- `Vector{Float64}`: The weights for the different projectors.
    The first three elements correspond to the four-bin interference and refere to the full
    coincidence, species coincidence, and signal-idler coincidence, respectively. The last
    two elements correspond to the Hong-Ou-Mandel interference and refer to the full
    coincidence and signal-idler coincidence, respectively.

See also [`combined_weights_pops_4bins_all`](@ref), [`combined_weights_pops_hom`](@ref),
[`j_out4bins`](@ref), [`j_out_hom`](@ref).

"""
function combined_projector_weights_auto(N, extraction_weights=[1, 1, 1, 1, 1])

    weights_single = 1 / 32 * reshape([
        [6, 6, 3, 12, 12];;
        [-2, -2, 3, 12, -4];;
        [4, 4, 6, 8, 0];;
        [0, 0, 8, 0, 0];;
        [8, 0, 4, 0, 0];;
        ], (5,5)
    )
	# columns meaning left to right:
    # weight of (wo) dimer-no-dimer (dnd) coherences with two different bins,
    # wo dnd coherences with three different bins,
    # wo dimer-dimer (dd) coherences with two different bins,
    # wo non-dimer-non-dimer (ndnd) coherences with four different bins,
    # wo ndnd coherences with three different bins.
	# first three lines: coming from 4bin configuration. Some terms appear in more than one
    # 4bin setting. Thus, their weights are multiplied by a scaling factor, which we compute
    # in the following. last two lines: from compound configuration. These appear once each
    # and thus have no multiple appearence problem.

    scaling_dd = scaling_dnd_one_diff = binomial(N - 2, 2)
            # number of 4bin configurations that include a given dd coherence
	scaling_ndnd_one_diff = scaling_dnd_all_diff = binomial(N - 3, 1)
        # number of 4bin configurations that include a given ndnd coherence with three
        # different bins
	scaling_ndnd_all_diff = 1
        # number of 4bin configurations that include a given ndnd coherence with 4 different
        # bins

	weights_single[1, 1:3] *= scaling_dnd_one_diff
	weights_single[2, 1:3] *= scaling_dnd_all_diff
	weights_single[3, 1:3] *= scaling_dd
	weights_single[4, 1:3] *= scaling_ndnd_all_diff
	weights_single[5, 1:3] *= scaling_ndnd_one_diff

	weights_single_inv = inv(weights_single)
	projector_weigths_auto = weights_single_inv * extraction_weights
    return projector_weigths_auto
end

"""
    combs_4bin(N)

Generates all possible combinations of 4 bins out of N bins, such that the bins are
ordered in increasing order. The output is a 4xN_comb matrix, where each column
corresponds to a combination of 4 bins.

"""
function combs_4bin(N)
	@argcheck N ≥ 4
	N_comb = binomial(N, N - 4)
	bin_idxs = zeros(Int64, 4, N_comb)
	i = 1
	for a in 0:N - 4, b in a + 1:N - 3, c in b + 1:N - 2, d in c + 1:N - 1
		bin_idxs[:, i] = [a, b, c, d]
		i += 1
	end
	return bin_idxs
end
"""
    combined_weights_pops_4bins_all(N, extraction_composition_weights; phases=true)

Compute the combined weights and populations for the three four-bin interference settings,
with the desired weights for the different projectors given in
`extraction_composition_weights`. Return both the combined weights and the total final-state
populations.

See also [`combined_weights_pops_hom`](@ref), ['angles4bins](@ref).

"""
function combined_weights_pops_4bins_all(
    N, ρ_mixed, extraction_composition_weights; phases=true
    )

	projector_weights_4b = [
        extraction_composition_weights[1],
        extraction_composition_weights[1],
        extraction_composition_weights[2],
        extraction_composition_weights[2],
        extraction_composition_weights[3]

    ] # projector weights for the four-bin interference
	d_local_hs_bl = N * (2 * N + 1)
	d_full_hs_bl = d_local_hs_bl ^ 2

    combined_weights_4b = spzeros(ComplexF64, d_full_hs_bl^2)
    # Agnostic to the actual state ρ_mixed, contains the weights of the coherences of
    # ρ_mixed. No noise here, as this is the ideal planned DTQW.
    pop_fs_4b = 0.0
    # Above term, but with explicit values of ρ plugged in. In experiment, this is subjected
^   # to imperfect DTQW trajectories and sampling noise.


	bin_combinations = combs_4bin(N)
	N_comb = binomial(N, N - 4) # only ever let 4 bins interfere
	phase_args = zeros(N, 4)
    for j in 1:N_comb
        a, b, c, d = bin_combinations[:, j]
        angles_all = angles4bins(N, a, b, c, d)
        phase_args .= 0
        if phases
            phase_args[a + 1, 1] = 1
            phase_args[b + 1, 2] = 1
            phase_args[c + 1, 3] = 1
            phase_args[d + 1, 4] = 1
        end
        phases_all = cispi.(phase_args)
        for angles in angles_all
            N_post = N + length(angles)
            j_out = j_out4bins(N_post, d, d + 1)
            for phases in eachcol(phases_all)
                j1_arr, j2_arr, weights = explicit_fs_coherence_map_identical(
                    j_out, angles, projector_weights_4b, phases
                )
                pop_fs_4b += explicit_fs_pop_identical(
                    ρ_mixed, j_out, angles, projector_weights_4b, phases
                )
                for i in eachindex(j1_arr)
                    j1 = j1_arr[i]
                    j2 = j2_arr[i]
                    j_comb = lm2j(d_full_hs_bl, j1 - 1, j2 - 1)
                    combined_weights_4b[j_comb] += weights[i]
                end
            end
		end
	end

	return combined_weights_4b, pop_fs_4b
end

"""
    combined_weights_pops_hom(N, ρ_mixed, extraction_composition_weights)

    Compute the combined weights and populations for the Hong-Ou-Mandel interference
    setting, with the desired weights for the different projectors given in
    `extraction_composition_weights`. Return both the combined weights and the total
    final-state populations.

    See also [`combined_weights_pops_4bins_all`](@ref), [`angles_compound`](@ref).
"""
function combined_weights_pops_hom(N, ρ_mixed, extraction_composition_weights)

	pop_fs_hom = 0.0
	d_local_hs_bl = N * (2 * N + 1)
	d_full_hs_bl = d_local_hs_bl ^ 2

	projector_weights_hom = [
        extraction_composition_weights[4],
        extraction_composition_weights[5],
        extraction_composition_weights[5]
    ]
	combined_weights_hom = spzeros(ComplexF64, d_full_hs_bl^2)

	for l in 0:N - 2
		for m in l + 1:N - 1
			phase_args = zeros(N, 1)
			#phase_args[m + 1, 0] = 0.5 # once with phases, once without
			phases_hom_all = cispi.(phase_args)
			angles_distance_k = angles_kth_neighbor_interference(N, m - l)
			angles_lm = angles_distance_k[l + 1]
			N_post_lm = N + length(angles_lm)
			j_out_SSSS = lcmk2j_super_identical(N_post_lm, m, 0, m, 0, m, 0, m, 0)
			j_out_LLLL =
                lcmk2j_super_identical(N_post_lm, m + 1, 1, m + 1, 1, m + 1, 1, m + 1, 1)

			j_out_lm_hom = lcmk2j_super_identical(N_post_lm, m, 0, m + 1, 1, m, 0, m + 1, 1)
			j_out = [j_out_lm_hom, j_out_SSSS, j_out_LLLL]

			for phases in eachcol(phases_hom_all)
				j1_arr_hom, j2_arr_hom, weights_hom = explicit_fs_coherence_map_identical(
                    j_out, angles_lm, projector_weights_hom, phases
                )
				pop_fs_hom += explicit_fs_pop_identical(
                    ρ_mixed, j_out, angles_lm, projector_weights_hom, phases
                )
				for i in eachindex(j1_arr_hom)
					j1 = j1_arr_hom[i]
					j2 = j2_arr_hom[i]
					j_comb = lm2j(d_full_hs_bl, j1 - 1, j2 - 1)
					combined_weights_hom[j_comb] += weights_hom[i]
				end
			end
		end
	end

	return combined_weights_hom, pop_fs_hom
end
