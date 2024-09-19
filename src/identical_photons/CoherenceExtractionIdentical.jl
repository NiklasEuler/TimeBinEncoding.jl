export coherence_extraction_identical, combined_measurement_coherence_extraction_identical
export j_out4bins, j_out_hom

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

    Extract the coherences from the final state population `pop_fs` corresponding to the in-
    dex (or vector of indices) `j_out` after a given beam-splitter configuration `angles`.
    Only coherences indicated in `contr_j_tuples` are extracted, all others are bounded and
    then subtracted.

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

    if length(j_out) == 1
        projector_weights = projector_weights[1]
    end

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
		if (j1,j2) in contr_j_tuples &&
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

function j_out4bins(N_post, bin_1, bin_2)
    j_SSSS = lcmk2j_super_identical(N_post, bin_1, 0, bin_1, 0, bin_1, 0, bin_1, 0)
    j_LLLL = lcmk2j_super_identical(N_post, bin_2, 1, bin_2, 1, bin_2, 1, bin_2, 1)
    j_SSLL = lcmk2j_super_identical(N_post, bin_1, 0, bin_1, 0, bin_2, 1, bin_2, 1)
    j_LLSS = lcmk2j_super_identical(N_post, bin_2, 1, bin_2, 1, bin_1, 0, bin_1, 0)
    j_SLSL = lcmk2j_super_identical(N_post, bin_1, 0, bin_2, 1, bin_1, 0, bin_2, 1)

    j_out = [j_SSSS, j_LLLL, j_SSLL, j_LLSS, j_SLSL]
    return j_out
end

function j_out_hom(N_post, bin_1, bin_2)
    # using the Hong-Ou-Mandel interference pattern for coherence extraction
    j_out_lm_SS = lcmk2j_super_identical(N_post, bin_1, 0, bin_1, 0, bin_1, 0, bin_1, 0)
    j_out_lm_LL = lcmk2j_super_identical(N_post, bin_2, 1, bin_2, 1, bin_2, 1, bin_2, 1)
    j_out_lm_hom = lcmk2j_super_identical(N_post, bin_1, 0, bin_2, 1, bin_1, 0, bin_2, 1,)

    j_out = [j_out_lm_SS, j_out_lm_LL, j_out_lm_hom]
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
