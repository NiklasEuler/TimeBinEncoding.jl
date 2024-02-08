
function shift_timebins(state_vec::Vector)
    new_vec = Vector{ComplexF64}(undef, length(state_vec)+4)
    new_vec[2] = 0
    new_vec[end-1] = 0
    new_vec[1:2:end-3] = state_vec[1:2:end]
    new_vec[2:2:end-1] = state_vec[2:2:end]
    return new_vec
end

function beam_splitter_operator(θ)
    
    cs = cos(θ)
    sn = im*sin(θ)
    cols = [1,1,2,2]
    rows = [1,2,1,2]
    vals = [cs, sn, sn, cs]
   # sparse([1, 1, 2, 3], [1, 3, 2, 3], [0, 1, 2, 0])
function coin_operator