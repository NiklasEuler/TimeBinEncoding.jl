export N_LOOPS, N_LOOPS2, WEIGHT_CUTOFF, Θ_13, Θ_23, θ_pop_ref, _cols, _rows

"""
Number of fiber loops. Saved as const to avoid magic numbers.
"""
const global N_LOOPS = 2

"""
Number of fiber loops squared. Saved as const to avoid magic numbers.
"""
const global N_LOOPS2 = 4

"""
Numerical cutoff value to determine wether a certain coefficient is ≈ 0
"""
global const WEIGHT_CUTOFF = 1e-16

"""
Angle such that sin(Θπ) = 1/√3, used in simplified beam-splitter settings. Saved as const to
avoid magic numbers.
"""
global const Θ_13 = asin(sqrt(1 / 3)) / π

"""
Angle such that sin(Φπ) = 2/√3, used in simplified beam-splitter settings. Saved as const to
avoid magic numbers.
"""
global const Θ_23 = asin(sqrt(2 / 3)) / π

"""
Angle such that sin(Θπ) = √3/4, used for the population reference beam-splitter settings.
Saved as const to avoid magic numbers.
"""
global const θ_pop_ref = asin(sqrt(3 / 4)) / π

"""
beam-splitter colums and rows.
"""
global const _cols = [1, 1, 2, 2]
global const _rows = [1, 2, 1, 2]
