export N_LOOPS, N_LOOPS2, WEIGHT_CUTOFF, Θ, Φ

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
global const Θ = asin(sqrt(1/3))/π

"""
Angle such that sin(Φπ) = 2/√3, used in simplified beam-splitter settings. Saved as const to
avoid magic numbers.
"""
global const Φ = asin(sqrt(2/3))/π
