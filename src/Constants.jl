export N_LOOPS, N_LOOPS2, WEIGHT_CUTOFF

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
