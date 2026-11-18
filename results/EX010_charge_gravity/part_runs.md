# EX010 inline drivers (as run, 2026-11-18)
Part 1: solve "−" (base + minus(0.005)); fsolve over log(kf), log(eps) of base + plus(kf, eps) for Q+ = -Q-, mu_C+ = -mu_C-,
start (0.0731, 1.4); result written to part1_tuned.json.
Part 2: for kg in (0, 0.005, 0.05, 0.5) and type in (plus tuned, minus(0.005), neutral): solve base + type + grav(kg);
record Q = net production of C, M = -net production of G, drifts for C and G; written to part2_combined.json.
