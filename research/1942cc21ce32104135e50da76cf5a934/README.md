# Phase-only frontier replay

This package replays the seeded `replication-01` card at repository base
`d9db8b4f07222da18b89cd4a936cbcf5e915925f`. It checks the uniform
broadside control and the four saved held-out phase-only frontier witnesses.
The result is an AI agent's author check, not an independent-person review.

Run from the repository root:

```sh
python research/1942cc21ce32104135e50da76cf5a934/replay.py
python -m src.check
python -m pytest -q
```

`replay.py` imports only the Python standard library. It uses a complex
polynomial for the array factor, rather than importing the project's model or
evaluator. The uniform control is also checked against the Dirichlet-kernel
formula. For each witness it checks fixed positions, unit amplitudes, 6-bit
phases, target power, the two sidelobe/notch peak enclosures, weighted score,
beamwidth, pointing, and feasibility. The peak enclosure uses a denser
20,013-point nominal grid and the closed-form curvature bound for a 16-element
unit-amplitude half-wavelength array. The script asserts that its enclosures
overlap the saved ones and writes the measured values and input hashes to
`result.json`.

All five witnesses pass. The four held-out weighted-score intervals from the
independent replay are, in target order, `[0.032081245, 0.032081850]`,
`[0.032984134, 0.032984747]`, `[0.031715947, 0.031716554]`, and
`[0.048638022, 0.048638638]`. The uniform control's fixed-mask sidelobe is
about −13.147 dB. The repository's `src.check` passes all 72 saved
reconstructions, and its 13 tests pass.

The enclosure is a real-arithmetic bound evaluated with ordinary floating
point; it is not directed-rounding interval arithmetic. This replay does not
repeat optimization or perturbation trials, establish hardware performance,
prove global optimality, or review the literature. The frozen result files and
frontier are unchanged.
