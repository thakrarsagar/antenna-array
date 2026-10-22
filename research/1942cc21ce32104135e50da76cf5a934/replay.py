"""Independent standard-library replay of the saved phase-only frontier."""

import cmath
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def power(u, phases):
    # The centered element coordinates add only a unit-magnitude common phase.
    z = cmath.exp(1j * math.pi * u)
    acc = 0j
    for phase in reversed(phases):
        acc = acc * z + cmath.exp(1j * phase)
    return abs(acc / len(phases)) ** 2


def uniform_power(u, n=16):
    if abs(u) < 1e-12:
        return 1.0
    return (math.sin(n * math.pi * u / 2) / (n * math.sin(math.pi * u / 2))) ** 2


def peak(intervals, phases, grid_n=20013):
    n = len(phases)
    # Exact for equally spaced, unit-amplitude elements; independent of phase.
    curvature = math.pi**2 * (n * n - 1) / 6
    lo = hi = 0.0
    for a, b in intervals:
        samples = max(3, math.ceil((b - a) * (grid_n - 1) / 2) + 1)
        h = (b - a) / (samples - 1)
        top = max(power(a + j * h, phases) for j in range(samples))
        lo = max(lo, top)
        hi = max(hi, top + curvature * h * h / 8 + 1e-12)
    return lo, hi


def bisect_crossing(a, b, level, phases):
    fa = power(a, phases) - level
    for _ in range(55):
        mid = (a + b) / 2
        fm = power(mid, phases) - level
        if (fa > 0) == (fm > 0):
            a, fa = mid, fm
        else:
            b = mid
    return (a + b) / 2


def main_peak_and_width(u0, guard, phases):
    a, h = u0 - guard, 2 * guard / 2002
    us = [a + j * h for j in range(2003)]
    ps = [power(u, phases) for u in us]
    k = max(range(len(ps)), key=ps.__getitem__)
    assert 0 < k < len(ps) - 1
    left, right = us[k - 1], us[k + 1]
    for _ in range(60):
        m1 = left + (right - left) / 3
        m2 = right - (right - left) / 3
        if power(m1, phases) < power(m2, phases):
            left = m1
        else:
            right = m2
    peak_u = (left + right) / 2
    level = power(peak_u, phases) / 2
    left_bracket = next(j for j in range(k - 1, -1, -1) if ps[j] <= level)
    right_bracket = next(j for j in range(k + 1, len(ps)) if ps[j] <= level)
    cross_left = bisect_crossing(us[left_bracket], us[left_bracket + 1], level, phases)
    cross_right = bisect_crossing(us[right_bracket - 1], us[right_bracket], level, phases)
    return peak_u, cross_right - cross_left


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    config_path = ROOT / "configs/experiment.json"
    candidates_path = ROOT / "results/candidates.json"
    frontier_path = ROOT / "discover/frontier-table.json"
    config = json.loads(config_path.read_text())
    candidates = {row["id"]: row for row in json.loads(candidates_path.read_text())["rows"]}
    targets = {target["id"]: target for target in config["targets"]}
    ids = ["dev-broadside__uniform__None"] + [r["id"] for r in json.loads(frontier_path.read_text())]
    result = {"source_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in
                               (config_path, candidates_path, frontier_path)}, "rows": []}
    assert config["element_count"] == 16 and config["spacing_wavelengths"] == 0.5
    for row_id in ids:
        row = candidates[row_id]
        target = targets[row["target"]]
        phases = row["phases_rad"]
        assert row["class"] == "phase-only" and row["amplitudes"] == [1.0] * 16
        assert all(abs(x - (j - 7.5) / 2) < 1e-12 for j, x in enumerate(row["positions_lambda"]))
        assert all(abs(phase * 32 / math.pi - round(phase * 32 / math.pi)) < 1e-10 for phase in phases)
        u0 = math.sin(math.radians(target["steer_deg"]))
        guard = config["guard_u"]
        side = [(-1, u0 - guard), (u0 + guard, 1)]
        notches = [(math.sin(math.radians(a)), math.sin(math.radians(b)))
                   for a, b in target["notches_deg"]]
        p0 = power(u0, phases)
        side_lo, side_hi = peak(side, phases)
        notch_lo, notch_hi = peak(notches, phases)
        peak_u, width = main_peak_and_width(u0, guard, phases)
        claimed = row["metrics"]
        assert abs(p0 - claimed["target_power"]) < 1e-10
        assert abs(width - claimed["hpbw_u"]) < 2e-7
        assert abs(peak_u - claimed["peak_u"]) < 2e-6
        for kind, bounds in (("sidelobe", (side_lo, side_hi)), ("notch", (notch_lo, notch_hi))):
            old_lo = p0 * 10 ** (claimed[kind + "_db_lower"] / 10)
            old_hi = p0 * 10 ** (claimed[kind + "_db_upper"] / 10)
            assert bounds[0] <= old_hi + 1e-10 and old_lo <= bounds[1] + 1e-10
        score_lo = max(side_lo, config["notch_weight_power"] * notch_lo) / p0
        score_hi = max(side_hi, config["notch_weight_power"] * notch_hi) / p0
        assert score_lo <= claimed["score_upper"] + 1e-10
        assert claimed["score_lower"] <= score_hi + 1e-10
        assert p0 >= config["minimum_coherent_efficiency"]
        assert width <= config["maximum_hpbw_u"]
        assert abs(peak_u - u0) <= config["maximum_pointing_error_u"]
        if row_id == "dev-broadside__uniform__None":
            for u in (0, 0.1, 0.18, 0.42, 0.9):
                assert abs(power(u, phases) - uniform_power(u)) < 1e-12
            assert abs(10 * math.log10(side_lo / p0) + 13.146835422) < 0.001
        result["rows"].append({"id": row_id, "target_power": p0,
                               "side_power": [side_lo, side_hi],
                               "notch_power": [notch_lo, notch_hi],
                               "score": [score_lo, score_hi],
                               "peak_u": peak_u, "hpbw_u": width,
                               "feasible": True})
    (HERE / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"PASS: independent polynomial replay of {len(ids)} saved witnesses")


if __name__ == "__main__":
    main()
