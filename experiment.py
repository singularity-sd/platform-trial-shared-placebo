"""Platform trials in dermatology: shared controls, seasons and the price of borrowing. A pre-specified simulation.

From the Singularity v3 forecast T3v3-2026-0046 (platform trial -> Dermatology, view v3-classic).
Three treatments for atopic dermatitis enter a platform one after another and are compared with a shared placebo
arm on a binary response (for example EASI-75 at week 16). Disease activity and placebo response vary with the
season. We compare four ways to run or analyse the comparison:
  separate   one two-arm trial per treatment (1:1);
  concurrent platform, compare with controls randomised while the arm was open;
  pooled     platform, compare with all controls up to the end of the arm (non-concurrent included);
  adjusted   platform, all controls, logistic model with a factor for each period between arm entries/exits.
Fixed before the single run: recruitment 40/month, 150 patients per active arm, arms entering at months 0, 6, 12 (no recruitment in a month with no open arm),
placebo response 15% on average, effect on the logit scale 0 or 1.12 (15% -> 35%), four time-trend scenarios,
2000 simulated platforms per scenario and effect, two-sided alpha 0.05, seed.
Run: python experiment.py -> results.json, fig_platform.png
"""

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RNG = np.random.default_rng(20261105)
RATE, N_ARM, ENTRY, REPS = 40, 150, [0, 6, 12], 2000
P0 = 0.15
EFFECTS = {"null": 0.0, "effective": float(np.log(0.35 / 0.65) - np.log(P0 / (1 - P0)))}
SCEN = {
    "no time trend": (0.0, 0.0, False),
    "seasons (±5 points)": (0.4, 0.0, False),
    "seasons + drift": (0.4, 0.6, False),
    "seasons, placebo only": (0.4, 0.0, True),
}
logit = lambda p: np.log(p / (1 - p))
expit = lambda x: 1 / (1 + np.exp(-x))


def simulate_platform(effect, amp, drift, placebo_only):
    """Patient-level platform: in each month, control and open arms share recruitment equally."""
    rows, done, k_open = [], {k: 0 for k in range(3)}, set()
    t = 0.0
    month = 0
    while any(done[k] < N_ARM for k in range(3)):
        k_open = {k for k in range(3) if ENTRY[k] <= month and done[k] < N_ARM}
        if not k_open:
            month += 1
            continue
        groups = ["c"] + sorted(k_open)
        for i in range(RATE):
            g = groups[i % len(groups)]
            if g != "c" and done[g] >= N_ARM:
                continue
            season = amp * np.cos(2 * np.pi * month / 12)
            trend = drift * month / 36
            base = logit(P0) + trend + (0 if placebo_only and g != "c" else season)
            eta = base + (effect if g != "c" else 0.0)
            y = RNG.random() < expit(eta)
            rows.append((month, -1 if g == "c" else g, int(y)))
            if g != "c":
                done[g] += 1
        month += 1
    return np.array(rows)


def logistic_z(y, X, iters=30):
    """Maximum-likelihood logistic regression by Newton-Raphson; returns the Wald z of column 1."""
    b = np.zeros(X.shape[1])
    for _ in range(iters):
        p = 1 / (1 + np.exp(-(X @ b)))
        W = p * (1 - p)
        H = X.T @ (X * W[:, None]) + 1e-9 * np.eye(X.shape[1])
        step = np.linalg.solve(H, X.T @ (y - p))
        b += step
        if np.max(np.abs(step)) < 1e-8:
            break
    cov = np.linalg.inv(H)
    return b[1] / np.sqrt(cov[1, 1])


def wald(y, x, extra=None):
    X = np.column_stack([np.ones(len(y)), x] + ([extra] if extra is not None else []))
    try:
        return bool(abs(logistic_z(y.astype(float), X)) > 1.959964)
    except Exception:
        return False


def analyse(rows, k):
    m, g, y = rows[:, 0], rows[:, 1], rows[:, 2]
    arm = g == k
    a0, a1 = m[arm].min(), m[arm].max()
    ctrl = g == -1
    conc = ctrl & (m >= a0) & (m <= a1)
    pooled = ctrl & (m <= a1)
    out = {}
    s = arm | conc
    out["concurrent"] = (wald(y[s], arm[s].astype(float)), int(conc.sum()))
    s = arm | pooled
    out["pooled"] = (wald(y[s], arm[s].astype(float)), int(pooled.sum()))
    # periods: from one arm entry or exit to the next
    cuts = sorted(set(ENTRY + [int(rows[rows[:, 1] == j, 0].max()) + 1 for j in range(3)]))
    per = np.digitize(m[s], cuts)
    D = np.column_stack([(per == p).astype(float) for p in np.unique(per)[1:]]) if len(np.unique(per)) > 1 else None
    out["adjusted"] = (wald(y[s], arm[s].astype(float), D), int(pooled.sum()))
    return out


def separate(effect, amp, drift, placebo_only, k):
    """A two-arm trial for arm k, 1:1, starting at its entry month, same recruitment rate."""
    n = 2 * N_ARM
    month = ENTRY[k] + np.arange(n) // RATE
    trt = np.arange(n) % 2
    season = amp * np.cos(2 * np.pi * month / 12)
    base = logit(P0) + drift * month / 36 + np.where(placebo_only & (trt == 1), 0, season)
    y = (RNG.random(n) < expit(base + effect * trt)).astype(int)
    return wald(y, trt.astype(float)), N_ARM


def main():
    res = {"protocol": {"RATE": RATE, "N_ARM": N_ARM, "ENTRY": ENTRY, "REPS": REPS, "P0": P0, "effects": EFFECTS}, "results": {}}
    for sname, (amp, drift, po) in SCEN.items():
        res["results"][sname] = {}
        for ename, eff in EFFECTS.items():
            rej = {d: [] for d in ["separate", "concurrent", "pooled", "adjusted"]}
            ctrl = {d: [] for d in rej}
            for _ in range(REPS):
                rows = simulate_platform(eff, amp, drift, po)
                for k in (1, 2):          # the arms that enter later, where borrowing matters
                    a = analyse(rows, k)
                    for d in ("concurrent", "pooled", "adjusted"):
                        rej[d].append(a[d][0]); ctrl[d].append(a[d][1])
                    r, c = separate(eff, amp, drift, po, k)
                    rej["separate"].append(r); ctrl["separate"].append(c)
            res["results"][sname][ename] = {d: {"rejection": float(np.mean(rej[d])), "controls_used": float(np.mean(ctrl[d]))} for d in rej}
            print(sname, ename, {d: round(v["rejection"], 3) for d, v in res["results"][sname][ename].items()}, flush=True)
        (HERE / "results.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    rows = simulate_platform(0.0, 0.0, 0.0, False)
    res["platform_size"] = {"patients": int(len(rows)), "placebo": int((rows[:, 1] == -1).sum()), "months": int(rows[:, 0].max() + 1),
                            "separate_patients": 3 * 2 * N_ARM, "separate_placebo": 3 * N_ARM}
    (HERE / "results.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    figure(res)


def figure(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    D = ["separate", "concurrent", "pooled", "adjusted"]
    sc = list(res["results"])
    fig, axs = plt.subplots(1, 2, figsize=(12, 3.8))
    for ax, e, title in [(axs[0], "null", "False positives (no real effect)"), (axs[1], "effective", "Power (15% → 35%)")]:
        for j, d in enumerate(D):
            ax.bar(np.arange(len(sc)) + (j - 1.5) * 0.2, [res["results"][s][e][d]["rejection"] for s in sc], 0.2, label=d)
        if e == "null":
            ax.axhline(0.05, color="0.3", ls="--", lw=1)
        ax.set_xticks(range(len(sc)))
        ax.set_xticklabels([s.replace(", ", ",\n").replace(" + ", " +\n") for s in sc], fontsize=8)
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.3, axis="y")
    axs[0].set_ylabel("share of comparisons declared significant")
    axs[1].legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(HERE / "fig_platform.png", dpi=160)


if __name__ == "__main__":
    main()
