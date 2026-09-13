"""Figures for the emphasis audit + validation technical report.
Palette: Calypso Pastel (~/System/Viz registry, categorical)."""
import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CAL = ['#178E8E', '#E08D79', '#E8B25C', '#9B8EC4', '#8FAE8B', '#9AA5B1']
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Palatino", "DejaVu Serif"],
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "figure.dpi": 200, "savefig.bbox": "tight",
})

REPO = "/Users/pancho/Code/emphasis"

# ---- Figure 1: audit outcome by theme ------------------------------------
res = json.load(open(f"{REPO}/dev/audit/audit-results.json"))
def theme_of(i):
    n = int(i[1:])
    return ("Likelihood / IS" if n <= 19 else "MCEM / CEM logic" if n <= 43 else
            "Augmentation" if n <= 61 else "Numerics / threads" if n <= 72 else
            "API consistency" if n <= 93 else "Tests / docs")
rows = {}
for v in res["verified"]:
    r = v.get("replication")
    sev = r["severity_final"] if (r and v["verdict"] == "confirmed") else v["severity_rerated"]
    st = "refuted" if v["verdict"] == "refuted" or (r and r["verdict"] == "refuted") else sev
    rows.setdefault(theme_of(v["id"]), []).append(st)
for t in res["triage"]:
    st = "refuted" if t["verdict"] == "refuted" else t["severity_rerated"]
    rows.setdefault(theme_of(t["id"]), []).append(st)

order = ["Likelihood / IS", "MCEM / CEM logic", "Augmentation",
         "Numerics / threads", "API consistency", "Tests / docs"]
cats = [("critical", CAL[1]), ("major", CAL[2]), ("minor", CAL[4]),
        ("not-a-bug", CAL[5]), ("author-question", CAL[3]), ("refuted", "#D8DEE4")]
fig, ax = plt.subplots(figsize=(6.4, 2.9))
left = np.zeros(len(order))
for name, col in cats:
    vals = np.array([sum(1 for s in rows.get(t, []) if s == name) for t in order], float)
    ax.barh(order, vals, left=left, color=col, label=name, height=0.62,
            edgecolor="white", linewidth=0.7)
    for i, (v, l) in enumerate(zip(vals, left)):
        if v >= 2:
            ax.text(l + v/2, i, int(v), ha="center", va="center", fontsize=7.5,
                    color="#2a2a2a")
    left += vals
ax.set_xlabel("hypotheses"); ax.invert_yaxis()
ax.legend(ncol=3, fontsize=7.5, frameon=False, loc="lower right",
          bbox_to_anchor=(1.0, -0.42))
ax.set_title("Audit outcome by subsystem (105 hypotheses)", loc="left")
fig.savefig("figures/fig1-audit-outcomes.pdf"); plt.close(fig)

# ---- Figure 2: TBB thread binding, measured -------------------------------
fig, ax = plt.subplots(figsize=(3.5, 2.6))
req = [1, 4]
before, after = [2, 33], [2, 5]
x = np.arange(len(req)); w = 0.36
ax.bar(x - w/2, before, w, color=CAL[1], label="arena not entered", edgecolor="white")
ax.bar(x + w/2, after, w, color=CAL[0], label="arena entered", edgecolor="white")
for xi, (b, a) in enumerate(zip(before, after)):
    ax.text(xi - w/2, b + 0.7, b, ha="center", fontsize=8)
    ax.text(xi + w/2, a + 0.7, a, ha="center", fontsize=8)
ax.axhline(33, color=CAL[5], lw=0.8, ls=":", zorder=0)
ax.text(1.42, 33.6, "host cores + 1", fontsize=7, color="#666", ha="right")
ax.set_xticks(x); ax.set_xticklabels([f"num_threads = {r}" for r in req])
ax.set_ylabel("peak threads in one E-step"); ax.set_ylim(0, 38)
ax.legend(fontsize=7.5, frameon=False, loc="upper left")
ax.set_title("Requested vs realised parallelism (32-core host)", loc="left", fontsize=9)
fig.savefig("figures/fig2-thread-binding.pdf"); plt.close(fig)
print("wrote fig1-audit-outcomes.pdf, fig2-thread-binding.pdf")
