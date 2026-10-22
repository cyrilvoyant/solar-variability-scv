"""Figures of the sCV / F demonstration, drawn from the CSV files written by demo_scv.m.

    python make_figures.py

Writes figures/*.png and figures/*.pdf. Requires numpy, pandas and matplotlib.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RES, FIG = HERE / "results", HERE / "figures"
FIG.mkdir(exist_ok=True)
D = 0.2
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10.5, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#dfe5ec",
    "grid.linewidth": 0.8, "axes.axisbelow": True, "legend.frameon": False, "savefig.dpi": 200,
})
BLUE, RED, GOLD, GREEN, INK, GREY = "#2e86ab", "#d1495b", "#e0a03a", "#3c9d5d", "#1b2a3a", "#9aa8b6"


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# 1. Four skies on the same clear-sky cycle
s = pd.read_csv(RES / "skies_series.csv")
m = pd.read_csv(RES / "skies_metrics.csv").set_index("sky")
titles = {"clear": "Clear sky", "light": "Light, short-lived clouds", "broken": "Broken clouds",
          "reference": "Reference case: clear or diffuse floor, 50/50"}
fig, axes = plt.subplots(2, 2, figsize=(10, 5.6), sharex=True, sharey=True)
for ax, key in zip(axes.ravel(), titles):
    ax.plot(s["hour"], s["clear_sky"], color=INK, lw=1.2, label="Clear sky $I_{clr}$")
    ax.plot(s["hour"], D * s["clear_sky"], color=INK, lw=1, ls=":", label="Diffuse floor $d\\,I_{clr}$, d = 0.2")
    ax.fill_between(s["hour"], s[key], color=BLUE, alpha=0.35, lw=0)
    ax.plot(s["hour"], s[key], color=BLUE, lw=1.4, label="Measured $I$")
    r = m.loc[key]
    ax.set_title(titles[key])
    ax.text(0.335, 0.97, f"CV = {r.CV:.2f}\nsCV = {r.sCV:.2f}\nF (1 h) = {r.F_1h:.2f}", transform=ax.transAxes,
            ha="center", va="top", fontsize=9.5, bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#ccd6e0"))
    ax.set_xticks(range(0, 73, 12))
for ax in axes[1]:
    ax.set_xlabel("Hour")
for ax in axes[:, 0]:
    ax.set_ylabel("Irradiance (W/m²)")
h, lab = axes[0, 0].get_legend_handles_labels()
fig.legend(h, lab, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 0.99), fontsize=9.5)
fig.tight_layout()
save(fig, "fig1_four_skies")

# 2. Forecast error against F and against sCV, Monte Carlo
mc = pd.read_csv(RES / "montecarlo.csv")
fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)
cols = {1: BLUE, 4: GOLD, 10: RED}
for ax, xname, xlabel in [(axes[0], "F", "Forecastability F"), (axes[1], "sCV", "Variability sCV")]:
    for hz, c in cols.items():
        sub = mc[mc.horizon == hz]
        rho = sub[[xname, "nRMSE_SP"]].corr(method="spearman").iloc[0, 1]
        ax.scatter(sub[xname], sub["nRMSE_SP"], s=16, color=c, alpha=0.75,
                   label=f"{hz} step{'s' if hz > 1 else ''} ahead, Spearman ρ = {rho:+.2f}")
    ax.set_xlabel(xlabel)
    ax.legend(loc="upper right" if xname == "F" else "upper left", fontsize=8.5)
axes[0].set_ylabel("nRMSE of smart persistence")
axes[0].set_title("Error falls as F rises")
axes[1].set_title("Error grows with sCV, less sharply")
save(fig, "fig2_error_vs_F_and_sCV")

# 3. F with the horizon; sCV does not depend on it
g = mc.groupby("horizon")
fig, ax = plt.subplots(figsize=(7.6, 3.9))
for name, c, lab in [("F", BLUE, "Forecastability F"), ("sCV", RED, "Variability sCV")]:
    med, q1, q3 = g[name].median(), g[name].quantile(0.25), g[name].quantile(0.75)
    ax.plot(med.index, med, color=c, lw=2, marker="o", label=f"{lab} (median of 100 series)")
    ax.fill_between(med.index, q1, q3, color=c, alpha=0.18, lw=0)
ax.set_xlabel("Horizon n (time steps)")
ax.set_ylabel("Value")
ax.set_ylim(0, 1.02)
ax.set_xticks([1, 4, 7, 10])
ax.legend(ncol=2, loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=9)
save(fig, "fig3_F_and_sCV_by_horizon")
print("figures written to", FIG)
