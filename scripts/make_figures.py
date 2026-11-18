"""Make the README figures from the CSV files in results/.

Usage (from the repository root, after `run_upper_airway` in MATLAB/Octave):
    pip install numpy pandas matplotlib
    python scripts/make_figures.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = ROOT / "docs" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# Palette (categorical slots in fixed order) and ink
BLUE, ORANGE, AQUA, MAGENTA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#e87ba4", "#6250d6"
INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8986", "#e6e5e1", "#fcfcfb"
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#eef4fc", "#9cc3ef", "#2a78d6", "#0d3a75"])

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 10, "axes.titlesize": 11.5, "axes.titleweight": "bold",
    "axes.titlecolor": INK, "axes.titlelocation": "left", "legend.frameon": False,
    "lines.linewidth": 2,
})

ts = pd.read_csv(RES / "timeseries.csv")
gen = pd.read_csv(RES / "generations.csv")
summ = pd.read_csv(RES / "summary.csv").set_index("quantity")["value"]
co2 = pd.read_csv(RES / "co2_map.csv")
axis = pd.read_csv(RES / "co2_map_axis.csv")
t = ts["t_s"].to_numpy()
t_switch = t[np.argmax(ts["Q_Lps"].to_numpy()[1:] < 0) + 1]


def phase_shade(ax):
    ax.axvspan(t[0], t_switch, color=BLUE, alpha=0.05, lw=0)
    ax.set_xlim(t[0], t[-1])


# ---------------------------------------------------------------- 1. overview
fig, axs = plt.subplots(2, 2, figsize=(11, 7.2), constrained_layout=True)

ax = axs[0, 0]
phase_shade(ax)
ax.plot(t, ts["Q_Lps"] * 1e3, color=BLUE)
ax.axhline(0, color=MUTED, lw=1)
ax.set(title="Airflow at the trachea", xlabel="Time in breath (s)", ylabel="Flow (mL/s)")
ax.text(0.03, 0.06, "inspiration", transform=ax.transAxes, color=INK2)
ax.text(0.53, 0.94, "expiration", transform=ax.transAxes, color=INK2, va="top")

ax = axs[0, 1]
phase_shade(ax)
ax.plot(t, ts["drop_trachea_mmHg"], color=BLUE, label="Mouth to trachea (nose-larynx)")
ax.plot(t, ts["drop_gen4_out_mmHg"], color=ORANGE, label="Mouth to Gen 4 exit")
ax.plot(t, -ts["P_alv_gauge_mmHg"], color=AQUA, label="Mouth to alveoli")
ax.axhline(0, color=MUTED, lw=1)
ax.set(title="Pressure lost along the airway", xlabel="Time in breath (s)",
       ylabel="P$_{atm}$ − P (mmHg)")
ax.legend(loc="upper right", fontsize=8.5)

ax = axs[1, 0]
vol = axis["volume_from_mouth_mL"].to_numpy()
Z = co2.drop(columns="t_s").to_numpy().T
im = ax.pcolormesh(t, vol, Z, cmap=SEQ, shading="auto", vmin=0, vmax=40)
ax.grid(False)
v_tr = vol[axis["region"].to_numpy() == 1][0] - 0.375
v_g4 = vol[axis["region"].to_numpy() == 1][-1] + 0.375
for v in (v_tr, v_g4):
    ax.axhline(v, color=ORANGE, lw=1.2, ls="--")
ax.text(t[-1] - 0.05, v_tr / 2, "nose–larynx", ha="right", va="center", color="white", fontsize=8.5, fontweight="bold")
ax.text(t[-1] - 0.05, (v_tr + v_g4) / 2, "Gens 0–4", ha="right", va="center", color="white",
        fontsize=8.5, fontweight="bold")
ax.invert_yaxis()
ax.set(title="CO$_2$ washout front (space × time)", xlabel="Time in breath (s)",
       ylabel="Volume from the mouth (mL)")
cb = fig.colorbar(im, ax=ax, pad=0.01)
cb.set_label("P$_{CO_2}$ (mmHg)", color=INK2)
cb.outline.set_visible(False)

ax = axs[1, 1]
phase_shade(ax)
ax.plot(t[:-1], ts["PCO2_mouth"][:-1], color=BLUE, label="Mouth")
ax.plot(t[:-1], ts["PCO2_trachea"][:-1], color=ORANGE, label="Trachea (Gen 0)")
ax.plot(t[:-1], ts["PCO2_gen4_5"][:-1], color=AQUA, label="Gen 4 → Gen 5")
ax.set(title="Capnograms: CO$_2$ at three boundaries", xlabel="Time in breath (s)",
       ylabel="P$_{CO_2}$ (mmHg)", ylim=(-2, 44))
ax.legend(loc="center left", fontsize=8.5)

fig.suptitle("One quiet breath through Generations 0–4 (VT 500 mL, 15 breaths/min)",
             x=0.01, ha="left", fontsize=13, fontweight="bold", color=INK)
fig.savefig(FIG / "overview.png", dpi=170)
plt.close(fig)

# ------------------------------------------------- 2. resistance by generation
fig, ax = plt.subplots(figsize=(8, 4.2), constrained_layout=True)
x = gen["gen"].to_numpy()
w = 0.36
b1 = ax.bar(x - w / 2 - 0.01, gen["R_poiseuille"] * 1e3, w, color=MUTED, label="Poiseuille (steady laminar)")
b2 = ax.bar(x + w / 2 + 0.01, gen["R_mmHg_per_Lps"] * 1e3, w, color=BLUE,
            label="This model (Pedley entrance-flow correction)")
for xi, r in zip(x, gen["R_ratio"]):
    ax.text(xi + w / 2, gen.loc[gen["gen"] == xi, "R_mmHg_per_Lps"].iloc[0] * 1e3 + 0.8,
            f"×{r:.1f}", ha="center", color=INK, fontsize=9, fontweight="bold")
ax.set_xticks(x, ["Gen 0\ntrachea", "Gen 1\nmain bronchi", "Gen 2", "Gen 3", "Gen 4"])
ax.set(title="Resistance of each generation at peak inspiratory flow")
ax.set_ylabel("Resistance (10⁻³ mmHg·s/L)")
ax.set_ylim(0, gen["R_mmHg_per_Lps"].max() * 1e3 * 1.22)
ax.legend(loc="upper left", fontsize=8.5)
ax.grid(axis="x", visible=False)
fig.savefig(FIG / "resistance_by_generation.png", dpi=170)
plt.close(fig)

# --------------------------------------------- 3. cross-section and velocity
d_cm = [1.80, 1.22, 0.83, 0.56, 0.45, 0.35, 0.28, 0.23, 0.186, 0.154, 0.130, 0.109,
        0.095, 0.082, 0.074, 0.066, 0.060]
g = np.arange(17)
area = (2.0 ** g) * np.pi * (np.array(d_cm) / 2) ** 2           # cm^2
q_pk = summ["peak_flow_mL_s"]
vel = q_pk / area                                                 # cm/s
fig, axs = plt.subplots(1, 2, figsize=(11, 3.8), constrained_layout=True)
cols = [BLUE if k <= 4 else "#c9dcf3" for k in g]
axs[0].bar(g, area, color=cols, width=0.8)
axs[0].set_yscale("log")
axs[0].set(title="Total airway cross-section by generation (Weibel)", xlabel="Generation",
           ylabel="Total area (cm², log)")
axs[1].bar(g, vel, color=cols, width=0.8)
axs[1].set(title="Mean air velocity at peak flow", xlabel="Generation", ylabel="Velocity (cm/s)")
for a in axs:
    a.set_xticks(g)
    a.grid(axis="x", visible=False)
axs[1].text(6.5, vel.max() * 0.95, "Gens 0–4: area barely changes,\nso air stays fast here",
            color=INK, fontsize=9, va="top")
fig.savefig(FIG / "area_and_velocity.png", dpi=170)
plt.close(fig)

# ----------------------------------------------------------- 4. schematic
fig, ax = plt.subplots(figsize=(11, 2.6))
ax.set_xlim(0, 100); ax.set_ylim(0, 26); ax.axis("off")
boxes = [
    (1, 15, "Mouth / nose", "air in at 760 mmHg", MUTED),
    (19, 20, "Nose → larynx", f"{summ['dead_space_nose_larynx_mL']:.0f} mL, resistance K₁\n(calibrated to Guyton)", MUTED),
    (42, 27, "Generations 0–4  (this repo)", f"{summ['dead_space_gen0_4_mL']:.1f} mL in 30 × 0.75 mL segments\nWeibel geometry, Womersley + Pedley R", BLUE),
    (72, 14, "Gens 5–16", f"{summ['dead_space_gen5_16_mL']:.1f} mL\n(boundary)", "#9cc3ef"),
    (89, 10, "Alveoli", "(boundary)", "#9cc3ef"),
]
for x0, wd, title, sub, col in boxes:
    ax.add_patch(FancyBboxPatch((x0, 6), wd, 14, boxstyle="round,pad=0.3,rounding_size=1.2",
                                fc=col, ec="none", alpha=1 if col == BLUE else 0.35))
    tc = "white" if col == BLUE else INK
    ax.text(x0 + wd / 2, 15.5, title, ha="center", va="center", fontsize=9.5, fontweight="bold", color=tc)
    ax.text(x0 + wd / 2, 10.5, sub, ha="center", va="center", fontsize=8, color=tc)
for x0 in (16.5, 39.5, 69.5, 86.5):
    ax.annotate("", xy=(x0 + 2.2, 13), xytext=(x0, 13),
                arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.4))
ax.text(50, 2.2, "Total anatomical dead space = 72 + 22.5 + 55.5 = 150 mL   •   "
        "flow reverses each breath: upwind (donor-cell) transport", ha="center", color=INK2, fontsize=8.5)
fig.savefig(FIG / "model_schematic.png", dpi=170, bbox_inches="tight")
plt.close(fig)

print("Figures written to", FIG)
