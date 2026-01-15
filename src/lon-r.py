import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- reuse your THETA_DEG / STYLE / wrap funcs ---
def wrap_to_pm180(phi_deg):
    phi = np.asarray(phi_deg, dtype=float)
    return (phi + 180.0) % 360.0 - 180.0

def cabrera_to_CO_frame(phi_cab_deg):
    return 180.0 - np.asarray(phi_cab_deg, dtype=float)

def phi_twisted_cepheids(R):
    phi0, Rt, beta = 0.9, 12.6, 8.0
    R = np.asarray(R, dtype=float)
    return np.where(R <= Rt, phi0, phi0 + beta*(R - Rt))

# ---- load the 3 curves (Cabrera frame: Sun=180; columns: R, LON) ----
cab  = pd.read_csv("LON_R_Cabrera_Gadea_2024.csv")
chen = pd.read_csv("LON_R_Chen_2019.csv")
dehn = pd.read_csv("LON_R_Dehnen_2023.csv")

# Convert to CO frame (Sun=0) then wrap for compact plot
cab_phi  = wrap_to_pm180(cabrera_to_CO_frame(cab["LON"].to_numpy()))
chen_phi = wrap_to_pm180(cabrera_to_CO_frame(chen["LON"].to_numpy()))
dehn_phi = wrap_to_pm180(cabrera_to_CO_frame(dehn["LON"].to_numpy()))

# Poggio twisted LON (already CO frame)
R_grid = np.linspace(9.5, 16.5, 400)
pog_phi = wrap_to_pm180(phi_twisted_cepheids(R_grid))

# Your CO LON points (already CO frame)
R_co = np.array([10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5])
phi_co = np.array([10.1138, -8.5029, -8.4242, -11.5317, 2.2867, 3.0608, -2.1427])
phi_co_plot = wrap_to_pm180(phi_co)

# ---- Figure: match amplitude plot sizing & legend style ----
fig, ax = plt.subplots(figsize=(6.5, 4.5))

# Use SAME styles as amplitude plot
ax.plot(cab["R"], cab_phi,  color="#555555", lw=2.5, ls="-", alpha=0.8, label="Cabrera+2024 (Cepheids)")
ax.plot(dehn["R"], dehn_phi, color="g", lw=2.5, ls="-",  alpha=0.8, label="Dehnen+2023 (Cepheids)")
ax.plot(chen["R"], chen_phi, color="#1E90FF", lw=2.5, ls="--", alpha=0.8, label="Chen+2019 (Cepheids)")

ax.plot(R_grid, pog_phi, color="b", lw=2.5, ls="--", alpha=0.9, label="Poggio+2025 (Cepheids)")

# CO points: keep consistent with your “This work” highlight (open markers)
ax.scatter(R_co, phi_co_plot, s=70, facecolors="grey", edgecolors="r",
           linewidths=2.0, zorder=10, label="CO (this work)")

# Axis labels consistent with amplitude plot
ax.set_xlabel("R (kpc)", fontsize=11, fontweight="bold")
ax.set_ylabel(r"$\phi_{\rm LON}$ (deg)", fontsize=11, fontweight="bold")

ax.grid(True, alpha=0.3)
ax.set_xlim(9.5, 16.6)       # match amplitude’s x-range if you want
ax.set_ylim(-40, 40)     # tune if needed

# Legend style: copy your amplitude plot style
ax.legend(fontsize=9, loc="upper left", framealpha=0.95, ncol=2)
ax.text(-0.08, 1., "b", fontsize=11,transform=ax.transAxes, va="top", ha="right",
         fontweight="bold")
plt.tight_layout()
plt.savefig("multi_warp_LON.png", dpi=300, bbox_inches="tight")
plt.show()

