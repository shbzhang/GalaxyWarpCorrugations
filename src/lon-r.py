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

#v = [+13.76, +6.405, -5.883, -7.884, -6.507, -4.916, -3.726, -1.362, +1.862]
#u = [+1.140, +0.610, +0.573, +0.670, +0.690, +0.713, +0.861, +1.716, +0.269]
#l = [-1.200, -0.616, -0.597, -0.695, -0.720, -0.776, -0.936, -1.763, -0.273]
v= [13.952136186266202, 6.970468027824097, -6.158314334143681, -7.358465264798522, -6.418554401960492, -1.7002470794191076, 0.7600196965323002, 0.8879777423578942, 4.545418596118302] 
l= [-1.176560001224772, -0.601040896573814, -0.5408510069066195, -0.5164368558979816, -5.479260919555037, -0.6056238596871395, -0.9734774064661539, -2.8055894181621457, -2.8515504083368963] 
u= [10.394607812456364, 5.1433493345113925, 1.2484348211487681, 1.150509999088195, 0.7559567782123796, 3.239083445674101, 2.1634062547989457, 1.3297327950252762, 1.0299301264634655]
# Your CO LON points (already CO frame)
R_co = np.array([9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5])[1:-1]
phi_co = np.array(v)[1:-1]#[23.1469, 10.1138, -8.5029, -8.4242, -11.5317, 2.2867, 3.0608, -2.1427, 1.8616])[1:-1]
phi_co_plot = wrap_to_pm180(phi_co)
phi_co_error_hi = np.array(u)[1:-1]#[ +6.2352, +3.6379, +4.4849, +4.0715, +4.3975,  +5.0159, +5.5255, +5.3171, +5.9002])[1:-1]
phi_co_error_lo = -np.array(l)[1:-1]#[-11.3540, -3.9882, -4.6857, -6.7652, -4.0458, -10.9955, -9.5311, -5.7419, -5.8337])[1:-1]


# ---- Figure: match amplitude plot sizing & legend style ----
fig, ax = plt.subplots(figsize=(6.5, 4.5))
plt.subplots_adjust(left=0.11, right=0.97, top=0.95, bottom=0.1)

# Use SAME styles as amplitude plot
ax.plot(cab["R"], cab_phi,  color="#555555", lw=2.5, ls="-", alpha=0.8, label="Cabrera+2024 (Cepheids)")
ax.plot(dehn["R"], dehn_phi, color="g", lw=2.5, ls="-",  alpha=0.8, label="Dehnen+2023 (Cepheids)")
ax.plot(chen["R"], chen_phi, color="#1E90FF", lw=2.5, ls="--", alpha=0.8, label="Chen+2019 (Cepheids)")

ax.plot(R_grid, pog_phi, color="b", lw=2.5, ls="--", alpha=0.9, label="Poggio+2025 (Cepheids)")

# CO points: keep consistent with your “This work” highlight (open markers)
ax.scatter(R_co, phi_co_plot, s=70, facecolors="grey", edgecolors="r",
           linewidths=2.0, zorder=10, label="CO (this work)")
ax.errorbar(R_co, phi_co_plot, yerr=[phi_co_error_lo, phi_co_error_hi], linestyle='none', ecolor='r', elinewidth=2, capsize=4, capthick=2)

# Axis labels consistent with amplitude plot
ax.set_xlabel("R (kpc)", fontsize=11, fontweight="bold")
ax.set_ylabel(r"$\phi_{\rm LON}$ (deg)", fontsize=11, fontweight="bold")

ax.grid(True, alpha=0.3)
ax.set_xlim(9.5, 16.6)       # match amplitude’s x-range if you want
ax.set_ylim(-40, 40)     # tune if needed

# Legend style: copy your amplitude plot style
ax.legend(fontsize=9, loc="upper left", framealpha=0.95, ncol=2)
#ax.text(-0.08, 1., "b", fontsize=11,transform=ax.transAxes, va="top", ha="right",
#         fontweight="bold")
ax.text(-0.1, 0.96, 'b', fontsize=18, transform=ax.transAxes)

plt.tight_layout()
plt.savefig("fig/multi_warp_LON.pdf", dpi=300, bbox_inches="tight")
plt.savefig("fig/multi_warp_LON.png", dpi=300, bbox_inches="tight")
plt.show()

