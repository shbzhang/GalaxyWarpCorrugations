import glob, os
import numpy as np
import matplotlib.pyplot as plt
from shared import *
from matplotlib.collections import LineCollection

if __name__ == '__main__':
	figscale = 0.75
	figwidth = textwidth*figscale

	fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(figwidth, figwidth*1.2))
	plt.subplots_adjust(left=0.11, right=0.95, top=0.95, bottom=0.11, wspace=0.5, hspace=0.5)

	#AX0 maxZ
	def warp_cheng(R):
		#Calculate maximum Cheng2020 warp amplitude
		R1 = 8.87
		R2 = 17.78
		R_w = 1.20
		alpha = 1.53

		h = np.zeros_like(R)
		
		mask1 = R <= R1
		mask2 = (R > R1) & (R <= R2)
		mask3 = R > R2
		
		h[mask1] = 0
		h[mask2] = (R_w / (R2 - R1)) * (R[mask2] - R1)**alpha
		h_R2 = (R_w / (R2 - R1)) * (R2 - R1)**alpha
		dh_dR_at_R2 = (R_w * alpha / (R2 - R1)) * (R2 - R1)**(alpha - 1)
		h[mask3] = h_R2 + dh_dR_at_R2 * (R[mask3] - R2)
		
		return h
	
	def warp_khanna2025(R):
		#Calculate maximum Khanna2025 warp amplitude
		R_w = 8.79
		h_w0 = 0.083
		a_w = 2.0
		phi_w = 178
		h = cal_warp_RC(R, -178-90)
		return h

	def warp_chen2019(R):
		#Calculate maximum Chen2019 warp amplitude by scanning over theta
		zm1,_ = cal_warp_Ceph(R,17.4+90)
		_,zm2 = cal_warp_Ceph(R,17.5+90)
		return zm1

	def warp_poggioYG(R, **kws):
		#Calculate maximum Poggio2025 warp amplitude by scanning over theta
		R_w = 5.5
		h_w0 = 0.012
		alpha_w = 1.9
		phi_w = 9.9

		h = np.zeros_like(R)
		mask = R > R_w
		h[mask] = h_w0 * (R[mask] - R_w)**alpha_w
		return h

	def warp_poggioCeph(R):
		#Calculate maximum Poggio2025 warp amplitude by scanning over theta
		R_w = 7.7
		h_w0 = 0.057
		alpha_w = 1.3
		phi_w = 14.0

		h = np.zeros_like(R)
		mask = R > R_w
		h[mask] = h_w0 * (R[mask] - R_w)**alpha_w
		return h

	def warp_hi(R):
		h = []
		for phi in np.arange(-90,180,0.5):
			h.append(cal_warp_HI(R, phi, scale=0.95))
		return np.max(h, axis=0)

	def warp_co(R, comp=1):
		if comp==1:
			return function_warp((R, p_1comp[-1]+90), p=p_1comp)
		else:
			h = []
			for phi in np.arange(-90,180,0.5):
				h.append(function_warp((R, phi), p=p_2comp))
			return np.max(h, axis=0)

	R = np.linspace(7, 20, 200)

	# 1. 老恒星群体 - 使用暖色调和实线
	ax[0].plot(R, warp_khanna2025(R), color='#8B4513', linewidth=2.5, label='Red clump (Khanna+2025)', 
			 alpha=0.8, linestyle='-')
	ax[0].plot(R, warp_cheng(R), color='#D2691E', linewidth=2.5, label='Gaia DR2/APOGEE (Cheng+2020)', 
			 alpha=0.8, linestyle='-')

	# 2. 年轻恒星群体 
	ax[0].plot(R, warp_chen2019(R), color='#1E90FF', linewidth=2.5, label='Cepheids (Chen+2019)', 
			 alpha=0.8, linestyle='-')

	ax[0].plot(R, warp_poggioYG(R), color='#4682B4', linewidth=2.5, label='Young Giants (Poggio+2025)', 
			 alpha=0.8, linestyle='--')

	ax[0].plot(R, warp_poggioCeph(R), color='b', linewidth=2.5, label='Cepheids (Poggio+2025)', 
			 alpha=0.8, linestyle='--')

	# 3. 原子气体 - 使用绿色系和点划线
	ax[0].plot(R, warp_hi(R), color='#2E8B57', linewidth=2.5, label='HI (Levine+2006)', 
			 alpha=0.8, linestyle='-.')

	# 4. 分子气体 - 特别突出显示（我们的结果）
	ax[0].plot(R, warp_co(R), color='r', linewidth=2, label='CO m=1 - This work', 
			 alpha=0.9, linestyle='-', zorder=100)

	# realizations
	modelFiles = glob.glob('oneCompExc/steps_errD_*pc_*mass.npy')
	steps = [np.load(f)[-100:] for f in modelFiles]
	steps = np.vstack(steps)
	zmodels = []
	for s in steps:
		a,R0,b,_ = s
		z = np.zeros_like(R)
		idx = R>R0
		z[idx] = a*(R[idx]-R0)**b
		zmodels.append(np.array((R,z)).T)
	lc = LineCollection(zmodels, color='#ffcccc', alpha=0.2, lw=1, zorder=0)
	lc.set_rasterized(True)
	ax[0].add_collection(lc)


	# 1. 添加单一浅灰色背景填充
	ax[0].fill_between(R, 0, 2.5, alpha=0.05, color='#f5f5f5', edgecolor='none')

	# 设置标签
	#ax[0].set_xticklabels([' ']*4)
	#ax[0].set_xlabel('R (kpc)')
	ax[0].set_ylabel('Maximum Z (kpc)')
	ax[0].grid(True, alpha=0.3)
	ax[0].legend(fontsize=16, loc='upper left', framealpha=0.95, ncol=2)
	ax[0].set_ylim(-0.1, 2.5)
	ax[0].set_xlim(9.5, 16.6)  # 调整x轴范围

	# 调整群体标注框位置
	ax[0].text(14, 1.3, 'OLD \nPOPULATIONS', fontsize=18, color='#8B4513', 
			 ha='center', va='center', weight='bold',
			 bbox=dict(boxstyle="round,pad=0.3", facecolor='#FAEBD7', alpha=0.8))

	# 年轻恒星群体标注 
	ax[0].text(14, 0.2, 'YOUNG \nPOPULATIONS', fontsize=18, color='#1E90FF', 
			 ha='center', va='center', weight='bold',
			 bbox=dict(boxstyle="round,pad=0.3", facecolor='#F0F8FF', alpha=0.8))


	# 打印分组统计信息，特别突出CO结果
	mR = max(R)
	print("=" * 100)
	print("GALACTIC WARP AMPLITUDES BY POPULATION TYPE")
	print("=" * 100)
	print("\n--- OLD STELLAR POPULATIONS ---")
	print(f"Red clump (Khanna+2025):      Max = {warp_khanna2025(mR):.3f} kpc at R = {mR:.1f} kpc")
	print(f"Gaia DR2/APOGEE (Cheng+2020): Max = {warp_cheng(mR):.3f} kpc at R = {mR:.1f} kpc")

	print("\n--- YOUNG STELLAR POPULATIONS ---")
	print(f"Young Giants (Poggio+2025):   Max = {warp_poggioYG(mR):.3f} kpc at R = {mR:.1f} kpc")
	print(f"Cepheids (Poggio+2025):       Max = {warp_poggioCeph(mR):.3f} kpc at R = {mR:.1f} kpc")
	print(f"Cepheids (Chen+2019):         Max = {warp_chen2019(mR):.3f} kpc at R = {mR:.1f} kpc")

	print("\n--- GASEOUS COMPONENTS ---")
	print(f"HI (Levine+2006):             Max = {max(warp_hi([mR])):.3f} kpc at R = {mR:.1f} kpc")
	print(warp_co([mR]))
	print(f"CO m=1 - This work:           Max = {max(warp_co([mR])):.3f} kpc at R = {mR:.1f} kpc")
	print(f"CO m=1,2 - This work:         Max = {max(warp_co([mR], comp=2)):.3f} kpc at R = {mR:.1f} kpc ★")



	#AX1: node
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


	# Use SAME styles as amplitude plot
	ax[1].plot(cab["R"], cab_phi,  color="#555555", lw=2.5, ls="-", alpha=0.8, label="Cabrera+2024 (Cepheids)")
	ax[1].plot(dehn["R"], dehn_phi, color="g", lw=2.5, ls="-",  alpha=0.8, label="Dehnen+2023 (Cepheids)")
	ax[1].plot(chen["R"], chen_phi, color="#1E90FF", lw=2.5, ls="--", alpha=0.8, label="Chen+2019 (Cepheids)")

	ax[1].plot(R_grid, pog_phi, color="b", lw=2.5, ls="--", alpha=0.9, label="Poggio+2025 (Cepheids)")

	# CO points: keep consistent with your “This work” highlight (open markers)
	ax[1].scatter(R_co, phi_co_plot, s=70, facecolors="grey", edgecolors="r",
	           linewidths=2.0, zorder=10, label="CO (this work)")
	ax[1].errorbar(R_co, phi_co_plot, yerr=[phi_co_error_lo, phi_co_error_hi], linestyle='none', ecolor='r', elinewidth=2, capsize=4, capthick=2)

	# Axis labels consistent with amplitude plot
	ax[1].set_xlabel("R (kpc)")
	ax[1].set_ylabel(r"$\phi_{\rm LON}$ (deg)")

	ax[1].grid(True, alpha=0.3)
	ax[1].set_xlim(9.5, 16.6)       # match amplitude’s x-range if you want
	ax[1].set_ylim(-40, 40)     # tune if needed

	# Legend style: copy your amplitude plot style
	ax[1].legend(fontsize=16, loc="upper left", framealpha=0.95, ncol=2)

	ax[0].text(-0.07, 0.96, 'a', color='black', font=subfigureIndexFont, transform=ax[0].transAxes)
	ax[1].text(-0.07, 0.96, 'b', color='black', font=subfigureIndexFont, transform=ax[1].transAxes)

	plt.tight_layout()

	plt.savefig("fig/multi_warp.pdf", dpi=300, bbox_inches="tight")
	plt.savefig("fig/multi_warp.png", dpi=300, bbox_inches="tight")
	plt.show()






