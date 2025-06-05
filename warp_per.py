### plot Az-Z and Az-(Z-warp model) of MCs in Perseus arms

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import scipy.interpolate as sp_interp
from scipy.stats import norm
import scipy.optimize as sp_opt
from shared import *

if __name__ == '__main__':
	col_mc = '#4169e1'#'#6fa01f'

	data = np.loadtxt('per_para.txt',comments='#')
	ll = data[:,0]
	bb = data[:,1]
	vv = data[:,2]
	rr = data[:,11]
	zz = data[:,14]
	az = data[:,19]
	mass = data[:,7]/data[:,8]

	### calculate arm
	azmin = az.min()-3
	azmax = az.max()+3
	PHI = np.linspace(azmin, azmax, 300)
	R = function_arm(PHI, best_per)
	print('Az range:', az.min(), az.max())
	print('Az plot range:', azmin, azmax)

	fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=[9,8])
	plt.subplots_adjust(wspace=0.1,hspace=0.1)
	plt.rcParams['xtick.direction'] = 'in'

	### plot cloud
	scatter_size,scatter_color = mass_distribute(mass)
	ax[0].scatter(az, zz, s=scatter_size*2, **arm_kws_mc)
	ax[0].text(0.02, 0.95, 'Perseus', transform=ax[0].transAxes, **arm_kws_text)

	### plot binned average
	azcen, azrms, zcen, zrms, vcen, vrms, rcen, rrms = cal_zcen_zrms(az, zz, vv, rr, weights=mass, binsize=4, nbin=60)
	ax[0].plot(azcen, zcen, '-', **arm_kws_bin)

	# connect broken with dashed
	idx = np.isfinite(zcen)
	ax[0].plot(azcen[idx], zcen[idx], '--', **arm_kws_bin)

	### plot warp models
	zh = cal_warp(R, PHI)
	ax[0].plot(PHI, zh, **arm_kws_hi)

	z1, z2, z4, z5 =  cal_warpc(R, PHI)
	ax[0].plot(PHI, z1, **arm_kws_ceph) #Chen b=1')
	ax[0].plot(PHI, z4, **arm_kws_co1)
	ax[0].plot(PHI, z5, **arm_kws_co2)

	ax[0].legend()
	ax[0].set_yticks(np.arange(-3, 3, 0.2))
	ax[0].minorticks_on()
	ax[0].set_xlim(azmin, azmax+20)
	ax[0].set_ylim(-0.4, 0.6)

	ax[0].tick_params(right=True, direction='in', labelsize=12)#, labelsize=1000/self.dpi)
	ax[0].tick_params(which='minor', right=True, direction='in')#, labelsize=1000/self.dpi)

	# upper tick
	import scipy.interpolate as sp_interp
	phiAxis = np.linspace(azmin, azmax+20, 200)
	lAxis = []
	for phi in phiAxis:
		lAxis.append(cal_arm_length(azmin, phi, p=best_per))
	### interpolate arm length over phi
	f = sp_interp.interp1d(lAxis, phiAxis, fill_value='extrapolate')

	upper = ax[0].twiny()
	# minor
	lTicks = np.arange(0, 45, 1)
	phiTicks = f(lTicks)
	upper.set_xticks(phiTicks, minor=True)
	# major
	lTicks = np.arange(0, 25, 5)
	phiTicks = f(lTicks)
	upper.set_xticks(phiTicks)
	lTickLabels = lTicks.astype(str)
	lTickLabels[-1] += ' kpc'
	upper.set_xticklabels(lTickLabels, fontsize=12) # add kpc at the end
	upper.set_xlim(azmin, azmax+20)

	### plot residual
	#zz_res = zz - function_warp((rr, az))
	zz_res = zz - function_warp( (function_arm(az, best_per), az) )
	ax[1].scatter(az, zz_res, s=scatter_size*2, **arm_kws_mc)

	### plot binned average
	azcen, azrms, zcen_res, zrms_res, vcen, vrms, rcen, rrms = cal_zcen_zrms(az, zz_res, vv, rr, weights=mass, binsize=3.5, nbin=60)
	ax[1].plot(azcen, zcen_res, '-', **arm_kws_bin)

	### connect broken
	idx = np.isfinite(zcen_res)
	ax[1].plot(azcen[idx], zcen_res[idx], '--', **arm_kws_bin)

	### plot H line
	arm_kws_co1['zorder']=0
	ax[1].plot([PHI[0], PHI[-1]], [0, 0], **arm_kws_co1)

	ax[1].set_yticks(np.arange(-3, 3, 0.2))
	ax[1].minorticks_on()
	ax[1].set_ylim([-0.4,0.4])
	#ax[1].yticks([-500,0,500,1000,1500],['-500','0','500','1000','1500'],fontsize=14,fontweight=1.8)
	#ax[1].xticks([-25,0,25,50,75,100,125,150],['-25','0','25','50','75','100','125','150'],fontsize=14,fontweight=1.8)
	ax[1].grid(True, ls='--', alpha=0.4)
	ax[1].set_xlabel('Galactocentric Azimuth (deg)', fontsize=15, fontweight='bold')
	ax[1].set_ylabel('Z (kpc)', fontsize=15, fontweight='bold')

	ax[1].tick_params(top=True, right=True, direction='in', labelsize=12)#, labelsize=1000/self.dpi)
	ax[1].tick_params(which='minor', top=True, right=True, direction='in')#, labelsize=1000/self.dpi)

	insert_arm_plot(ax[0], [0.74, 0.15, 0.36, 0.36], boldArm='per', boldRange=[azmin,azmax])


	plt.savefig('fig/per_warp_corrugation.png',format='png',bbox_inches='tight', dpi=400) 
	plt.show()
