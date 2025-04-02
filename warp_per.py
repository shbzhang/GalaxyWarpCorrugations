### plot Az-Z and Az-(Z-warp model) of MCs in Perseus arms

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sp_interp
from scipy.stats import norm
import scipy.optimize as sp_opt
from shared import *

if __name__ == '__main__':
	col_mc = '#6fa01f'

	data = np.loadtxt('per_para.txt',comments='#')
	ll = data[:,0]
	bb = data[:,1]
	vv = data[:,2]
	rr = data[:,11]
	zz = data[:,14]
	az = data[:,19]
	mass = data[:,7]/data[:,8]
	print('Az range:', az.min(), az.max())

	### calculate arm
	PHI = np.linspace(az.min()-3, az.max()+3, 300)
	R = function_arm(PHI, best_per)

	fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=[9,8])
	plt.subplots_adjust(wspace=0.1,hspace=0.1)
	plt.rcParams['xtick.direction'] = 'in'

	### plot cloud
	scatter_size,scatter_color = mass_distribute(mass)
	ax[0].scatter(az, zz, s=scatter_size, c=col_mc, edgecolor='none', alpha=0.3, label='Perseus')

	### plot binned average
	azcen, azrms, zcen, zrms, vcen, vrms, rcen, rrms = cal_zcen_zrms(az, zz, vv, rr, weights=mass, binsize=3.5, nbin=60)
	ax[0].plot(azcen, zcen, '-', c=darker_hex(col_mc), lw=2, alpha=0.8)

	### plot warp models
	zh = cal_warp(R, PHI)
	ax[0].plot(PHI, zh, sty_hi, c=col_hi,lw=2, label='HI')

	z1, z2, z4, z5 =  cal_warpc(R, PHI)
	ax[0].plot(PHI, z1, sty_ceph, c=col_ceph, lw=2, label='Cepheids') #Chen b=1')
	ax[0].plot(PHI, z4, sty_co1, c=col_co,lw=2, label='CO 1comp')
	ax[0].plot(PHI, z5, sty_co2, c=col_co,lw=2, label='CO 2comp')
	#z1=np.array(per_zcen)+np.array(per_zrms)*2.355/2.
	#z2=np.array(per_zcen)-np.array(per_zrms)*2.355/2.
	ax[0].legend()
	ax[0].set_yticks(np.arange(-3, 3, 0.2))
	ax[0].minorticks_on()
	ax[0].set_xlim(-30, 100)
	ax[0].set_ylim(-0.4, 0.6)

	# upper tick
	import scipy.interpolate as sp_interp
	phiAxis = np.linspace(-30, 100, 200)
	lAxis = []
	for phi in phiAxis:
		lAxis.append(cal_arm_length(-30, phi, p=best_per))
	### interpolate arm length over phi
	f = sp_interp.interp1d(lAxis, phiAxis, fill_value='extrapolate')

	lTicks = np.arange(0, 25, 5)
	phiTicks = f(lTicks)
	upper = ax[0].twiny()
	upper.set_xticks(phiTicks)
	lTickLabels = lTicks.astype(str)
	lTickLabels[-1] += ' kpc'
	upper.set_xticklabels(lTickLabels) # add kpc at the end
	upper.set_xlim(-30, 100)

	### plot residual
	#zz_res = zz - function_warp((rr, az))
	zz_res = zz - function_warp( (function_arm(az, best_per), az) )
	ax[1].scatter(az, zz_res, s=scatter_size, c=col_mc, edgecolor='none', alpha=0.3, label='Residual')

	### plot binned average
	azcen, azrms, zcen_res, zrms_res, vcen, vrms, rcen, rrms = cal_zcen_zrms(az, zz_res, vv, rr, weights=mass, binsize=3.5, nbin=60)
	ax[1].plot(azcen, zcen_res, '-', c=darker_hex(col_mc), lw=2, alpha=0.8)

	### plot H line
	ax[1].plot([-30, 150], [0, 0], '--', color=col_co)
	ax[1].set_yticks(np.arange(-3, 3, 0.2))
	ax[1].minorticks_on()
	ax[1].set_ylim([-0.4,0.4])
	#ax[1].yticks([-500,0,500,1000,1500],['-500','0','500','1000','1500'],fontsize=14,fontweight=1.8)
	#ax[1].xticks([-25,0,25,50,75,100,125,150],['-25','0','25','50','75','100','125','150'],fontsize=14,fontweight=1.8)
	ax[1].grid(True, ls='--', alpha=0.4)
	ax[1].set_xlabel('Galactocentric Azimuth (degree)', fontsize=18, fontweight=1.8)
	ax[1].set_ylabel('Z (kpc)', fontsize=18, fontweight=1.8)
	plt.savefig('fig/per_warp_corrugation.png',format='png',bbox_inches='tight', dpi=400) 
	plt.show()
