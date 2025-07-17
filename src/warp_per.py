### plot Az-Z and Az-(Z-warp model) of MCs in Perseus arms

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import scipy.interpolate as sp_interp
from scipy.stats import norm
import scipy.optimize as sp_opt
from shared import *

if __name__ == '__main__':
	figscale = 0.45
	figwidth = textwidth*figscale
	#col_mc = '#4169e1'#'#6fa01f'

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
	azmax = az.max()+5
	PHI = np.linspace(azmin, azmax, 300)
	R = function_arm(PHI, best_per)
	print('Az range:', az.min(), az.max())
	print('Az plot range:', azmin, azmax)

	fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(figwidth, figwidth*0.95))
	plt.subplots_adjust(left=0.13, right=0.98, wspace=0.1, hspace=0.06)
	plt.rcParams['xtick.direction'] = 'in'

	### plot cloud
	scatter_size,scatter_color = mass_distribute(mass)
	ax[0].scatter(az, zz, s=scatter_size*2, **arm_kws_mc)
	ax[0].text(0.02, 0.95, 'Perseus', transform=ax[0].transAxes, **arm_kws_text)

	### plot binned average
	azcen, azrms, zcen, zrms = bin_data(az, zz, mass, grid=np.arange(160, -50, -2), width=4, method='gauss')
	#azcen, azrms, zcen, zrms, vcen, vrms, rcen, rrms = cal_zcen_zrms(az, zz, vv, rr, weights=mass, binsize=4)
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

	ax[0].legend(frameon=False, borderpad=0.2, labelspacing=0.1)
	ax[0].set_yticks(np.arange(-3, 3, 0.2))
	ax[0].set_xlim(azmin, azmax+20)
	ax[0].set_ylim(-0.4, 0.6)
	ax[0].set_ylabel('Z (kpc)')

	# upper tick
	import scipy.interpolate as sp_interp
	phiAxis = np.linspace(azmin, azmax+20, 200)
	lenAxis = []
	for phi in phiAxis:
		lenAxis.append(cal_arm_length(azmin, phi, p=best_per))
	### interpolate arm length over phi
	f = sp_interp.interp1d(lenAxis, phiAxis, fill_value='extrapolate')

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
	upper.set_xticklabels(lTickLabels) # add kpc at the end
	upper.set_xlim(azmin, azmax+20)

	### plot residual
	#zz_res = zz - function_warp((rr, az))
	zz_res = zz - function_warp( (function_arm(az, best_per), az), p=p_1comp)
	ax[1].scatter(az, zz_res, s=scatter_size*2, **arm_kws_mc)

	### plot binned average
	azcen, azrms, zcen_res, zrms_res = bin_data(az, zz_res, mass, grid=np.arange(160, -50, -2), width=4, method='gauss')
	#azcen, azrms, zcen_res, zrms_res, vcen, vrms, rcen, rrms = cal_zcen_zrms(az, zz_res, vv, rr, weights=mass, binsize=3.5)
	ax[1].plot(azcen, zcen_res, '-', **arm_kws_bin)

	### connect broken
	idx = np.isfinite(zcen_res)
	ax[1].plot(azcen[idx], zcen_res[idx], '--', **arm_kws_bin)

	### plot H line
	arm_kws_co1['zorder']=0
	ax[1].plot([PHI[0], PHI[-1]], [0, 0], **arm_kws_co1)

	ax[1].set_yticks(np.arange(-3, 3, 0.2))
	ax[1].set_ylim([-0.4,0.4])
	#ax[1].yticks([-500,0,500,1000,1500],['-500','0','500','1000','1500'],fontsize=14,fontweight=1.8)
	#ax[1].xticks([-25,0,25,50,75,100,125,150],['-25','0','25','50','75','100','125','150'],fontsize=14,fontweight=1.8)
	ax[1].grid(True, ls='--', alpha=0.4)
	ax[1].set_xlabel('Galactocentric Azimuth (deg)')
	ax[1].set_ylabel('$\mathbf{\Delta}$Z (kpc)')

	insert_arm_plot(ax[0], [0.74, 0.15, 0.36, 0.36], boldArm='per', boldRange=[azmin,azmax])

	if 0:
		ax[1].set_ylim([-0.4, 0.65])
		### insert a power spectrum
		#powerAxes = ax[1].inset_axes([0.21, 0.12, 0.21, 0.21])
		powerAxes = ax[1].inset_axes([0.77, 0.76, 0.21, 0.21])
		from astropy.timeseries import LombScargle
		f = sp_interp.interp1d(phiAxis, lenAxis, fill_value='extrapolate')
		length = f(az)
		ls = LombScargle(length, zz_res, dy=1/np.sqrt(mass), normalization='standard')
		frequency, power = ls.autopower(minimum_frequency=1/50, maximum_frequency=1/1)
		powerAxes.plot(frequency, power, color='#4169e1')

		# Find the dominant period
		peak = np.max(power)
		bestFreq = frequency[np.argmax(power)]
		powerAxes.plot([bestFreq, bestFreq, bestFreq+0.1], [peak+0.01, peak+0.03, peak+0.03], 'k', lw=0.5, alpha=0.8)
		powerAxes.text(bestFreq+0.11, peak+0.03, '%.1f kpc' % (1/bestFreq), ha='left', va='center')
		# find the second period
		#powerAxes.plot([frequency[20],frequency[20]], [0,1])
		peak = np.max(power[20:])
		bestFreq = frequency[np.argmax(power[20:])+20]
		powerAxes.plot([bestFreq, bestFreq, bestFreq+0.1], [peak+0.01, peak+0.03, peak+0.03], 'k', lw=0.5, alpha=0.8)
		powerAxes.text(bestFreq+0.11, peak+0.03, '%.1f kpc' % (1/bestFreq), ha='left', va='center')

		#powerAxes.set_xscale('log')
		powerAxes.set_xlim(0, 1)
		powerAxes.set_ylim(-0.01, np.max(power)*1.6)
		powerAxes.tick_params(top=False, left=True, direction='in', pad=3)
		powerAxes.minorticks_on()
		powerAxes.set_xlabel('frequency (kpc$^{-1}$)', labelpad=-2)
		powerAxes.set_ylabel('Norm\nPower', labelpad=0)
		powerAxes.spines['top'].set_visible(False)
		#powerAxes.spines['left'].set_visible(False)
		powerAxes.spines['right'].set_visible(False)
		#powerAxes.yaxis.set_visible(False)
		powerAxes.patch.set_alpha(0.0)


	plt.savefig('fig/per_warp_corrugation.%s' % (mpl.rcParams['savefig.format']), bbox_inches='tight')
	plt.savefig('fig/per_warp_corrugation.png', bbox_inches='tight')
	plt.show()
