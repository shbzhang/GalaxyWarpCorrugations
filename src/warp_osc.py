### plot Az-Z and Az-(Z-warp model) of MCs in OSC arms

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

	data = np.loadtxt('osc2_para.txt',comments='#')
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
	R = function_arm(PHI, best_osc)
	print('Az range:', az.min(), az.max())
	print('Az plot range:', azmin, azmax)


	fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(figwidth, figwidth*0.95))
	plt.subplots_adjust(left=0.13, right=0.98, wspace=0.1, hspace=0.06)
	plt.rcParams['xtick.direction'] = 'in'



	### plot cloud
	scatter_size,scatter_color = mass_distribute(mass)
	ax[0].scatter(az, zz, s=scatter_size*2, **arm_kws_mc)
	ax[0].text(0.02, 0.95, 'OSC', transform=ax[0].transAxes, **arm_kws_text)

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

	ax[0].legend(loc = (0.52, 0.02), frameon=False, borderpad=0.2, labelspacing=0.1)
	ax[0].set_yticks(np.arange(-3, 3, 0.5))
	ax[0].set_xlim(azmin, azmax)
	ax[0].set_ylim(-0.7, 1.5)
	ax[0].set_ylabel('Z (kpc)')

	# upper tick in kpc
	import scipy.interpolate as sp_interp
	phiAxis = np.linspace(azmin, azmax, 200)
	lenAxis = []
	for phi in phiAxis:
		lenAxis.append(cal_arm_length(azmin, phi, p=best_osc))
	### interpolate arm length over phi
	f = sp_interp.interp1d(lenAxis, phiAxis, fill_value='extrapolate')

	upper = ax[0].twiny()
	# minor
	lTicks = np.arange(0, 70, 2)
	phiTicks = f(lTicks)
	upper.set_xticks(phiTicks, minor=True)
	# major
	lTicks = np.arange(0, 60, 10)
	phiTicks = f(lTicks)
	upper.set_xticks(phiTicks)
	lTickLabels = lTicks.astype(str)
	lTickLabels[-1] += ' kpc'
	upper.set_xticklabels(lTickLabels) # add kpc at the end
	upper.set_xlim(azmin, azmax)



	### plot cloud
	#zz_res = zz - function_warp((rr, az))
	zz_res = zz - function_warp( (function_arm(az, best_osc), az), p=p_1comp)
	ax[1].scatter(az, zz_res, s=scatter_size*2, **arm_kws_mc)

	### plot binned average
	azcen, azrms, zcen_res, zrms_res = bin_data(az, zz_res, mass, grid=np.arange(160, -50, -2), width=4, method='gauss')
	#azcen, azrms, zcen_res, zrms_res, vcen, vrms, rcen, rrms = cal_zcen_zrms(az, zz_res, vv, rr, weights=mass, binsize=4)
	ax[1].plot(azcen, zcen_res, '-', **arm_kws_bin)

	### connect broken
	idx = np.isfinite(zcen_res)
	ax[1].plot(azcen[idx], zcen_res[idx], '--', **arm_kws_bin)

	### plot H line
	arm_kws_co1['zorder']=0
	ax[1].plot([PHI[0], PHI[-1]], [0, 0], **arm_kws_co1)

	ax[1].set_yticks(np.arange(-3, 3, 0.5))
	ax[1].set_ylim([-0.8, 0.8])
	#ax[1].yticks([-500,0,500,1000,1500],['-500','0','500','1000','1500'],fontsize=14,fontweight=1.8)
	#ax[1].xticks([-25,0,25,50,75,100,125,150],['-25','0','25','50','75','100','125','150'],fontsize=14,fontweight=1.8)
	ax[1].grid(True, ls='--', alpha=0.4)
	ax[1].set_xlabel('Galactocentric Azimuth (deg)')
	ax[1].set_ylabel('$\mathbf{\Delta}$Z (kpc)')

	insert_arm_plot(ax[0], [0.25, 0.02, 0.36, 0.36], boldArm='osc', boldRange=[azmin,azmax])

	if 0:
		ax[1].set_ylim([-0.8, 1.2])
		### insert a power spectrum
		#powerAxes = ax[1].inset_axes([0.21, 0.12, 0.21, 0.21])
		powerAxes = ax[1].inset_axes([0.77, 0.76, 0.21, 0.21])
		from astropy.timeseries import LombScargle
		f = sp_interp.interp1d(phiAxis, lenAxis, fill_value='extrapolate')
		length = f(az)
		ls = LombScargle(length, zz_res, dy=1/np.sqrt(mass), normalization='standard')
		frequency, power = ls.autopower(minimum_frequency=1/100, maximum_frequency=1/2)
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
		powerAxes.set_xlim(0, 0.5)
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


	plt.savefig('fig/osc_warp_corrugation.%s' % (mpl.rcParams['savefig.format']), bbox_inches='tight')
	plt.savefig('fig/osc_warp_corrugation.png', bbox_inches='tight')
	plt.show()

'''
def function(x,a,b):

        y =  a * np.exp(-1*x/b)

        return(y)

def mass_distribute(mass):

    N = len(mass)
    mass_dist = np.zeros(N)
    colors = []
    for i in range(N):
    	if mass[i] < 100:
    		mass_dist[i] = 0.01
    		colors.append('purple')
    	elif mass[i] >= 100 and mass[i] < 1000:
    		mass_dist[i] = 0.1
    		colors.append('blue')
    	elif mass[i] >= 1000 and mass[i] < 10000:
    		mass_dist[i] = 1
    		colors.append('green')
    	elif mass[i] >= 10000 and mass[i] < 100000:
    		mass_dist[i] = 10
    		colors.append('red')
    	else:
    		mass_dist[i] = 20
    		colors.append('cyan')
    return mass_dist * 10,colors
def weighted_avg_and_std(values, weights):
	average = np.average(values, weights=weights)
	# Fast and numerically precise:
	variance = np.average((values-average)**2, weights=weights)
	return (average, np.sqrt(variance))
def cal_zcen_zrms(az, zz, vv, rr, mass, bi):
    zcen, zrms, vcen, vrms, rcen, rrms, thea, z_err, v_err, r_err, zcen_err, vcen_err, rcen_err, xcen, ycen = [], [], [], [], [], [], [], [], [], [], [], [], [], [], []
    for i in range(60):
        thea1 = 160 - i * bi
        thea2 = 160 - (i + 1) * bi
        ind = np.where((az > thea2) & (az < thea1))[0]  # Extract the array of indices
        if len(ind) > 5:
            temp_zcen, temp_zrms = weighted_avg_and_std(zz[ind], mass[ind])
            temp_vcen, temp_vrms = weighted_avg_and_std(vv[ind], mass[ind])
            temp_rcen, temp_rrms = weighted_avg_and_std(rr[ind], mass[ind])
            zcen.append(temp_zcen)
            zrms.append(temp_zrms)
            vcen.append(temp_vcen)
            vrms.append(temp_vrms)
            rcen.append(temp_rcen)
            rrms.append(temp_rrms)
            thea.append((thea2 + bi / 2.))
            z_err.append(temp_zrms / np.sqrt(len(ind)))
            v_err.append(temp_vrms / np.sqrt(len(ind)))
            r_err.append(temp_rrms / np.sqrt(len(ind)))
            zcen_err.append(np.abs(temp_zcen) / np.sqrt(len(ind)))
            vcen_err.append(temp_vcen / np.sqrt(len(ind)))
            rcen_err.append(temp_rcen / np.sqrt(len(ind)))
            xcen.append(temp_rcen * np.sin(np.deg2rad(thea2 + bi / 2.)))
            ycen.append(temp_rcen * np.cos(np.deg2rad(thea2 + bi / 2.)))
    return (zcen, zrms, vcen, vrms, rcen, rrms, thea, z_err, v_err, r_err, zcen_err, vcen_err, rcen_err, xcen, ycen)
def function(x, free_params):
#	width, PHIkink, Rkink, pitchAngle1, pitchAngle2 = fp2p(free_params)
	PHI = x

	#ln(R/Rkink) = -(beta_rad - beta_kink_rad)*tan(pa_rad)
	if len(free_params)>3:	###use kink
		PHIkink, Rkink, pitchAngle1, pitchAngle2 = free_params
		dPHI = (PHI - PHIkink)/180*np.pi
		lnRtoRkink = np.empty(PHI.shape)

		index = PHI >= PHIkink
		lnRtoRkink[index] = -dPHI[index] * np.tan(pitchAngle1/180*np.pi)

		index = PHI < PHIkink
		lnRtoRkink[index] = -dPHI[index] * np.tan(pitchAngle2/180*np.pi)
	else:
		PHIkink, Rkink, pitchAngle1 = free_params
		dPHI = (PHI - PHIkink)/180*np.pi
		lnRtoRkink = -dPHI * np.tan(pitchAngle1/180*np.pi)

	lnR = lnRtoRkink + np.log(Rkink)

	return np.exp(lnR)
def cal_warpc(rg,thea):
	p1=[9.26,17.4,0.148]
	p2=[7.72,17.5,0.06,1.33]
#	p3=[-2.97397887e-02,1.72161063e-02,-1.66075673e+01,7.93265614e-03,6.05519348e+01]
	p4=[-0.07616907,  0.10783739,  7.79992154,  0.9348509,  -9.02005785]  ##1st notcorrect
	p5=[-0.069110292892,0.10257483877,7.868629069,1.0,-8.5999172236, 0.21143504143391126, 16.42369138978075, 1.0, -98.64681680983567]##a0, a1, Rw1, bw1, PHIw1, a2, Rw2, bw2, PHIw2
	zw1, zw2, zw4,zw5 =[], [], [], []
#	thea = np.linspace(-30,170,200)
	zw1=p1[2]*(rg-p1[0])*np.sin(np.deg2rad(thea-p1[1]))   # chenxiaodian all b=1
	zw2=p2[2]*((rg-p2[0])**p2[3])*np.sin(np.deg2rad(thea-p2[1]))  # chenxiaodian all b=1.33
#	zw2=-p3[0]+(rg - 8)**2*(p3[1]*np.sin(np.deg2rad(thea-p3[2]))+p3[3]*np.sin(2.0*(thea-p3[4])/180.*np.pi))  
	zw4=p4[0]+p4[1] * (rg-p4[2])**p4[3] * np.sin((thea-p4[4])/180*np.pi) #mwisp first component
	zw5=np.zeros(thea.size)                                               ##mwisp two component
	if (rg>p5[6]): 
		zw5 = p5[0]+p5[1] * (rg-p5[2])**p5[3] * np.sin((thea-p5[4])/180*np.pi)+p5[5] * (rg-p5[6])**p5[7] * np.sin(2*(thea-p5[8])/180*np.pi)#
	elif (rg<p5[6]):
		zw5 = p5[0]+p5[1] * (rg-p5[2])**p5[3] * np.sin((thea-p5[4])/180*np.pi)  
	return(zw1*1000.,zw2*1000.,zw4*1000.,zw5*1000.)
def cal_warp(rg,thea):
#	rg=rg+0.5    # considering the different rotaion curve
	w1p=[9,197,10,-3.1] # [k0,k1,rk1,k2]
	w1=w1p[0]+w1p[1]*(rg-w1p[2])+w1p[3]*(rg-w1p[2])**2
	w0p=[-66,150,15,-0.47] # [k0,k1,rk1,k2]
	w0=w0p[0]+w0p[1]*(rg-w0p[2])+w0p[3]*(rg-w0p[2])**2
	w2p=[-70,171,15,-5.3] # [k0,k1,rk1,k2]
	w2=w2p[0]+w2p[1]*(rg-w2p[2])+w2p[3]*(rg-w2p[2])**2
	w = [], []
	if ((w2 >= 150) & (rg > 15)):
		w= w0+w1*np.sin(np.deg2rad(thea))+w2*np.sin(np.deg2rad(2.*thea))
	elif ((w2 < 150) & (rg > 15)):
		w= w0+w1*np.sin(np.deg2rad(thea))
	else:
		w= w1*np.sin(np.deg2rad(thea))
	return(w)
if __name__ == '__main__':
	osc_data = np.loadtxt('osc2_para.txt',comments='#')
	osc_mass = osc_data[:,7]/osc_data[:,8]
	osc_rr = osc_data[:,11]
	osc_zz = osc_data[:,14]*1000.
	osc_az = osc_data[:,19]
	osc_vv = osc_data[:,2]
	osc_bb = osc_data[:,1]
	osc_ll = osc_data[:,0]
	best_osc=[47, 16.190233947534733, 12.330164850703794] ##47, 16.060484876948312, 11.996988806408883] #mass
	PHI_osc = np.linspace(-25, 165, 300)
	R_osc = function(PHI_osc, best_osc)
	osc_zh,osc_z1,osc_z2,osc_z4,osc_z5 = [],[],[],[],[]
	for i in range(len(R_osc)):
		z1,z2,z4,z5 = cal_warpc(R_osc[i],PHI_osc[i])
		zh = cal_warp(R_osc[i],PHI_osc[i])
		osc_zh.append(zh)
		osc_z1.append(z1)
		osc_z2.append(z2)
		osc_z4.append(z4)
		osc_z5.append(z5)
    
	osc_zcen,osc_zrms,osc_vcen,osc_vrms,osc_rcen,osc_rrms,osc_thea,osc_zerr,osc_verr,osc_rerr,osc_zcerr,osc_vcerr,osc_rcerr,osc_xcen,osc_ycen=cal_zcen_zrms(osc_az,osc_zz,osc_vv,osc_rr,osc_mass,4)
	maxval = int(np.max(PHI_osc))
#        maxval = 148.
	oscc_thea = np.linspace(-32,maxval,maxval*25+1,endpoint=True)#np.arange(maxval+1)
	f = sp_interp.interp1d(PHI_osc, osc_z4, fill_value='extrapolate')
	oscc_zcen = f(oscc_thea)
	print(oscc_thea)
	osc_zzz = []
	osc_mm = []
	for i in range(len(osc_zz)):
		ind2 = np.where(np.abs(osc_az[i]-oscc_thea)<=0.026)[0]
		if (len(ind2>=1)):
			osc_zzz.append(osc_zz[i]-oscc_zcen[ind2[0]])
			osc_mm.append(osc_mass[i])
			wei_osc = osc_mm/np.linalg.norm(osc_mm)
	osc_dist,osc_colors = mass_distribute(osc_mass)
	pl.rcParams['xtick.direction'] = 'in'
#	pl.tick_params(width)
	fig1, ax = pl.subplots(nrows=2,ncols=1,sharex=True, figsize=[10,13])
	pl.subplots_adjust(wspace=0.1,hspace=0.1)
	pl.subplot(2,1,1)
	pl.scatter(osc_az,osc_zz,s=osc_dist,c='b',alpha=0.2,label='OSC')
	pl.plot(osc_thea,osc_zcen,'-',c='r',lw=4)
	pl.plot(PHI_osc,osc_z1,'-',c='k',label='Chen b=1')
	pl.plot(PHI_osc,osc_zh,'-',c='cyan',label='HI')
	pl.plot(PHI_osc,osc_z4,'-',c='magenta',label='MWISP 1comp')
	pl.plot(PHI_osc,osc_z5,'-',c='orange',label='MWISP 2comp')
	z1=np.array(osc_zcen)+np.array(osc_zrms)*2.355/2.
	z2=np.array(osc_zcen)-np.array(osc_zrms)*2.355/2.
#	pl.plot(osc_thea[3:],osc_zcen[3:],'-',c='g',lw=4)
	z1=np.array(osc_zcen)+np.array(osc_zrms)*2.355/2.
	z2=np.array(osc_zcen)-np.array(osc_zrms)*2.355/2.
	pl.ylim([-700,1500])
	pl.legend()
	pl.subplot(2,1,2)
	print(len(osc_zzz))
	print(len(osc_az))
	osc_zcen,osc_zrms,osc_vcen,osc_vrms,osc_rcen,osc_rrms,osc_thea,osc_zerr,osc_verr,osc_rerr,osc_zcerr,osc_vcerr,osc_rcerr,osc_xcen,osc_ycen=cal_zcen_zrms(osc_az,np.transpose(osc_zzz),osc_vv,osc_rr,osc_mass,4)
	pl.scatter(osc_az,osc_zzz,s=osc_dist,c='b',alpha=0.2,label='out')
	pl.plot(osc_thea,osc_zcen,'-',c='r',lw=4)
	pl.plot([-30,150],[0,0],'--',color='b')
	pl.xlabel('Galactocentric Azimuth (degree)',fontsize=18,fontweight=1.8)
	pl.ylabel('Z scale (pc)',fontsize=18,fontweight=1.8)
	pl.xlim([-30,160])
	pl.ylim([-600,800])
#	pl.yticks([-500,0,500,1000,1500],['-500','0','500','1000','1500'],fontsize=14,fontweight=1.8)
#	pl.xticks([-25,0,25,50,75,100,125,150],['-25','0','25','50','75','100','125','150'],fontsize=14,fontweight=1.8)
	pl.grid(True,ls='--',alpha=0.4)
	fig1.savefig('osc_warp_corrugation.png',format='png',bbox_inches='tight') 
	pl.show()
'''
