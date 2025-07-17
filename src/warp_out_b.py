### plot GL-Z and GL-(Z-warp model) of MCs in Outer arms

import numpy as np
import pylab as pl
import scipy.interpolate as sp_interp
from scipy.stats import norm
import scipy.optimize as sp_opt
from  Distance import KinUnit, Distance
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from shared import *



def getKD(l, b, v):
	skyc = SkyCoord(l = l, b = b, frame = 'galactic', unit = 'deg')
	d = Distance(skyc, v, np.zeros(l.size)).Kdist()
	return d[1]

d2r = np.pi/180
def g2gc(l, b, d):
	#galactic to galactocentric
	X = d * np.cos(b*d2r) * np.cos(l*d2r) - 8.15
	Y = d * np.cos(b*d2r) * np.sin(l*d2r)
	Z = d * np.sin(b*d2r)

	Rgal = np.sqrt(X**2 + Y**2)
	Az = np.arctan2(Y, -X) / d2r #az=0 toward sun and counterclockwise
	return Rgal, Az, Z

def gc2g(Rgal, Az, Z):
	Xs = Rgal * np.cos(Az*d2r) - 8.15
	Ys = Rgal * np.sin(Az*d2r)

	l = np.arctan2(Ys, -Xs) / d2r
	d = np.sqrt(Xs**2 + Ys**2 + Z**2)
	b = np.arcsin(Z/d) / d2r

	l = l%360
	return l,b,d

### get velocity range for baseline
def _getVRange(l):
	if 0<=l<=6: vrange = [-220,+360]
	elif l<=10: vrange = [-150,+280]
	elif l<=20: vrange = [-150,+250]
	elif l<=30: vrange = [-160,+210]
	elif l<=40: vrange = [-170,+200]
	elif l<=50: vrange = [-180,+170]
	elif l<=60: vrange = [-190,+150]
	elif l<=70: vrange = [-200,+120]
	elif l<=80: vrange = [-200,+120]
	elif l<=90: vrange = [-200,+110]
	elif l<=100: vrange = [-200,+110]
	elif l<=110: vrange = [-200,+110]
	elif l<=120: vrange = [-200,+110]
	elif l<=130: vrange = [-200,+110]
	elif l<=140: vrange = [-190,+100]
	elif l<=150: vrange = [-180,+100]
	elif l<=160: vrange = [-160,+100]
	elif l<=170: vrange = [-140,+100]
	elif l<=180: vrange = [-130,+100]
	elif l<=190: vrange = [-120,+110]
	elif l<=200: vrange = [-110,+120]
	elif l<=210: vrange = [-100,+130]
	elif l<=220: vrange = [-100,+150]
	elif l<=230: vrange = [-90,+160]
	elif l<=240: vrange = [-90,+180]
	elif l<=250: vrange = [-90,+190]
	else: vrange = [-100, 100]
	return vrange
from scipy.optimize import root_scalar
from scipy.optimize import fsolve

def inverseKD(Rgal, Az, Z):
	# get l, b, v, d from Rgal, Az, Z
	l, b, d = gc2g(Rgal, Az, Z)
	skyc = SkyCoord(l = l, b = b, frame = 'galactic', unit = 'deg')


	def diffKD(v, skyc, d):
		kd = Distance(skyc, v, 0).Kdist()
		return kd[1]-d

	v = []
	for aSkyc, aD in zip(skyc, d):
		try:
			sol = root_scalar(diffKD, args=(aSkyc, aD), x0=0, method='brentq', bracket=_getVRange(aSkyc.l.degree))
			v.append(sol.root)
		except:
			v.append(np.nan)
	v = np.array(v)

	return l, b, v, d


if __name__ == '__main__':
	figscale = 0.45
	figwidth = textwidth*figscale
	#col_mc = '#4169e1'#'#2f7fa0'

	data = np.loadtxt('out_para.txt',comments='#')
	ll = data[:,0]
	bb = data[:,1]
	vv = data[:,2]
	rr = data[:,11]
	zz = data[:,14]
	az = data[:,19]
	mass = data[:,7]/data[:,8]
	print('Az range:', az.min(), az.max())


	### calculate arm
	azmin = az.min()-3
	azmax = az.max()+15
	PHI = np.linspace(azmin, azmax, 300)
	R = function_arm(PHI, best_out)

	fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(figwidth, figwidth*0.95))
	plt.subplots_adjust(left=0.13, right=0.98, wspace=0.1,hspace=0.06)
	pl.rcParams['xtick.direction'] = 'in'

	### convert r, theta, z of an arm back to l,b,v

	### plot cloud
	scatter_size, scatter_color = mass_distribute(mass)
	ax[0].scatter(ll, bb, s=scatter_size*2, **arm_kws_mc)#, label='Outer')
	ax[0].text(0.02, 0.95, 'Outer', transform=ax[0].transAxes, **arm_kws_text)

	### plot binned average
	lcen, lrms, bcen, brms = bin_data(ll, bb, mass, grid=np.arange(12.5, 250, 2.5), width=5, method='gauss')
	#lcen, lrms, bcen, brms = cal_zcen_zrms(ll, bb, weights=mass, binsize=-5, bin0=12.5)
	ax[0].plot(lcen, bcen, '-', **arm_kws_bin)

	# connect broken with dashed
	idx = np.isfinite(bcen)
	ax[0].plot(lcen[idx], bcen[idx], '--', **arm_kws_bin)

	### calculate warp
	### HI
	Z_hi = cal_warp(R, PHI)
	ll_hi, bb_hi, dd_hi = gc2g(R, PHI, Z_hi)
	ax[0].plot(ll_hi, bb_hi, **arm_kws_hi)
	### warp models
	z1, z2, z4, z5 = cal_warpc(R, PHI)
	l1, b1, d1 = gc2g(R, PHI, z1)
	ax[0].plot(l1, b1, **arm_kws_ceph) #Chen b=1')
	l4, b4, d4 = gc2g(R, PHI, z4)
	ax[0].plot(l4, b4, **arm_kws_co1)
	l5, b5, d5 = gc2g(R, PHI, z5)
	ax[0].plot(l5, b5, **arm_kws_co2)

	ax[0].set_ylim(-5.25, 5.25)
	ax[0].set_xlim(l4.max(), l4.min())
	ax[0].legend(loc = (0.54, 0.02), frameon=False, borderpad=0.2, labelspacing=0.1)
	ax[0].set_ylabel('  (deg)')
	ax[0].text(-0.125, 0.35, 'b', va='center', ha='center', fontsize=20, fontweight='bold', fontfamily='Georgia', fontstyle='italic', transform=ax[0].transAxes, rotation=90)


	# upper tick
	import scipy.interpolate as sp_interp
	phiAxis = np.linspace(azmin, azmax, 200)
	lenAxis = []
	for phi in phiAxis:
		lenAxis.append(cal_arm_length(azmin, phi, p=best_out))
	### interpolate arm length over phi
	f = sp_interp.interp1d(lenAxis, phiAxis, fill_value='extrapolate')

	upper = ax[0].twiny()
	# minor
	lTicks = np.arange(0, 40, 1)
	phiTicks = f(lTicks)
	R = function_arm(phiTicks, best_out)
	lenTick,_,_ = gc2g(R, phiTicks, 0)
	upper.set_xticks(lenTick, minor=True)
	# major
	lTicks = np.arange(0, 36, 5)
	phiTicks = f(lTicks)
	R = function_arm(phiTicks, best_out)
	lenTick,_,_ = gc2g(R, phiTicks, 0)
	upper.set_xticks(lenTick)
	lTickLabels = lTicks.astype(str)
	lTickLabels[-1] = '  '+lTickLabels[-1]+' kpc'
	upper.set_xticklabels(lTickLabels) # add kpc at the end
	upper.set_xlim(l4.max(), l4.min())


	### plot residual
	#zz_res = zz - function_warp((rr, az))
	zz_res = zz - function_warp( (function_arm(az, best_out), az) )
	ll_res, bb_res, dd_res = gc2g(rr, az, zz_res)
	ax[1].scatter(ll_res, bb_res, s=scatter_size*2, **arm_kws_mc)

	### plot binned average
	lcen, lrms, bcen_res, brms_res = bin_data(ll, bb_res, mass, grid=np.arange(12.5, 250, 2.5), width=5, method='gauss')
	#lcen, lrms, bcen_res, brms_res = cal_zcen_zrms(ll_res, bb_res, weights=mass, binsize=-5, bin0=12.5)
	ax[1].plot(lcen, bcen_res, '-', **arm_kws_bin)

	# connect broken with dashed
	idx = np.isfinite(bcen)
	ax[1].plot(lcen[idx], bcen_res[idx], '--', **arm_kws_bin)

	### plot H line
	arm_kws_co1['zorder']=0
	ax[1].plot([l4[0], l4[-1]], [0, 0], **arm_kws_co1)

	#ax[1].set_xlabel('Galactic Longitude (deg)', fontsize=15, fontweight='bold')#, fontfamily='Georgia', fontstyle='italic')
	#ax[1].set_ylabel('Galactic Latitude (deg)     ', fontsize=15, fontweight='bold')
	ax[1].set_xlabel('  (deg)')
	ax[1].text(0.43, -0.17, 'l', va='center', ha='center', fontsize=20, fontweight='bold', fontfamily='Georgia', fontstyle='italic', transform=ax[1].transAxes)

	ax[1].set_ylabel('$\mathbf{\Delta}$   (deg)')
	ax[1].text(-0.125, 0.40, 'b', va='center', ha='center', fontsize=20, fontweight='bold', fontfamily='Georgia', fontstyle='italic', transform=ax[1].transAxes, rotation=90)

	ax[1].set_ylim([-5.25, 5.25])
	ax[1].grid(True, ls='--', alpha=0.4)


	### insert a arm plot
	insert_arm_plot(ax[0], [0.26, 0.02, 0.36, 0.36], boldArm='out', boldRange=[azmin,azmax])


	if 0:
		ax[1].set_ylim([-4.75, 6.75])
		### insert a power spectrum
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
		powerAxes.plot([bestFreq, bestFreq, bestFreq+0.1], [peak+0.01, peak+0.03, peak+0.03], 'k')
		powerAxes.text(bestFreq+0.11, peak+0.03, '%.1f kpc' % (1/bestFreq), ha='left', va='center')
		# find the second period
		peak = np.max(power[20:])
		bestFreq = frequency[np.argmax(power[20:])+20]
		powerAxes.plot([bestFreq, bestFreq, bestFreq+0.1], [peak+0.01, peak+0.03, peak+0.03], 'k')
		powerAxes.text(bestFreq+0.11, peak+0.03, '%.1f kpc' % (1/bestFreq), ha='left', va='center')

		#powerAxes.set_xscale('log')
		#powerAxes.set_xlim(2, 50)
		powerAxes.set_ylim(0, np.max(power)*1.5)
		powerAxes.tick_params(top=False, left=True, direction='in', pad=2)
		powerAxes.minorticks_on()
		powerAxes.set_xlabel('frequency (kpc$^{-1}$)', labelpad=-2)
		powerAxes.set_ylabel('Norm\nPower', labelpad=0)
		powerAxes.spines['top'].set_visible(False)
		#powerAxes.spines['left'].set_visible(False)
		powerAxes.spines['right'].set_visible(False)
		#powerAxes.yaxis.set_visible(False)
		powerAxes.patch.set_alpha(0.0)


	ax[0].text(-0.15, 1.0, 'a', ha='left', va='top', color='black', font=subfigureIndexFont, transform=ax[0].transAxes)
	ax[1].text(-0.15, 1.0, 'b', ha='left', va='top', color='black', font=subfigureIndexFont, transform=ax[1].transAxes)

	plt.savefig('fig/out_warp_corrugation_b.%s' % (mpl.rcParams['savefig.format']), bbox_inches='tight')
	plt.savefig('fig/out_warp_corrugation_b.png', bbox_inches='tight')
	plt.show()

'''
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

def cal_zcen_zrms(az, zz, bb, vv, rr, mass, bi):
    zcen, zrms, bcen,brms,vcen, vrms, rcen, rrms, thea, z_err, b_err,v_err, r_err, zcen_err, vcen_err, rcen_err, xcen, ycen = [],[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []
    for i in range(60):
        thea1 = 160 - i * bi
        thea2 = 160 - (i + 1) * bi
        ind = np.where((az > thea2) & (az < thea1))[0]  # Extract the array of indices
        if len(ind) > 5:
            temp_zcen, temp_zrms = weighted_avg_and_std(zz[ind], mass[ind])
            temp_bcen, temp_brms = weighted_avg_and_std(bb[ind], mass[ind])
            temp_vcen, temp_vrms = weighted_avg_and_std(vv[ind], mass[ind])
            temp_rcen, temp_rrms = weighted_avg_and_std(rr[ind], mass[ind])
            zcen.append(temp_zcen)
            zrms.append(temp_zrms)
            bcen.append(temp_bcen)
            brms.append(temp_brms)            
            vcen.append(temp_vcen)
            vrms.append(temp_vrms)
            rcen.append(temp_rcen)
            rrms.append(temp_rrms)
            thea.append((thea2 + bi / 2.))
            z_err.append(temp_zrms / np.sqrt(len(ind)))
            b_err.append(temp_brms / np.sqrt(len(ind)))            
            v_err.append(temp_vrms / np.sqrt(len(ind)))
            r_err.append(temp_rrms / np.sqrt(len(ind)))
            zcen_err.append(np.abs(temp_zcen) / np.sqrt(len(ind)))
            vcen_err.append(temp_vcen / np.sqrt(len(ind)))
            rcen_err.append(temp_rcen / np.sqrt(len(ind)))
            xcen.append(temp_rcen * np.sin(np.deg2rad(thea2 + bi / 2.)))
            ycen.append(temp_rcen * np.cos(np.deg2rad(thea2 + bi / 2.)))
    return (zcen, zrms, bcen,brms,vcen, vrms, rcen, rrms, thea, z_err, b_err,v_err, r_err, zcen_err, vcen_err, rcen_err, xcen, ycen)
def cal_vcen_vrms(ll,bb,vv,zz,mass,bi):
    vcen,vrms,bcen,brms,zcen,zrms,lon,berr,verr,zerr = [],[],[],[],[],[],[],[],[],[]
    for i in range(44):
        l1 = 12.5+i*bi
        l2 = 12.5+(i+1)*bi
        ind = np.where((ll>l1) & (ll<l2))
        if (len(np.transpose(ind))>=1):   ##>= 3 in each bin
            temp_bcen,temp_brms = weighted_avg_and_std(bb[ind],mass[ind])
            bcen.append(temp_bcen)
            brms.append(temp_brms)
            berr.append(temp_brms/np.sqrt(len(np.transpose(ind))))
            temp_vcen,temp_vrms = weighted_avg_and_std(vv[ind],mass[ind])
            vcen.append(temp_vcen)
            vrms.append(temp_vrms)
            lon.append((l1+bi/2.))
            verr.append(temp_vrms/np.sqrt(len(np.transpose(ind))))
            temp_zcen,temp_zrms = weighted_avg_and_std(zz[ind],mass[ind])
            zcen.append(temp_zcen)
            zrms.append(temp_zrms)
            zerr.append(temp_zrms/np.sqrt(len(np.transpose(ind))))            
    return(lon,bcen,brms,berr,vcen,vrms,verr,zcen,zrms,zerr)

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
	return(zw1,zw2,zw4,zw5)
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
	return(w/1000.)
'''
'''
	### plot 
	plt.plot(l,b,'--',c=col_co2,lw=1.5)
	plt.plot(out_lon,out_bcen,'-',c='r',lw=4)
#	plt.plot(PHI_out,out_z1,'-',c='k',label='Chen b=1')
#	plt.plot(PHI_out,out_zh,'-',c='cyan',label='HI')
#	plt.plot(PHI_out,out_z4,'-',c='magenta',label='MWISP 1comp')
#	plt.plot(PHI_out,out_z5,'-',c='orange',label='MWISP 2comp')
	b1=np.array(out_bcen)+np.array(out_brms)*2.355/2.
	b2=np.array(out_bcen)-np.array(out_brms)*2.355/2.

#	pl.plot(out_thea[3:],out_zcen[3:],'-',c='g',lw=4)
	b1=np.array(out_bcen)+np.array(out_brms)*2.355/2.
	b2=np.array(out_bcen)-np.array(out_brms)*2.355/2.

	pl.ylim([-5.25,5.25])
	pl.xlim([230,15])
	pl.legend()


	
	best_per=[30.03, 10.062122711773025, 9.839300621559302]  #mass	
	best_out=[20.25195158643863, 13.259894850233614, 11.05257961109539, 3.500917359458265]  ##mass
	best_osc=[47, 16.190233947534733, 12.330164850703794] ##47, 16.060484876948312, 11.996988806408883] #mass
	PHI_per = np.linspace(-22, 80, 200)
	PHI_out = np.linspace(-25, 155, 300)
	PHI_osc = np.linspace(-25, 155, 300)
	R_per = function(PHI_per, best_per)
	R_out = function(PHI_out, best_out)
	R_osc = function(PHI_osc, best_osc)
	out_zh,out_z1,out_z2,out_z4,out_z5 = [],[],[],[],[]
	for i in range(len(R_out)):
		z1,z2,z4,z5 = cal_warpc(R_out[i],PHI_out[i])
		zh = cal_warp(R_out[i],PHI_out[i])
		out_zh.append(zh)
		out_z1.append(z1)
		out_z2.append(z2)
		out_z4.append(z4)
		out_z5.append(z5)
    


	out_zcen,out_zrms,out_bcen,out_brms,out_vcen,out_vrms,out_rcen,out_rrms,out_thea,out_zerr,out_berr,out_verr,out_rerr,out_zcerr,out_vcerr,out_rcerr,out_xcen,out_ycen=cal_zcen_zrms(out_az,out_zz,out_bb,out_vv,out_rr,out_mass,4)
	out_lon,out_bcen,out_brms,out_berr,out_vcen,out_vrms,out_verr,out_zcen,out_zrms,out_zerr=cal_vcen_vrms(out_ll,out_bb,out_vv,out_zz,out_mass,5)
	maxval = int(np.max(PHI_out))
	outt_thea = np.linspace(-30,maxval,maxval*25+1,endpoint=True)#np.arange(maxval+1)
	outt_ll = np.linspace(15,230,230*25+1,endpoint=True)#np.arange(maxval+1)
	l, b, v, d = inverseKD(np.array(R_out), np.array(PHI_out), np.array(out_z5))
	f = sp_interp.interp1d(l, b, fill_value='extrapolate')
	outt_bb = f(outt_ll)
	out_bbb = []
	out_mm = []
	for i in range(len(out_bb)):
		ind2 = np.where(np.abs(out_ll[i]-outt_ll)<=0.0243)[0]
		if (len(ind2>=1)):
			out_bbb.append(out_bb[i]-outt_bb[ind2[0]])
			out_mm.append(out_mass[i])
			wei_out = out_mm/np.linalg.norm(out_mm)
	out_dist,out_colors = mass_distribute(out_mass)

	pl.rcParams['xtick.direction'] = 'in'
#	pl.tick_params(width)
	fig1, ax = pl.subplots(nrows=2,ncols=1,sharex=True, figsize=[10,13])
	pl.subplots_adjust(wspace=0.1,hspace=0.1)
	pl.subplot(2,1,1)
	pl.scatter(out_ll,out_bb,s=out_dist,c='g',alpha=0.2,label='Out')
	
	pl.plot(l,b,'--',c=col_co2,lw=1.5)
	pl.plot(out_lon,out_bcen,'-',c='r',lw=4)
#	pl.plot(PHI_out,out_z1,'-',c='k',label='Chen b=1')
#	pl.plot(PHI_out,out_zh,'-',c='cyan',label='HI')
#	pl.plot(PHI_out,out_z4,'-',c='magenta',label='MWISP 1comp')
#	pl.plot(PHI_out,out_z5,'-',c='orange',label='MWISP 2comp')
	b1=np.array(out_bcen)+np.array(out_brms)*2.355/2.
	b2=np.array(out_bcen)-np.array(out_brms)*2.355/2.

#	pl.plot(out_thea[3:],out_zcen[3:],'-',c='g',lw=4)
	b1=np.array(out_bcen)+np.array(out_brms)*2.355/2.
	b2=np.array(out_bcen)-np.array(out_brms)*2.355/2.

	pl.ylim([-5.25,5.25])
	pl.xlim([230,15])
	pl.legend()
	


	pl.subplot(2,1,2)
	#out_zcen,out_zrms,out_vcen,out_vrms,out_rcen,out_rrms,out_thea,out_zerr,out_verr,out_rerr,out_zcerr,out_vcerr,out_rcerr,out_xcen,out_ycen=cal_zcen_zrms(out_az,np.transpose(out_zzz),out_vv,out_rr,out_mass,4)
	out_lon,out_bcen,out_brms,out_berr,out_vcen,out_vrms,out_verr,out_zcen,out_zrms,out_zerr=cal_vcen_vrms(out_ll,np.transpose(out_bbb),out_vv,out_zz,out_mass,5)
	pl.scatter(out_ll,out_bbb,s=out_dist,c='g',alpha=0.2,label='out')
	pl.plot(out_lon,out_bcen,'-',c='r',lw=4)
	pl.plot([230,15],[0,0],'--',color='b')

	pl.xlabel('Galactic Longitude (degree)',fontsize=18,fontweight=1.8)
	pl.ylabel('Galactic Latitude (degree)',fontsize=18,fontweight=1.8)
	pl.xlim([230,15])
	pl.ylim([-5.25,5.25])
#	pl.yticks([-500,0,500,1000,1500],['-500','0','500','1000','1500'],fontsize=14,fontweight=1.8)
#	pl.xticks([-25,0,25,50,75,100,125,150],['-25','0','25','50','75','100','125','150'],fontsize=14,fontweight=1.8)
	pl.grid(True,ls='--',alpha=0.4)
	#fig1.savefig('out_warp_corrugation_b.png',format='png',bbox_inches='tight') 
	
	pl.show()
'''