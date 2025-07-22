### params shared by other scripts

import numpy as np
import matplotlib as mpl
### Warp Params
p_1comp = [0,  0.09362826, 8.56849936, 1.05015227, -0.7049789 ]##1st exclude 160-200
#[-0.0717748888594908, 0.1074529579441738, 7.8409214840721795, 0.9370789056657302, -8.059413168263955]   
p_2comp = [0, 1.19228266e-01, 9.01450681e+00, -3.53824563e+00, -1.24059941e-02, 1.27147257e+01, 2.03108539e+00, -2.04459975e+01]##2nd exclude 160-200
#[0, 1.19210237e-01, 9.01434283e+00, -3.53962133e+00, -1.24379884e-02, 1.27166892e+01, 2.03007431e+00, -2.04479290e+01]## not properly burn-in


p_sin1comp = [0.20052978, 7.83808911, 2.76825318, 0.17047776, 0.01625255]	# 1comp+sin
#[0.2012454, 7.8246938, 2.78618375, 0.16944783, 0.01766535] # not properly burn-in
#p_sin1comp = [0.19725642, 7.49363834, 3.69774322, 0.01476242] #1comp+sin(fix P1=0)
p_sin2comp = [0.19866392, 7.95872642, 2.61382166, 0.20077682, 0.02240027] # 2comp+sin
#[0.19869108, 7.958983, 2.61323087, 0.20085326, 0.02240435] # use not properly burn-in 2comp result

### Arm Params
best_per=[30.03, 10.062122711773025, 9.839300621559302]  #mass  
best_out=[20.25195158643863, 13.259894850233614, 11.05257961109539, 3.500917359458265]  ##mass
best_osc=[47, 16.190233947534733, 12.330164850703794] ##47, 16.060484876948312, 11.996988806408883] #mass


colorText = '#1e3f66'
#HI
colorHI = 'gray'
styleHI = ':'
labelHI = 'HI'
#Cepheids
colorCeph = 'gray'
styleCeph = (0,(6,3))#'--'
labelCeph = 'Cepheids'
#CO 1comp
colorCO1 = 'coral'
styleCO1 = ':'
labelCO1 = 'CO $m=1$'
#CO 2comp
colorCO2 = 'coral'
styleCO2 = (0,(6,3))#'--'
labelCO2 = 'CO $m=1,2$'


def darkerColor(color_name, factor=0.5):
	### 将一个颜色变暗
	import matplotlib.colors as mcolors
	r, g, b = mcolors.to_rgb(color_name)
	r = max(0, min(1, r*factor))
	g = max(0, min(1, g*factor))
	b = max(0, min(1, b*factor))
	return (r, g, b)

### kws in ring plots: warp1.py / warp1_xf.py / warp_resi.py
ring_kws_mc   = dict(color='#2074b0', edgecolors='none', alpha=0.2, zorder=5)
ring_kws_bin  = dict(color=darkerColor('#2074b0'),  lw=2, alpha=0.8, zorder=20)#5ABEF0
ring_kws_hi   = dict(color=colorHI,   ls=styleHI,   lw=1,   label=labelHI,   zorder=11)
ring_kws_ceph = dict(color=colorCeph, ls=styleCeph, lw=1,   label=labelCeph, zorder=12)
ring_kws_co1  = dict(color=colorCO1,  ls=styleCO1,  lw=0.8, label=labelCO1,  zorder=13)
ring_kws_co2  = dict(color=colorCO2,  ls=styleCO2,  lw=0.8, label=labelCO2,  zorder=14)
ring_kws_text = dict(color=colorText, ha='left', va='top', fontsize=20, fontweight='bold')

### kws in radial plots: corrugation1b.py, corrugation_resi.py, corrugation_resi_xf.py
rad_kws_mc   = dict(color='#008080', edgecolors='none', alpha=0.3, zorder=5)
rad_kws_err  = dict(color=darkerColor('#008080'), fmt='.', markersize=5., elinewidth = 1, capsize=2, zorder=10)
rad_kws_hi   = dict(color=colorHI,   ls=styleHI,   lw=1.5, label=labelHI,   zorder=11)
rad_kws_ceph = dict(color=colorCeph, ls=styleCeph, lw=1.5, label=labelCeph, zorder=12)
rad_kws_co1  = dict(color=colorCO1,  ls=styleCO1,  lw=1,   label=labelCO1,  zorder=13)
rad_kws_co2  = dict(color=colorCO2,  ls=styleCO2,  lw=1,   label=labelCO2,  zorder=14)
rad_kws_sin  = dict(color=colorCO1,  ls=':',       lw=3,   label='SIN comp',zorder=11)
rad_kws_text = dict(color=colorText, ha='left', va='top', fontsize=20, fontweight='bold')

### kws in arm plots: warp_out.py / warp_per.py / warp_osc.py / warp_osc_b.py
arm_kws_mc   = dict(color='#4169e1', edgecolors='none', alpha=0.2, zorder=5)
arm_kws_bin  = dict(color=darkerColor('#4169e1'),  lw=2, alpha=0.8, zorder=20)
arm_kws_hi   = dict(color=colorHI,   ls=styleHI,   lw=1.5, label=labelHI,   zorder=11)
arm_kws_ceph = dict(color=colorCeph, ls=styleCeph, lw=1.5, label=labelCeph, zorder=12)
arm_kws_co1  = dict(color=colorCO1,  ls=styleCO1,  lw=1,   label=labelCO1,  zorder=13)
arm_kws_co2  = dict(color=colorCO2,  ls=styleCO2,  lw=1,   label=labelCO2,  zorder=14)
arm_kws_text = dict(color=colorText, ha='left', va='top', fontsize=20, fontweight='bold')

textwidth=18 #inches, width of text in latex
subfigureIndexFont = dict(family="Arial Black", size=24)
mpl.rcParams['savefig.format'] = 'pdf'
mpl.rcParams['savefig.dpi'] = 400
mpl.rcParams['axes.linewidth'] = 0.8
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.labelweight'] = 'bold'
mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True
mpl.rcParams['legend.fontsize'] = 15

def funciton_loglinear(x,a,b):
	y =  np.log(a)-x/b
	return(y)


def function_arm(PHI, p):
	### arm function, PHI in degree, return R in kpc
	PHI = np.array(PHI)
	if len(p)>3:  ###use kink
		PHIkink, Rkink, pitchAngle1, pitchAngle2 = p
		dPHI = (PHI - PHIkink)/180*np.pi
		lnRtoRkink = np.empty(PHI.shape)

		index = PHI >= PHIkink
		lnRtoRkink[index] = -dPHI[index] * np.tan(pitchAngle1/180*np.pi)

		index = PHI < PHIkink
		lnRtoRkink[index] = -dPHI[index] * np.tan(pitchAngle2/180*np.pi)
	else:
		PHIkink, Rkink, pitchAngle1 = p
		dPHI = (PHI - PHIkink)/180*np.pi
		lnRtoRkink = -dPHI * np.tan(pitchAngle1/180*np.pi)

	lnR = lnRtoRkink + np.log(Rkink)

	return np.exp(lnR)


def function_warp(x, p=p_1comp, sin=None):
	### warp function, return Z
	R, PHI = x
	a0 = p[0]

	if isinstance(R, (int, float)):
		R = np.ones(PHI.shape) * R
	elif isinstance(PHI, (int, float)):
		PHI = np.ones(R.shape) * PHI
	z = np.ones(R.shape) * a0

	if len(p)<=5: #1 comp
		a1, Rw1, bw1, PHIw1 = p[1:]
		idx = R>Rw1
		z[idx] += a1 * (R[idx]-Rw1)**bw1 * np.sin(np.deg2rad(PHI[idx]-PHIw1))
	else: #2comp
		a1, Rw1, PHIw1, a2, Rw2, bw2, PHIw2 = p[1:]
		idx = R>Rw1
		z[idx] += a1 * (R[idx]-Rw1) * np.sin(np.deg2rad(PHI[idx]-PHIw1))
		idx = R>Rw2
		z[idx] += a2 * (R[idx]-Rw2)**bw2 * np.sin(2*np.deg2rad(PHI[idx]-PHIw2))

	if sin is not None:
		PHIsin, PHIsinwid = 33.5, 4
		asin, Rsin, Period0, Period1, a4 = sin
		z[idx] += (a4*(R[idx] - Rsin) + asin*np.sin( (R[idx] - Rsin) / (Period0 + Period1*(R[idx]-Rsin)) * 2*np.pi )) * np.exp(-(PHI[idx]-PHIsin)**2/ 2 / PHIsinwid**2)
	return z


def mass_distribute(mass):
	### regulate mass as size and color of scatters
	normSize = np.zeros(len(mass))
	normColor = np.empty(len(mass), dtype=str)

	idx = mass<100
	normSize[idx] = 0.5
	normColor[idx] = 'purple'

	idx = (mass >= 100) & (mass < 1000)
	normSize[idx] = 1
	normColor[idx] = 'blue'

	idx = (mass >= 1000) & (mass < 10000)
	normSize[idx] = 10
	normColor[idx] = 'green'

	idx = (mass >= 10000) & (mass < 100000)
	normSize[idx] = 100
	normColor[idx] = 'red'

	idx = mass >= 100000
	normSize[idx] = 200
	normColor[idx] = 'cyan'

	return normSize, normColor


def cal_warp(R = np.linspace(7.5,22,200), PHI = np.linspace(-30,170,200), scale=0.95):
	# regulate R, PHI *95%
	if isinstance(R, (int, float)):
		R = np.ones(PHI.shape) * R
	else:
		PHI = np.ones(R.shape) * PHI

	Rs=R/scale

	#   R=R+0.5    # considering the different rotaion curve
	w1p = [9, 197, 10, -3.1] # [k0,k1,rk1,k2]
	w1 = w1p[0] + w1p[1]*(Rs-w1p[2]) + w1p[3]*(Rs-w1p[2])**2
	w0p = [-66, 150, 15, -0.47] # [k0,k1,rk1,k2]
	w0 = w0p[0] + w0p[1]*(Rs-w0p[2]) + w0p[3]*(Rs-w0p[2])**2
	w2p = [-70,171,15,-5.3] # [k0,k1,rk1,k2]
	w2 = w2p[0] + w2p[1]*(Rs-w2p[2]) + w2p[3]*(Rs-w2p[2])**2

	w = w1*np.sin(np.deg2rad(PHI))

	idx = (w2 >= 150) & (Rs > 15)
	w[idx] = (w0 + w1*np.sin(np.deg2rad(PHI)) + w2*np.sin(np.deg2rad(2.*PHI)))[idx]

	idx = (w2 < 150) & (Rs > 15)
	w[idx] = (w0+w1*np.sin(np.deg2rad(PHI)))[idx]
	'''
	if (w2 >= 150) & (R > 15):
		w = w0 + w1*np.sin(np.deg2rad(PHI)) + w2*np.sin(np.deg2rad(2.*PHI))
	elif (w2 < 150) & (R > 15):
		w = w0+w1*np.sin(np.deg2rad(PHI))
	else:
		w = w1*np.sin(np.deg2rad(PHI))
	'''
	return w/1000. *scale


def cal_warpc(R = np.linspace(7.5,22,200), PHI = np.linspace(-30,170,200)):
	### calculate z at certain R and PHI in kpc
	p1=[9.26, 17.4, 0.148]
	p2=[7.72, 17.5, 0.06, 1.33]

	zw1, zw2, zw4,zw5 =[], [], [], []
	zw1=p1[2] * (R-p1[0]) * np.sin(np.deg2rad(PHI-p1[1]))   # chenxiaodian all b=1
	zw2=p2[2] * (R-p2[0])**p2[3] * np.sin(np.deg2rad(PHI-p2[1]))  # chenxiaodian all b=1.33
	#zw3=-p3[0]+(R - 8)**2*(p3[1]*np.sin(np.deg2rad(PHI-p3[2]))+p3[3]*np.sin(2.0*(PHI-p3[4])/180.*np.pi))  
	zw4=function_warp((R, PHI), p_1comp) ##mwisp first component
	zw5=function_warp((R, PHI), p_2comp) #mwisp two component
	#if (R>p5[6]):
	#    zw5 = p5[0]+p5[1] * (R-p5[2])**p5[3] * np.sin((PHI-p5[4])/180*np.pi)+p5[5] * (R-p5[6])**p5[7] * np.sin(2*(PHI-p5[8])/180*np.pi)#
	#elif (R<p5[6]):
	#    zw5 = p5[0]+p5[1] * (R-p5[2])**p5[3] * np.sin((PHI-p5[4])/180*np.pi)  
	return zw1, zw2, zw4, zw5


def cal_warpmwisp(R):
	a0, a, Rw, PHIw, b=[1.01,2.62,8.8,120.38,1.34] 
	zw=a0 + a * (R - Rw)**b * np.sin((PHI-PHIw) / 180 * np.pi)
	return(PHI,zw)


def weighted_avg_and_std(values, weights):
	average = np.average(values, weights=weights)
	# Fast and numerically precise:
	variance = np.average((values-average)**2, weights=weights)
	return (average, np.sqrt(variance))


def bin_data(x, y, w, grid=None, width=1, method='median'):
	xmean = []
	xstd = []
	ymean = []
	ystd = []
	for xc in grid:
		idx = (x>=xc-width/2) & (x<xc+width/2)
		if idx.sum() >= 3:
			if method == 'median':
				wc = w * idx.astype(int)
			elif method == 'gauss':
				wc = w * np.exp(-(x - xc)**2 / 2 / (width/2.355)**2)
			xm = xc#np.average(x, weights=wc)
			xs = 0#np.sqrt( np.average((x-xm)**2, weights=wc) )
			ym = np.average(y, weights=wc)
			ys = np.sqrt( np.average((y-ym)**2, weights=wc) )
		else:
			xm = xc
			xs = 0
			ym = np.nan
			ys = np.nan
		xmean.append( xm )
		xstd.append( xs )
		ymean.append( ym )
		ystd.append( ys )
	return np.array(xmean), np.array(xstd), np.array(ymean), np.array(ystd)


def cal_zcen_zrms(az, *zz, weights=1, binsize = 6, binstep = None, bin0=160):
	### bin zz, positive binsize means counting downward
	bin_az, bin_zz= [], []

	if binstep is None: binstep = binsize
	if binstep>0: nbin = (bin0 - np.nanmin(az)) // binstep +1
	else: nbin = (bin0 - np.nanmax(az)) // binstep +1
	bincen = bin0 - binstep*np.arange(nbin)

	zz_mean = []
	zz_std = []

	for b in bincen:
		PHI1 = b - abs(binsize)/2
		PHI2 = b + abs(binsize)/2
		idx = (az >= PHI1) & (az < PHI2)
		#idx = (az >=min(PHI1, PHI2)) & (az < max(PHI1, PHI2))
		#idx = (az >= PHI2) & (az < PHI1)

		if idx.sum() > 2:
			bin_az.append(((PHI1+PHI2)/2, 0))
			#bin_az.append(weighted_avg_and_std(az[idx], mass[idx]))
			bin_zz.append([])
			for z in zz:
				avg, std = weighted_avg_and_std(z[idx], weights[idx])
				bin_zz[-1].append(avg)
				bin_zz[-1].append(std)
		else:
			bin_az.append(((PHI1+PHI2)/2, 0))
			bin_zz.append([np.nan]*2*len(zz))

	bin_az = np.array(bin_az).T
	bin_zz = np.array(bin_zz).T
	return *bin_az, *bin_zz#azcen, azrms, zcen, zrms, vcen, vrms, rcen, rrms
	'''
	### bi is the binsize
	zcen,zrms,vcen,vrms,rcen,rrms,azcen,azrms,z_err,v_err,r_err,az_err,zcen_err,vcen_err,rcen_err,azcen_err = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
	for i in range(42):
		PHI1 = 160-i*bi
		PHI2 = 160-(i+1)*bi


		ind = np.where((az>PHI2) & (az<PHI1))
		if (len(np.transpose(ind))>3):
			temp_zcen,temp_zrms = weighted_avg_and_std(zz[ind],mass[ind])
			temp_vcen,temp_vrms = weighted_avg_and_std(vv[ind],mass[ind])
			temp_rcen,temp_rrms = weighted_avg_and_std(rr[ind],mass[ind])
			temp_azcen,temp_azrms = weighted_avg_and_std(az[ind],mass[ind])
#            temp_zcen = np.average(zz[ind],weights=mass[ind])
#            temp_zrms = np.sqrt(np.average((zz[ind]-temp_zcen)**2, weights=mass[ind]))
			zcen.append(temp_zcen)
			zrms.append(temp_zrms)
			vcen.append(temp_vcen)
			vrms.append(temp_vrms)
			rcen.append(temp_rcen)
			rrms.append(temp_rrms)
			azcen.append(temp_azcen)
			azrms.append(temp_azrms)
#            PHI.append((PHI2+bi/2.))
			z_err.append(temp_zrms/np.sqrt(len(np.transpose(ind))))
			v_err.append(temp_vrms/np.sqrt(len(np.transpose(ind))))
			r_err.append(temp_rrms/np.sqrt(len(np.transpose(ind))))
			az_err.append(temp_azrms/np.sqrt(len(np.transpose(ind))))
			zcen_err.append(np.abs(temp_zcen)/np.sqrt(len(np.transpose(ind))))
			vcen_err.append(temp_vcen/np.sqrt(len(np.transpose(ind))))
			rcen_err.append(temp_rcen/np.sqrt(len(np.transpose(ind))))
			azcen_err.append(temp_azcen/np.sqrt(len(np.transpose(ind))))
#            xcen.append(temp_rcen*np.sin(np.deg2rad(PHI2+bi/2.)))
#            ycen.append(temp_rcen*np.cos(np.deg2rad(PHI2+bi/2.)))
	return(zcen,zrms,vcen,vrms,rcen,rrms,azcen,azrms,z_err,v_err,r_err,az_err,zcen_err,vcen_err,rcen_err,azcen_err)
	'''

def cal_arm_length(phi0, phi1, p=best_out):
	from scipy.integrate import quad
	def dR_dPHI(PHI, h=1e-5):
		return (function_arm(PHI + h, p) - function_arm(PHI - h, p)) / (2 * h)
	
	def dL(PHI):
		return np.sqrt(dR_dPHI(PHI)**2 + function_arm(PHI, p)**2) / 180*np.pi

	return quad(dL, phi0, phi1)[0]


def insert_arm_plot(ax, bbox, boldArm='per', boldRange=[-50, 100]):
	### a sketch of arm in ax[0]
	from matplotlib.patches import Ellipse
	sk = ax.inset_axes(bbox)

	def plotArm(ax, phiRange, p, **kws):
		phi = np.linspace(*phiRange, 1400)
		r = function_arm(phi, p)
		x = r*np.sin(phi/180*np.pi)
		y = r*np.cos(phi/180*np.pi)
		ax.plot(x, y, **kws)

	arm_kws = dict(lw=0.5)
	plotArm(sk, (-170, 220), [0, 6.5, 12], color='k', **arm_kws) #sag
	plotArm(sk, ( -50, 45), [0, 8.15, 12], color='pink', **arm_kws) #loc
	plotArm(sk, ( -50, 270), [0, 10, 12], color=(0,0.2,0), **arm_kws) #per
	plotArm(sk, ( -48, 420), [0, 13, 12], color=(0.4,0,0), **arm_kws) #out
	plotArm(sk, ( -45, 450), [0, 19, 12], color=(0,0,0.4), **arm_kws) #osc

	if boldArm == 'per': p = [0, 10, 12]
	elif boldArm == 'out': p = [0, 13, 12]
	elif boldArm == 'osc': p = [0, 19, 12]
	plotArm(sk, boldRange, p, color=arm_kws_mc['color'], lw=3, alpha=0.7)

	ell = Ellipse(xy = (0,0), width=3, height=7, angle=-30.5, facecolor='gray', edgecolor='none')
	sk.add_patch(ell)
	sk.plot(0, 8.15, 'r.', ms=4)
	sk.plot(0, 0, 'k.', ms=4)

	sk.set_aspect('equal')
	sk.axis('off')
	sk.set_xlim(-15, 15)
	sk.set_ylim(-12, 20)


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


### import OB star catalog
from astropy.io import fits
hdu = fits.open('OBstars.fits')
data = hdu[1].data
starL = data['GLON']
starB = data['GLAT']
starD = data['DistBJ']/1e3
idx = np.abs(starB)<5.25
starL = starL[idx]
starB = starB[idx]
starD = starD[idx]
starR, starAz, starZ = g2gc(starL, starB, starD)
starX = starD * np.cos(starB*d2r) * np.cos(starL*d2r) - 8.15
starY = starD * np.cos(starB*d2r) * np.sin(starL*d2r)


### digit to unicode of subscript
def to_subscript(n):
    sub_map = str.maketrans('0123456789', '₀₁₂₃₄₅₆₇₈₉')
    return str(n).translate(sub_map)

'''
def readCepheidCat(cat='allGalCep.listID'):
	import pandas as pd
	df = pd.read_csv(cat, comment='#', sep='\s+')
	l = df[]
	return df


a = readCepheidCat()
print(a.shape)
'''