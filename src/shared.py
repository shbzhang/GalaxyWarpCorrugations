### params shared by other scripts

import numpy as np
import matplotlib as mpl
import pandas as pd
### Warp Params
#p_1comp = [0,  0.09362826, 8.56849936, 1.05015227, -0.7049789 ]##1st exclude 160-200, #1 revision
p_1comp = [0,  0.11025648, 8.6446515, 0.97143317, -1.39623657] # 2 revision


#p_2comp = [0, 1.19228266e-01, 9.01450681e+00, -3.53824563e+00, -1.24059941e-02, 1.27147257e+01, 2.03108539e+00, -2.04459975e+01] # 2nd exclude 160-200, #1 revision
#[0, 1.19210237e-01, 9.01434283e+00, -3.53962133e+00, -1.24379884e-02, 1.27166892e+01, 2.03007431e+00, -2.04479290e+01]## not properly burn-in
p_2comp = [0,  0.13201435, 9.13772156, -5.21004583, -0.01892544, 12.20359468, 1.83408649, -18.6826663] # 2 revision

# sin with damp period
p_sin1comp = [0.20052978, 7.83808911, 2.76825318, 0.17047776, 0.01625255]	# 1comp+sin
#[0.2012454, 7.8246938, 2.78618375, 0.16944783, 0.01766535] # not properly burn-in
#p_sin1comp = [0.19725642, 7.49363834, 3.69774322, 0.01476242] #1comp+sin(fix P1=0)
p_sin2comp = [0.19866392, 7.95872642, 2.61382166, 0.20077682, 0.02240027] # 2comp+sin
#[0.19869108, 7.958983, 2.61323087, 0.20085326, 0.02240435] # use not properly burn-in 2comp result

#radcliffe wave + free phi width
#p_rad1comp = ([0.2000, 11.3651, 3.8522, 4.7776, 0.0968, 31.7963, 9.6802], 
#	[-0.1466, 14.2452, 4.5423, -17.1834], [-0.0997, 10.7216, 4.2233, 61.9456], [-0.1033, 10.0907, 3.5992, 106.0447], [0.2852, 11.9379, 7.9774, 141.9115])
p_rad1comp = (\
	[0.1642051497224096, 11.414854787892217, 4.128583627294581, 5.225774753530958, 0.09565408147605053, 31.930404841538362, 9.517598869442107],\
	[-0.0752974320964544, 14.738579782282665, 6.270996498318865, -18.0477289625138],\
	[-0.0751491985230566, 10.638797575560321, 4.162115683827265, 61.598374168264655],\
	[-0.08233119781797832, 10.074165882338802, 3.8639715704036672, 104.8403903270003],\
	[0.0990740265830077, 12.023991841711359, 7.883899877704808, 144.88903338908486]\
	)

#p_rad2comp = ([0.1798, 11.3697, 6.8567, 4.5007, 0.0800, 31.6410, 9.2270],
#	[-0.1359, 14.1529, 4.3494,- 16.3676], [-0.1291, 10.7583, 4.4029, 61.6866], [-0.1204, 10.1934, 3.5576, 106.0012], [0.2537, 11.9445, 7.6762, 143.0137])
p_rad2comp = (\
	[0.13820738517301873, 11.347985842308768, 6.340582229618389, 4.4619858178452745, 0.06782841802800504, 31.94008898314093, 8.47732501926159],\
	[-0.07366799068508818, 14.280756002073709, 4.292107182706715, -18.104910091680452],\
	[-0.1122249790458754, 10.736459606294025, 4.4298851050408015, 61.092745362061486],\
	[-0.1065529170374176, 10.276317477103312, 3.802771325549306, 103.99392236131952],\
	[0.08087599862431931, 12.016131284829374, 8.057626049588588, 147.03947811305233]
	)
#p_circ1comp = [0.1673, 22.3272, 52.3849, 12.7471, 0.8703]
p_circ1comp = [0.12022708146730038, 21.403601891560587, 52.5683095317009, 12.70490222708815, 0.9438855120388356]
#[0.1180277698730274, 21.494414320521308, 52.33155122498118, 12.712047956550592, 0.93938384269732]
#p_circ2comp = [0.1618, 21.6950, 52.8728, 12.6713, 0.8555]
p_circ2comp = [0.11817596030526281, 20.49938444782939, 53.29284840274242, 12.550864375494976, 0.8821148522344949]
#[0.11392508987721547, 20.456378814518096, 53.03621039757925, 12.547505887478895, 0.8765020569282277]

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
#RC giant
colorRC = 'gray'
styleRC = (0,(5,2,1,2))#'--'
labelRC = 'Red clump'
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
ring_kws_rc   = dict(color=colorRC,   ls=styleRC,   lw=1,   label=labelRC,   zorder=13)
ring_kws_co1  = dict(color=colorCO1,  ls=styleCO1,  lw=0.8, label=labelCO1,  zorder=14)
ring_kws_co2  = dict(color=colorCO2,  ls=styleCO2,  lw=0.8, label=labelCO2,  zorder=15)
ring_kws_sin  = dict(color=colorCO1,  ls=(0, (1.5,1)),lw=3,   label='SIN comp',zorder=16)
ring_kws_text = dict(color=colorText, ha='left', va='top', fontsize=20, fontweight='bold')

### kws in radial plots: corrugation1b.py, corrugation_resi.py, corrugation_resi_xf.py
rad_kws_mc   = dict(color='#008080', edgecolors='none', alpha=0.3, zorder=5)
rad_kws_err  = dict(color=darkerColor('#008080'), fmt='.', markersize=5., elinewidth = 1, capsize=2, zorder=10)
rad_kws_hi   = dict(color=colorHI,   ls=styleHI,   lw=1.5, label=labelHI,   zorder=11)
rad_kws_ceph = dict(color=colorCeph, ls=styleCeph, lw=1.5, label=labelCeph, zorder=12)
rad_kws_rc   = dict(color=colorRC,   ls=styleRC,   lw=1.5, label=labelRC,   zorder=13)
rad_kws_co1  = dict(color=colorCO1,  ls=styleCO1,  lw=1,   label=labelCO1,  zorder=14)
rad_kws_co2  = dict(color=colorCO2,  ls=styleCO2,  lw=1,   label=labelCO2,  zorder=15)
rad_kws_sin  = dict(color=colorCO1,  ls=':',       lw=3,   label='SIN comp',zorder=11)
rad_kws_text = dict(color=colorText, ha='left', va='top', fontsize=20, fontweight='bold')
#rad_sep = [155, 143, 131, 119, 107, 95, 83, 71, 59, 47, 38, 29, 20, 11, -11, -17, -27]
rad_sep = [155, 143, 131, 119, 107, 95, 83, 71, 59, 47, 38, 29, 20, 10, -7, -16, -27]

### kws in arm plots: warp_out.py / warp_per.py / warp_osc.py / warp_osc_b.py
arm_kws_mc   = dict(color='#4169e1', edgecolors='none', alpha=0.2, zorder=5)
arm_kws_bin  = dict(color=darkerColor('#4169e1'),  lw=2, alpha=0.8, zorder=20)
arm_kws_hi   = dict(color=colorHI,   ls=styleHI,   lw=1.5, label=labelHI,   zorder=11)
arm_kws_ceph = dict(color=colorCeph, ls=styleCeph, lw=1.5, label=labelCeph, zorder=12)
arm_kws_rc   = dict(color=colorRC,   ls=styleRC,   lw=1.5, label=labelRC,   zorder=13)
arm_kws_co1  = dict(color=colorCO1,  ls=styleCO1,  lw=1,   label=labelCO1,  zorder=14)
arm_kws_co2  = dict(color=colorCO2,  ls=styleCO2,  lw=1,   label=labelCO2,  zorder=15)
arm_kws_text = dict(color=colorText, ha='left', va='top', fontsize=20, fontweight='bold')

textwidth=18 #inches, width of text in latex
subfigureIndexFont = dict(family="Arial Black", size=28)
mpl.rcParams['savefig.format'] = 'pdf'
mpl.rcParams['savefig.dpi'] = 280
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
	### p = (PHIkink, Rkink, pitch1, pitch2, ...)
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


def function_warp(x, p=p_1comp, rad=None, circ=None):
	### warp function, return Z
	R, PHI = x
	a0 = p[0]

	R, PHI = np.broadcast_arrays(R, PHI)
	#if isinstance(R, (int, float)):
	#	R = np.ones(PHI.shape) * R
	#elif isinstance(PHI, (int, float)):
	#	PHI = np.ones(R.shape) * PHI
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

	if rad is not None:
		# radcliffe wave function
		Amp, Rrad, sigmarad, Period, Amp0, PHIcirc, sigmacirc = rad[0]
		d = R-Rrad
		z += np.exp(-d**2/2/sigmarad**2) * (Amp * np.sin(2 * np.pi * d / Period) + Amp0) * np.exp(-(PHI-PHIcirc)**2/ 2 / sigmacirc**2)

		for i in [1,2,3]:
			Amp, Rrad, Period, PHIcirc = rad[i]
			d = R-Rrad
			z += np.exp(-d**2/2/4**2) * (Amp * np.sin(2 * np.pi * d / Period)) * np.exp(-(PHI-PHIcirc)**2/ 2 / 9.5**2)
		Amp, Rrad, Period, PHIcirc = rad[4]
		d = R-Rrad
		z += np.exp(-d**2/2/2**2) * (Amp * np.sin(2 * np.pi * d / Period)) * np.exp(-(PHI-PHIcirc)**2/ 2 / 9.5**2)
		

	if circ is not None:
		Amp, PHIcirc, Period, Rrad, sigmarad = circ
		z += Amp * np.sin(2 * np.pi * (PHI-PHIcirc) / Period ) * np.exp(-(R-Rrad)**2/ 2 / sigmarad**2)

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

	normSize = np.clip(np.sqrt(mass) / 3, 5, 120)
	return normSize, normColor


def cal_warp_HI(R = np.linspace(7.5,22,200), PHI = np.linspace(-30,170,200), scale=0.95):
	# regulate R, PHI *95%
	R, PHI = np.broadcast_arrays(R, PHI)
	#if isinstance(R, (int, float)):
	#	R = np.ones(PHI.shape) * R
	#else:
	#	PHI = np.ones(R.shape) * PHI

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
	R, PHI = np.broadcast_arrays(R, PHI)
	p1=[9.26, 17.4, 0.148]
	p2=[7.72, 17.5, 0.06, 1.33]

	zw1, zw2, zw4, zw5 =[], [], [], []
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


def cal_warp_Ceph(R = np.linspace(7.5,22,200), PHI = np.linspace(-30,170,200)):
	### calculate z at certain R and PHI in kpc, using ChenXiaodian's model
	#b=1
	R, PHI = np.broadcast_arrays(R, PHI)

	R_w, phi_w, Amp = [9.26, 17.4, 0.148]
	zw1 = np.zeros_like(R)
	idx = R>R_w
	zw1[idx] += Amp * (R[idx]-R_w) * np.sin(np.deg2rad(PHI[idx]-phi_w))

	#b=1.33
	R_w, phi_w, Amp, b_w=[7.72, 17.5, 0.06, 1.33]
	zw2 = np.zeros_like(R)
	idx = R>R_w
	zw2[idx] += Amp * (R[idx]-R_w)**b_w * np.sin(np.deg2rad(PHI[idx]-phi_w))

	return zw1, zw2


def cal_warp_MWSIP(R = np.linspace(7.5,22,200), PHI = np.linspace(-30,170,200), p=p_1comp):
	### calculate z at certain R and PHI in kpc
	R, PHI = np.broadcast_arrays(R, PHI)
	zw1=function_warp((R, PHI), p_1comp) ##mwisp first component
	zw2=function_warp((R, PHI), p_2comp) #mwisp two component
	return zw1, zw2


def cal_warp_RC(R = np.linspace(7.5,22,200), PHI = np.linspace(-30,170,200)):
	### calculate z at certain R and PHI in kpc, using Khanna2025 model
	R, PHI = np.broadcast_arrays(R, PHI)
	R_w = 8.79
	Amp = 0.083
	a_w = 2.0
	phi_w = 180-178

	zw = np.zeros_like(R)
	idx = R>R_w
	zw[idx] += Amp * (R[idx]-R_w)**a_w * np.sin(np.deg2rad(PHI[idx]-phi_w))
	return zw


def weighted_avg_and_std(values, weights):
	average = np.average(values, weights=weights)
	# Fast and numerically precise:
	variance = np.average((values-average)**2, weights=weights)
	return (average, np.sqrt(variance))


def load_data(arm='all', excluded=True):
	data = []
	if (arm=='all') | (arm=='out'):
		out_data = np.loadtxt('out_para.txt', comments='#')
		data.append(out_data)
	if (arm=='all') | (arm=='osc'):
		osc_data = np.loadtxt('osc2_para.txt', comments='#')
		data.append(osc_data)
	if (arm=='all') | (arm=='per'):
		per_data = np.loadtxt('per_para.txt', comments='#')
		data.append(per_data)

	data = np.vstack(data).T#(osc_data, out_data, per_data)).T

	l, b, v, _, d, \
	ed, sz, mass, complete, sd, \
	vir, rgal, x, y, z, \
	index, l_rms, b_rms, v_rms, az, \
	T_int, area, tpeak, tex, nh2, \
	pint, fill, lp, bp, vp = data

	if excluded:
		### exclude distant but massive clouds within 195~200
		idx = (mass > 0) & ~((l>195) & (l<199) & (x**2+y**2>16**2))# & (xx**2+yy**2<16**2)# & (mass>1e4))
	else:
		idx = (mass > 0)
	data = data[:,idx]

	l, b, v, lwidth, distance, \
	err_dist, size_pc, mass, complete, surface_density, \
	viral_para, r, x, y, z,\
	index, l_rms, b_rms, v_rms, az, \
	T_int, area, tpeak, tex, nh2, \
	pint, fill, lp, bp, vp = data

	cosl = np.cos(np.deg2rad(l))
	cosb = np.cos(np.deg2rad(b))
	sinb = np.sin(np.deg2rad(b))
	err2_Z = (err_dist * sinb)**2
	err2_Z[err2_Z < 1e-12] = 1e-12 #avoid dZ=0
	#err2_R = (err_dist * np.abs((distance * cosb**2 - Rsun * cosl * cosb) / r))**2

	# normalize mass
	norm_mass = np.sqrt(mass)
	norm_mass = norm_mass/np.mean(norm_mass)

	# weight
	H = 150 # kpc
	weight = norm_mass / (err2_Z + H**2)

	return r, az, x, y, z, l, b, mass, distance, err_dist, weight



def bin_data(x, y, w, grid=None, width=1, method='median', minbin=10):
	xmean = []
	xstd = []
	ymean = []
	ystd = []
	for xc in grid:
		idx = (x>=xc-width/2) & (x<xc+width/2)
		#if idx.sum()>0: print(xc, idx.sum())
		if idx.sum() >= minbin:
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
def g2gc(l, b, d, xyz=False):
	#galactic to galactocentric
	X = d * np.cos(b*d2r) * np.sin(l*d2r)
	Y = - d * np.cos(b*d2r) * np.cos(l*d2r) + 8.15
	Z = d * np.sin(b*d2r)
	if xyz: return X, Y, Z

	Rgal = np.sqrt(X**2 + Y**2)
	Az = np.arctan2(X, Y) / d2r #az=0 toward sun and counterclockwise
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
obDat = hdu[1].data
starL = obDat['GLON']
starB = obDat['GLAT']
starD = obDat['DistBJ']/1e3
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



### Poggio warp (Poggio, E., et al. 2025 A&A 699, 199)
def function_PoggioWarp(x, cepheids=True, straight=False):
	### warp function, return Z
	R, PHI = x

	if isinstance(R, (int, float)):
		R = np.ones(PHI.shape) * R
	elif isinstance(PHI, (int, float)):
		PHI = np.ones(R.shape) * PHI

	if cepheids:
		R_w = 7.7 #kpc
		h_w0 = 0.057 #kpc**(1-aw)
		a_w = 1.3
		phi_lon0 = 0.9 #deg
		R_t = 12.6 #kpc
		beta_w = 8.0 #kpc
		phi_straightlon = 9.9 #deg
	else:
		#young giants
		R_w = 5.5 #kpc
		h_w0 = 0.012 #kpc**(1-aw)
		a_w = 1.9
		phi_lon0 = 0.1 #deg
		R_t = 12.1 #kpc
		beta_w = 9.9 #kpc
		phi_straightlon = 14.0 #deg

	def h_w(R):
		cond = [R<=R_w, R>R_w]
		func = [lambda x: h_w0, lambda x: h_w0 * (R-R_w)**a_w]
		return np.piecewise(R, cond, func)
	def phi_lon(R, straight=False):
		if straight:
			return phi_lon0 * np.ones(R.shape)
		else:
			cond = [R<=R_t, R>R_t]
			func = [lambda x: phi_lon0, lambda x: phi_lon0 + beta_w * (x-R_t)]
			return np.piecewise(R, cond, func)

	Z_w = h_w(R) * np.sin(np.deg2rad(PHI - phi_lon(R, straight=straight)))
	return Z_w

#load Poggio Cepheids catalog
from astropy.io import ascii
cepheid = ascii.read("star/CepheidTable")['glon', 'glat', 'd-all', 'e_d-all']
cepheid.rename_columns(['glon', 'glat', 'd-all', 'e_d-all'], ['l', 'b', 'd', 'ed'])
for col in cepheid.columns.values():
	col.unit = None
cepheid['d'] = cepheid['d']/1e3 #to kpc
cepheid['ed'] = cepheid['ed']/1e3 #to kpc
cepheid = cepheid[cepheid['ed'] < cepheid['d']*0.1]
cepheid['w'] = 1/cepheid['ed']**2
cepheid['x'], cepheid['y'], cepheid['z'] = g2gc(cepheid['l'], cepheid['b'], cepheid['d'], xyz=True)
cepheid['r'], cepheid['az'], cepheid['z'] = g2gc(cepheid['l'], cepheid['b'], cepheid['d'], xyz=False)

#load Poggio Young giant catalog
ygiant = ascii.read("star/YoungGiantTable")['col2', 'col3', 'col4', 'col5']
ygiant.rename_columns(['col2', 'col3', 'col4', 'col5'], ['l', 'b', 'd', 'ed'])
ygiant = ygiant[ygiant['ed'] < ygiant['d']*0.5]
ygiant['w'] = 1/ygiant['ed']**2
ygiant['x'], ygiant['y'], ygiant['z'] = g2gc(ygiant['l'], ygiant['b'], ygiant['d'], xyz=True)
ygiant['r'], ygiant['az'], ygiant['z'] = g2gc(ygiant['l'], ygiant['b'], ygiant['d'], xyz=False)


def arm_mask(az, r):
	### interpolate a mask for outer arm
	from scipy.stats import binned_statistic
	def envelope(arm_file, bins = np.arange(-40, 180, 3)):
		data = np.loadtxt(arm_file, comments='#')
		r = data[:,11]
		#z = data[:,14]
		az = data[:,19]
		ma, edges, num = binned_statistic(az, r, statistic='max', bins=bins)
		mi, edges, num = binned_statistic(az, r, statistic='min', bins=bins)
		edges = (edges[:-1]+edges[1:])/2
		idx = ((edges>-26) & (edges<-6)) | ((edges>4) & (edges<153))
		return edges[idx], ma[idx], mi[idx]

	e, perMax, _ = envelope('per_para.txt')
	e, outMax, outMin = envelope('out_para.txt')
	e, _, oscMin = envelope('osc2_para.txt')

	import scipy.interpolate as sp_interp
	outMin = np.nanmean([perMax, outMin], axis=0)
	outMax = np.nanmean([outMax, oscMin], axis=0) 

	fMin = sp_interp.interp1d(e, outMin, fill_value='extrapolate')
	fMax = sp_interp.interp1d(e, outMax, fill_value='extrapolate')

	mask = (r >= fMin(az)) & (r <= fMax(az))
	return mask

'''
def readCepheidCat(cat='allGalCep.listID'):
	import pandas as pd
	df = pd.read_csv(cat, comment='#', sep='\s+')
	l = df[]
	return df

'''