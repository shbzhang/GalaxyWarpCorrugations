### params shared by other scripts

import numpy as np

### Warp Params
p_1comp = [0,  0.09362826, 8.56849936, 1.05015227, -0.7049789 ]#
#[-0.0717748888594908, 0.1074529579441738, 7.8409214840721795, 0.9370789056657302, -8.059413168263955]   ##1st exclude 160-200
p_2comp = [0, 1.19210237e-01, 9.01434283e+00, -3.53962133e+00, -1.24379884e-02, 1.27166892e+01, 2.03007431e+00, -2.04479290e+01]
#[-0.08225466502024195, 0.10192102165531088, 7.779011564372266, 1.0, -11.000005387066379, -0.11935558287557957, 13.791351404542997, 1.0, -23.76326592401999]##a0, a1, Rw1, bw1, PHIw1, a2, Rw2, bw2, PHIw2 exclude 160-200

p_sin1comp = [0.2012454, 7.8246938, 2.78618375, 0.16944783, 0.01766535] # 1comp+sin

p_sin2comp = [0.19869108, 7.958983, 2.61323087, 0.20085326, 0.02240435] # 2comp+sin

### Arm Params
best_per=[30.03, 10.062122711773025, 9.839300621559302]  #mass  
best_out=[20.25195158643863, 13.259894850233614, 11.05257961109539, 3.500917359458265]  ##mass
best_osc=[47, 16.190233947534733, 12.330164850703794] ##47, 16.060484876948312, 11.996988806408883] #mass


### color and style
col_mc = '#008080'#'#1f8097'
col_co = '#f33a23'
sty_co1 = '--'
sty_co2 = '-'
col_hi = '#df8f0e'
sty_hi = ':'
col_ceph = '#8f00ff'#'#751fdf'#
sty_ceph = ':'
col_text = '#1e3f66'
col_err = '#A9A9A9'

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


def cal_warp(R = np.linspace(7.5,22,200), PHI = np.linspace(-30,170,200)):
	# regulate R, PHI
	if isinstance(R, (int, float)):
		R = np.ones(PHI.shape) * R
	else:
		PHI = np.ones(R.shape) * PHI

	#   R=R+0.5    # considering the different rotaion curve
	w1p = [9, 197, 10, -3.1] # [k0,k1,rk1,k2]
	w1 = w1p[0] + w1p[1]*(R-w1p[2]) + w1p[3]*(R-w1p[2])**2
	w0p = [-66, 150, 15, -0.47] # [k0,k1,rk1,k2]
	w0 = w0p[0] + w0p[1]*(R-w0p[2]) + w0p[3]*(R-w0p[2])**2
	w2p = [-70,171,15,-5.3] # [k0,k1,rk1,k2]
	w2 = w2p[0] + w2p[1]*(R-w2p[2]) + w2p[3]*(R-w2p[2])**2

	w = w1*np.sin(np.deg2rad(PHI))

	idx = (w2 >= 150) & (R > 15)
	w[idx] = (w0 + w1*np.sin(np.deg2rad(PHI)) + w2*np.sin(np.deg2rad(2.*PHI)))[idx]

	idx = (w2 < 150) & (R > 15)
	w[idx] = (w0+w1*np.sin(np.deg2rad(PHI)))[idx]
	'''
	if (w2 >= 150) & (R > 15):
		w = w0 + w1*np.sin(np.deg2rad(PHI)) + w2*np.sin(np.deg2rad(2.*PHI))
	elif (w2 < 150) & (R > 15):
		w = w0+w1*np.sin(np.deg2rad(PHI))
	else:
		w = w1*np.sin(np.deg2rad(PHI))
	'''
	return w/1000.


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


def cal_zcen_zrms(az, *zz, weights=1, binsize = 6, nbin=60, bin0=160):
	bin_az, bin_zz= [], []
	for i in range(nbin):
		PHI1 = bin0 - i * binsize
		PHI2 = bin0 - (i+1) * binsize
		idx = (az >=min(PHI1, PHI2)) & (az < max(PHI1, PHI2))
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


def darker_hex(hex_color, factor=0.4):
	"""
	将一个十六进制颜色变暗。
	
	:param hex_color: 十六进制颜色字符串，可以带或不带 '#' 前缀
	:param factor: 减少比例，默认是0.9，即减少10%
	:return: 变暗后的十六进制颜色字符串
	"""
	# 如果有 '#' 前缀，去掉它
	if hex_color.startswith('#'):
		hex_color = hex_color[1:]
	
	# 将hex转换为RGB元组
	r = int(hex_color[0:2], 16)
	g = int(hex_color[2:4], 16)
	b = int(hex_color[4:6], 16)
	
	# 计算新的RGB值
	r_dark = int(r * factor)
	g_dark = int(g * factor)
	b_dark = int(b * factor)
	
	# 确保RGB值在0到255之间
	r_dark = max(0, min(255, r_dark))
	g_dark = max(0, min(255, g_dark))
	b_dark = max(0, min(255, b_dark))
	
	# 将RGB转换回hex
	return f'#{r_dark:02x}{g_dark:02x}{b_dark:02x}'.upper()
