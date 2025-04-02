### get l, b, v, d from Rgal, Az, Z

from astropy.io import fits
import numpy as np
from Distance import KinUnit, Distance
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt

#ku = KinUnit( Ro = 8.15, To = 236, dTdR = -0.2, Uo = 10.6, Vo = 10.7, Wo = 7.6, Us = 6.0, Vs = -4.3, Ws = -3.0, err_v_lsr = 5.0) #
d2r = np.pi/180


def getKD(l, b, v):
	skyc = SkyCoord(l = l, b = b, frame = 'galactic', unit = 'deg')
	d = Distance(skyc, v, np.zeros(l.size)).Kdist()
	return d[1]

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
	fig,ax = plt.subplots(3,2, figsize=(10,8))
	best_out=[20.25195158643863, 13.259894850233614, 11.05257961109539, 3.500917359458265]  ##mass
	Az = np.linspace(-25, 155, 300)
	Rgal = function(Az, best_out)
	out_zh,out_z1,out_z2,out_z4,out_z5 = [],[],[],[],[]
	for i in range(len(Rgal)):
		z1,z2,z4,z5 = cal_warpc(Rgal[i],Az[i])
		zh = cal_warp(Rgal[i],Az[i])
		out_zh.append(zh)
		out_z1.append(z1)
		out_z2.append(z2)
		out_z4.append(z4)
		out_z5.append(z5)

#	Az = np.linspace(160, -30, 100)
#	Rgal = 21.275633909457124*np.exp(-0.00537711182911569*Az)
#	Z = np.sin(Az * d2r * 5)
	

	ax[0,0].plot(Az, Rgal, 'b.')
	ax[0,0].set_xlabel('Az')
	ax[0,0].set_ylabel('Rgal')

	ax[1,0].plot(Az, out_z5, 'b.')
	ax[1,0].set_xlabel('Az')
	ax[1,0].set_ylabel('Z')


	l, b, v, d = inverseKD(np.array(Rgal), np.array(Az), np.array(out_z5))

	ax[0,1].plot(l, d, 'r.')	
	ax[0,1].set_xlabel('l')
	ax[0,1].set_ylabel('distance')
	ax[1,1].plot(l, b, 'r.')
	ax[1,1].set_xlabel('l')
	ax[1,1].set_ylabel('b')

	ax[2,1].plot(l, v, 'r.')
	ax[2,1].set_xlabel('l')
	ax[2,1].set_ylabel('v')

	plt.show()
