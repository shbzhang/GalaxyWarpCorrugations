### plot Az-(Z-warp model) of MCs in different Rgal


import numpy as np
import pylab as pl
import scipy.interpolate as sp_interp
from scipy.stats import norm
import scipy.optimize as sp_opt
from shared import *
import matplotlib as mpl
import matplotlib.pyplot as plt



comp=1

if __name__ == '__main__':
	az_min, az_max = -30, 165

	out_data = np.loadtxt('out_para.txt',comments='#')
	osc_data = np.loadtxt('osc2_para.txt',comments='#')
	per_data = np.loadtxt('per_para.txt',comments='#')
	out_bb = out_data[:,1]
	out_ll = out_data[:,0]
	osc_bb = osc_data[:,1]
	osc_ll = osc_data[:,0]
	per_bb = per_data[:,1]
	per_ll = per_data[:,0]


	data = np.vstack((osc_data, out_data, per_data)).T
	az = data[19]
	rr = data[11]
	zz = data[14] - function_warp((rr, az), p=p_1comp if comp==1 else p_2comp)
	vv = data[2]
	mass = data[7]
	com = data[8]
	scatter_size, scatter_color = mass_distribute(mass)

	theta_axis = np.linspace(-30,170,200)

	fig1, ax = plt.subplots(nrows=2, ncols=5, sharex=False, sharey=True, figsize=[17,7])
	plt.subplots_adjust(wspace=0, hspace=0.2)
	ax = ax.ravel()

	gr_sep = [8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 18.5, 20.5]
	for i in range(len(gr_sep)-1):
		gr1, gr2 = gr_sep[i:i+2]
		gr = (gr1 + gr2) / 2
		radius = 2*np.pi*gr*(np.abs(az_min-az_max)/360.)

		### filter gr
		idx = (rr>=gr1) & (rr<gr2)

		### plot scatter
		s = ax[i].scatter(az[idx], zz[idx], s=scatter_size[idx], c=col_mc, alpha=0.3, edgecolors='none', zorder=1)

		### plot binned average
		azcen, azrms, zcen, zrms, vcen, vrms, rcen, rrms = cal_zcen_zrms(az[idx], zz[idx], vv[idx], rr[idx], weights=mass[idx], binsize=6, nbin=60)

		ax[i].plot(azcen, zcen, '-', color=darker_hex(col_mc), lw=2, alpha=0.8, zorder=10)
		### fill gap
		idx = np.isfinite(zcen)
		ax[i].plot(azcen[idx], zcen[idx], '--', c=darker_hex(col_mc), lw=2, alpha=0.8, zorder=10)

		### plot errorbar
		ax[i].errorbar(azcen[idx], zcen[idx], yerr=zrms[idx]*2.355/2., fmt='^', c='#777777', markersize=0.05, elinewidth=1, capsize=1.1, zorder=8)

		### plot H line
		ax[i].plot([az_min, az_max], [0, 0], sty_co1 if comp==1 else sty_co2, color=col_co, lw=2, zorder=0)

		'''
		### plot warp models
		# HI
		if i>0:
			w = cal_warp(gr, theta_axis)
			ax[i].plot(theta_axis, w, ls=sty_hi, lw=1.5, color=col_hi, label='HI')
		# cepheids, co
		zw1, zw2, zw4, zw5=cal_warpc(gr, theta_axis)
		ax[i].plot(theta_axis, zw1, ls=sty_ceph, lw=1.5, color=col_ceph, label='Cepheids')
		ax[i].plot(theta_axis, zw4, ls=sty_co1, lw=1.5, color=col_co, label='CO 1comp')
		ax[i].plot(theta_axis, zw5, ls=sty_co2, lw=1.5, color=col_co, label='CO 2comp')
		'''

		### plot gr text
		text = '%i' % gr if gr%1==0 else '%.1f' % gr
		if i==0: text += ' kpc'
		ax[i].text(0.02, 0.95, text, ha='left', va='top', fontsize=11, fontweight='bold', color=col_text, transform=ax[i].transAxes)

		### axes
		ax[i].set_xticks(np.arange(-50, 200, 50))
		ax[i].set_yticks(np.arange(-3, 3, 0.5))
		ax[i].minorticks_on()
		#ax[i].tick_params(axis='both', which='major', labelsize=10.5, width=1.8)
		ax[i].set_xlim(az_min, az_max)
		ax[i].set_ylim(-0.8,0.8)
		#ax[i].set_yticklabels(fontsize=10.5, fontweight=1.8)
		upper = ax[i].twiny()
		upper.set_xticks(np.arange(0, 100, 20))
		upper.set_xlim(0, radius)
		if i == 5:
			ax[i].set_xticklabels(['%i$^{\circ}$' % i for i in np.arange(-50, 200, 50)])
			ax[i].set_xlabel('Galactocentric Azimuth',fontsize=13,fontweight='bold')
			ax[i].set_ylabel('Z (kpc)', fontsize=13, fontweight='bold')  
		else:
			ax[i].set_xticklabels([])

		if i==3: upper.set_xticklabels(['0', '20', '40   ', None, None]) # avoid overlap
		if i==4: upper.set_xticklabels(['  0', '20', '40 kpc', None, None]) # add kpc at the end
		if i==9: upper.set_xticklabels(['0', '20', '40', '60 kpc', None]) # add kpc at the end
		#if i==9: ax[i].legend()

	fig1.savefig('fig/az_z_warp_resi_%icomp.png' % comp,format='png',bbox_inches='tight', dpi=400)
	plt.show()


'''
import numpy as np
import pylab as pl
import scipy.interpolate as sp_interp
from scipy.stats import norm
import scipy.optimize as sp_opt
from shared import *

def fp2p(free_params):
	#free params -> mixed params with free and fixed
	
	#return params[r'$a_0$'][0], params[r'$a_1$'][0], params[r'$R_{w0}$'][0], params[r'$b_{w0}$'], params[r'$\phi_{w0}$'], \
	#	free_params
	
	i = 0
	p = []
	for v in params.values():
		if v[1]=='free':
			p.append(free_params[i])
			i += 1
		else:
			p.append(v[0])
	return p
def function(x, free_params, warp=True, sin=True):
	a0, a1, Rw1, bw1, PHIw1, a2, Rw2, bw2, PHIw2, a3, Rsin, Period0, Period1, PHIsin, PHIsinwid, a4 = fp2p(free_params)
	#a0, a1, Rw, b, fR, Rs, fPHI, PHIs = fp2p(free_params)
	R, PHI = x
	
	#y = np.zeros(R.size)
	y = np.ones(R.shape)*a0

	if warp:
		###warp component 1
		index = R>Rw1
		y[index] += a1 * (R[index]-Rw1)**bw1 * np.sin((PHI[index]-PHIw1)/180*np.pi)

		###warp component 2
		index = R>Rw2
		y[index] += a2 * (R[index]-Rw2)**bw2 * np.sin(2*(PHI[index]-PHIw2)/180*np.pi)

	if sin:
		###sin component
		index = R>Rsin
		y[index] += ( a4*(R[index] - Rsin) + a3*np.sin( (R[index] - Rsin) / (Period0 + Period1*(R[index]-Rsin)) * 2*np.pi ) ) * np.exp(-(PHI[index]-PHIsin)**2/ 2 / PHIsinwid**2)

	return y
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
			colors.append(col_cen)
		elif mass[i] >= 10000 and mass[i] < 100000:
			mass_dist[i] = 10
			colors.append('red')
		else:
			mass_dist[i] = 20
			colors.append('cyan')
	return mass_dist * 10,colors
def cal_warp(rg):
#	rg=rg+0.5    # considering the different rotaion curve
	w1p=[9,197,10,-3.1] # [k0,k1,rk1,k2]
	w1=w1p[0]+w1p[1]*(rg-w1p[2])+w1p[3]*(rg-w1p[2])**2
	w0p=[-66,150,15,-0.47] # [k0,k1,rk1,k2]
	w0=w0p[0]+w0p[1]*(rg-w0p[2])+w0p[3]*(rg-w0p[2])**2
	w2p=[-70,171,15,-5.3] # [k0,k1,rk1,k2]
	w2=w2p[0]+w2p[1]*(rg-w2p[2])+w2p[3]*(rg-w2p[2])**2
	thea, w = [], []
	thea = np.linspace(-30,170,200)
	if ((w2 >= 150) & (rg > 15)):
		w= w0+w1*np.sin(np.deg2rad(thea))+w2*np.sin(np.deg2rad(2.*thea))
	elif ((w2 < 150) & (rg > 15)):
		w= w0+w1*np.sin(np.deg2rad(thea))
	else:
		w= w1*np.sin(np.deg2rad(thea))
	return(thea,w/1000.)
def cal_warpc(rg):
	p1=[9.26,17.4,0.148]
	p2=[7.72,17.5,0.06,1.33]
	p3=[-2.97397887e-02,1.72161063e-02,-1.66075673e+01,7.93265614e-03,6.05519348e+01]
#	p4=[-0.06307575167323597, 0.09669646781880904, 7.798302460855644, 1.0, -7.792941115015663]   ###first component
	p4=[-0.07623866250085598, 0.10786410974035177, 7.799981504630573, 0.9349929777917091, -9.018240871527818] ##1st notcorrect
	p5=[-0.07,0.104,7.91,1.0,79.7,-0.114,13.7,1.0,20.7,0.1356797] #a0, a1, Rw1, bw1, PHIw1, a2, Rw2, bw2, PHIw2, a3
	thea, zw1, zw2, zw3,zw5 =[], [], [], [], []
	thea = np.linspace(-30,170,200)
	zw1=p1[2]*(rg-p1[0])*np.sin(np.deg2rad(thea-p1[1]))
	zw2=p2[2]*((rg-p2[0])**p2[3])*np.sin(np.deg2rad(thea-p2[1]))
	zw2=-p3[0]+(rg - 8)**2*(p3[1]*np.sin(np.deg2rad(thea-p3[2]))+p3[3]*np.sin(2.0*(thea-p3[4])/180.*np.pi))  #mwisp
	zw3=p4[0]+p4[1] * (rg-p4[2])**p4[3] * np.sin((thea-p4[4])/180*np.pi) #mwisp first component
	zw5=np.zeros(thea.size)
	if (rg>p5[6]):
		zw5 = p5[0]+p5[1] * (rg-p5[2])**p5[3] * np.sin((thea-p5[4])/180*np.pi)+p5[5] * (rg-p5[6])**p5[7] * np.cos(2*(thea-p5[8])/180*np.pi)#
	elif (rg<p5[6]):
		zw5 = p5[0]+p5[1] * (rg-p5[2])**p5[3] * np.sin((thea-p5[4])/180*np.pi)  
	return(thea,zw1,zw2,zw3,zw5)
def cal_warpmwisp(rg):
	a0, a, Rw, PHIw, b=[1.01,2.62,8.8,120.38,1.34]
	zw=a0 + a * (rg - Rw)**b * np.sin((thea-PHIw) / 180 * np.pi)
	return(thea,zw)
def weighted_avg_and_std(values, weights):
	average = np.average(values, weights=weights)
	# Fast and numerically precise:
	variance = np.average((values-average)**2, weights=weights)
	return (average, np.sqrt(variance))
def cal_zcen_zrms(az,zz,vv,rr,mass,bi):
	zcen,zrms,vcen,vrms,rcen,rrms,azcen,azrms,z_err,v_err,r_err,az_err,zcen_err,vcen_err,rcen_err,azcen_err = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
	for i in range(42):
		thea1 = 160-i*bi
		thea2 = 160-(i+1)*bi
		ind = np.where((az>thea2) & (az<thea1))
#        ran = np.abs(np.max(az[ind])-np.min(az[ind]))
#        if ((len(np.transpose(ind)) > 1) & (np.abs(np.max(az[ind])-np.min(az[ind])>= 3))):
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
 #           thea.append((thea2+bi/2.))
			z_err.append(temp_zrms/np.sqrt(len(np.transpose(ind))))
			v_err.append(temp_vrms/np.sqrt(len(np.transpose(ind))))
			r_err.append(temp_rrms/np.sqrt(len(np.transpose(ind))))
			az_err.append(temp_azrms/np.sqrt(len(np.transpose(ind))))
			zcen_err.append(np.abs(temp_zcen)/np.sqrt(len(np.transpose(ind))))
			vcen_err.append(temp_vcen/np.sqrt(len(np.transpose(ind))))
			rcen_err.append(temp_rcen/np.sqrt(len(np.transpose(ind))))
			azcen_err.append(temp_azcen/np.sqrt(len(np.transpose(ind))))
	return(zcen,zrms,vcen,vrms,rcen,rrms,azcen,azrms,z_err,v_err,r_err,az_err,zcen_err,vcen_err,rcen_err,azcen_err)
if __name__ == '__main__':
	col,col1,col_hi,col_cepheid,col_co,col_cen,col_err1,col_err = ['teal','none','olive','blue','crimson','black','#5A5A5A','#A9A9A9']
	best =  [-0.07616907102729832, 0.1078373858220152, 7.799921542773368, 0.9348508994894277, -9.02005784678019, 0.0, 0.0, 0.0, 0.0, 0.0, 8.3, 1.695, 0.372, 31.22, 9.27, 0.0]
	out_data = np.loadtxt('out_para.txt',comments='#')
	osc_data = np.loadtxt('osc2_para.txt',comments='#')
	per_data = np.loadtxt('per_para.txt',comments='#')
	#residual = np.loadtxt('residual.razm')
	out_mass = out_data[:,7]/out_data[:,8]
	out_rr = out_data[:,11]
	out_zz = out_data[:,14]#*1000.
	out_az = out_data[:,19]
	out_vv = out_data[:,2]
	out_bb = out_data[:,1]
	out_ll = out_data[:,0]
	osc_mass = osc_data[:,7]/osc_data[:,8]
	osc_rr = osc_data[:,11]
	osc_zz = osc_data[:,14]#*1000.
	osc_az = osc_data[:,19]
	osc_vv = osc_data[:,2]
	osc_bb = osc_data[:,1]
	osc_ll = osc_data[:,0]
	per_mass = per_data[:,7]/per_data[:,8]
	per_rr = per_data[:,11]
	per_zz = per_data[:,14]#*1000.
	per_az = per_data[:,19]
	per_vv = per_data[:,2]
	per_bb = per_data[:,1]
	per_ll = per_data[:,0]
	az_min, az_max = -30, 165
	#az = np.loadtxt('az_out1.txt')   ##9 to 19kpc
	
	az = np.hstack((osc_az,out_az,per_az))
	rr = np.hstack((osc_rr,out_rr,per_rr))
	zz = np.hstack((osc_zz,out_zz,per_zz))
	zz = zz - function_warp((rr, az))  #
	vv = np.hstack((osc_vv,out_vv,per_vv))
	mass = np.hstack((osc_mass,out_mass,per_mass))
	all_mass = np.hstack((osc_mass*1, out_mass*1., per_mass))    #####################################pl mass use different scale
	all_dist,all_colors = mass_distribute(all_mass)
	bii=7
	fig1, ax = pl.subplots(nrows=3,ncols=4,sharex=True, figsize=[17,7])
	pl.subplots_adjust(wspace=0,hspace=0.2)
	pl.subplot(2,5,1)
	gr=9
	radius = 2*np.pi*gr*(np.abs(az_min-az_max)/360.)
	thea,w = cal_warp(gr)
	thea,zw1,zw2,zw3,zw5=cal_warpc(gr)
	ind = np.where((rr>=(gr-0.5)) & (rr<(gr+0.5)))[0]
	zcen_9,zrms_9,vcen_9,vrms_9,rcen_9,rrms_9,azcen_9,azrms_9,zerr_9,verr_9,rerr_9,azerr_9,zcerr_9,vcerr_9,rcerr_9,azcerr_9=cal_zcen_zrms(az[ind],zz[ind],vv[ind],rr[ind],mass[ind],bii)
	az_9, zz_9, mass_9, zzz_9, mm_9, rr_9 = [], [], [], [], [], []
	az_9 = az[ind]
	mass_9 = mass[ind]
	zz_9 = zz[ind]
	rr_9 = rr[ind]
	maxval = int(np.max(azcen_9))
	minval = int(np.min(azcen_9))
	thea_9 = np.linspace(minval,maxval,maxval*25+1,endpoint=True)#np.arange(maxval+1)
	f = sp_interp.interp1d(azcen_9, zcen_9, fill_value='extrapolate')
	zcen_99 = f(thea_9)
	for i in range(len(ind)):
		ind2 = np.where(np.abs(az_9[i]-thea_9)<=0.02383)[0]
		if (len(ind2>=1)):
			zzz_9.append(zz_9[i]-zcen_99[ind2[0]])
			mm_9.append(mass_9[i])
			wei_9 = mm_9/np.linalg.norm(mm_9)

	pl.scatter(az[ind],zz[ind],s=all_dist[ind],c=col,alpha=0.3,edgecolors=col1)
	pl.plot([165,-30],[0,0],'--')
	pl.errorbar(azcen_9,zcen_9,yerr=np.array(zrms_9)*2.355/2.,fmt='^',c=col_err,markersize=0.05,elinewidth = 1, capsize=1.1)
	pl.plot(thea_9,zcen_99,'-',color=col_cen)
	pl.text(25,-0.45,'%skpc' %(gr),fontsize=11,fontweight='bold',color='b')
	pl.xlim([-30,165])
	pl.ylim([-0.61,0.61])
	pl.yticks([-0.6,-0.4,-0.2,0,0.2,0.4,0.6],['-0.6','-0.4','-0.2','0','0.2','0.4','0.6'],fontsize=11,fontweight=1.8)
	pl.xticks([0,50,100,150],['','','',''],fontsize=11,fontweight=1.8)
	pl.twiny()
	pl.xlim(0,radius)
	pl.xticks([0,20],['0','20'],fontsize=10.5,fontweight=1.8)
	################################
	pl.subplot(2,5,2)
	gr=10
	radius = 2*np.pi*gr*(np.abs(az_min-az_max)/360.)
	thea,w = cal_warp(gr)
	thea,zw1,zw2,zw3,zw5=cal_warpc(gr)
	ind = np.where((rr>=(gr-0.5)) & (rr<(gr+0.5)))[0]
	zcen_0,zrms_0,vcen_0,vrms_0,rcen_0,rrms_0,azcen_0,azrms_0,zerr_0,verr_0,rerr_0,azerr_0,zcerr_0,vcerr_0,rcerr_0,azcerr_0=cal_zcen_zrms(az[ind],zz[ind],vv[ind],rr[ind],mass[ind],bii)
	az_0, zz_0, mass_0, zzz_10, mm_10, rr_10 = [], [], [], [], [], []
	az_0 = az[ind]
	mass_0 = mass[ind]
	zz_0 = zz[ind]
	rr_10 =rr[ind]
	maxval = int(np.max(azcen_0))
	minval = int(np.min(azcen_0))
	thea_10 = np.linspace(minval,maxval,maxval*25+1,endpoint=True)#np.arange(maxval+1)
	f = sp_interp.interp1d(azcen_0, zcen_0, fill_value='extrapolate')
	zcen_10 = f(thea_10)
	for i in range(len(ind)):
		ind2 = np.where(np.abs(az_0[i]-thea_10)<=0.02383)[0]
		if (len(ind2>=1)):
			zzz_10.append(zz_0[i]-zcen_10[ind2[0]])
			mm_10.append(mass_0[i])
			wei_10 = mm_10/np.linalg.norm(mm_10)

	pl.scatter(az[ind],zz[ind],s=all_dist[ind],c=col,alpha=0.3,edgecolors=col1)
	pl.plot([170,-40],[0,0],'--')
	pl.errorbar(azcen_0,zcen_0,yerr=np.array(zrms_0)*2.355/2.,fmt='^',c=col_err,markersize=0.05,elinewidth = 1, capsize=1.1)
	pl.plot(thea_10,zcen_10,'-',color=col_cen)
	pl.text(25,-0.45,'%skpc' %(gr),fontsize=11,fontweight='bold',color='b')
	pl.xlim([165,-30])
	pl.ylim([-0.61,0.61])
	pl.yticks([-0.6,-0.4,-0.2,0,0.2,0.4,0.6],['','','','','','',''],fontsize=11,fontweight='bold')
	pl.xticks([0,50,100,150],['','','',''],fontsize=11,fontweight=1.8)
	pl.twiny()
	pl.xlim(0,radius)
	pl.xticks([0,20],['0','20'],fontsize=10.5,fontweight=1.8)    
	gr=11
	pl.subplot(2,5,3)
	radius = 2*np.pi*gr*(np.abs(az_min-az_max)/360.)
	thea,w = cal_warp(gr)
	thea,zw1,zw2,zw3,zw5=cal_warpc(gr)
	ind = np.where((rr>=(gr-0.5)) & (rr<(gr+0.5)))[0]
	zcen_1,zrms_1,vcen_1,vrms_1,rcen_1,rrms_1,azcen_1,azrms_1,zerr_1,verr_1,rerr_1,azerr_1,zcerr_1,vcerr_1,rcerr_1,azcerr_1=cal_zcen_zrms(az[ind],zz[ind],vv[ind],rr[ind],mass[ind],bii)
	az_1, zz_1, mass_1, zzz_11, mm_11, rr_11 = [], [], [], [], [], []
	az_1 = az[ind]
	mass_1 = mass[ind]
	zz_1 = zz[ind]
	rr_11 =rr[ind]
	maxval = int(np.max(azcen_1))
	minval = int(np.min(azcen_1))
	thea_11 = np.linspace(minval,maxval,maxval*25+1,endpoint=True)#np.arange(maxval+1)
	f = sp_interp.interp1d(azcen_1, zcen_1, fill_value='extrapolate')
	zcen_11 = f(thea_11)
	for i in range(len(ind)):
		ind2 = np.where(np.abs(az_1[i]-thea_11)<=0.0216)[0]
		if (len(ind2>=1)):
			zzz_11.append(zz_1[i]-zcen_11[ind2[0]])
			mm_11.append(mass_1[i])
			wei_11 = mm_11/np.linalg.norm(mm_11)

	pl.scatter(az[ind],zz[ind],s=all_dist[ind],c=col,alpha=0.3,edgecolors=col1)
	pl.plot([170,-40],[0,0],'--')
	pl.errorbar(azcen_1,zcen_1,yerr=np.array(zrms_1)*2.355/2.,fmt='^',c=col_err,markersize=0.05,elinewidth = 1, capsize=1.1)
	pl.plot(thea_11,zcen_11,'-',color=col_cen)
	pl.text(25,-0.45,'%skpc' %(gr),fontsize=11,fontweight='bold',color='b')
	pl.xlim([-30,165])
	pl.ylim([-0.61,0.61])
	pl.yticks([-0.6,-0.4,-0.2,0,0.2,0.4,0.6],['','','','','','',''],fontsize=11,fontweight='bold')
	pl.xticks([0,50,100,150],['','','',''],fontsize=11,fontweight=1.8)
	pl.twiny()
	pl.xlim(0,radius)
	pl.xticks([0,20],['0','20'],fontsize=10.5,fontweight=1.8)    
	gr=12
	pl.subplot(2,5,4)
	radius = 2*np.pi*gr*(np.abs(az_min-az_max)/360.)
	thea,w = cal_warp(gr)
	thea,zw1,zw2,zw3,zw5=cal_warpc(gr)
	ind = np.where((rr>=(gr-0.5)) & (rr<(gr+0.5)))[0]
	zcen_2,zrms_2,vcen_2,vrms_2,rcen_2,rrms_2,azcen_2,azrms_2,zerr_2,verr_2,rerr_2,azerr_2,zcerr_2,vcerr_2,rcerr_2,az_cerr_2=cal_zcen_zrms(az[ind],zz[ind],vv[ind],rr[ind],mass[ind],bii)
	az_2, zz_2, mass_2, zzz_12, mm_12, rr_12 = [], [], [], [], [], []
	az_2 = az[ind]
	mass_2 = mass[ind]
	zz_2 = zz[ind]
	rr_12 = rr[ind]
	maxval = int(np.max(azcen_2))
	minval = int(np.min(azcen_2))
	thea_12 = np.linspace(minval,maxval,maxval*25+1,endpoint=True)#np.arange(maxval+1)
	f = sp_interp.interp1d(azcen_2, zcen_2, fill_value='extrapolate')
	zcen_12 = f(thea_12)
#    print(thea_12)
	for i in range(len(ind)):
		ind2 = np.where(np.abs(az_2[i]-thea_12)<=0.0222)[0]
		if (len(ind2>=1)):
			zzz_12.append(zz_2[i]-zcen_12[ind2[0]])
			mm_12.append(mass_2[i])
			wei_12 = mm_12/np.linalg.norm(mm_12)

	pl.scatter(az[ind],zz[ind],s=all_dist[ind],c=col,alpha=0.3,edgecolors=col1)
	pl.plot([170,-40],[0,0],'--')
	pl.errorbar(azcen_2,zcen_2,yerr=np.array(zrms_2)*2.355/2.,fmt='^',c=col_err,markersize=0.05,elinewidth = 1, capsize=1.1)
	pl.plot(thea_12,zcen_12,'-',color=col_cen)
	pl.text(25,-0.45,'%skpc' %(gr),fontsize=11,fontweight='bold',color='b')
	pl.xlim([-30,165])
	pl.ylim([-0.61,0.61])
	pl.yticks([-0.6,-0.4,-0.2,0,0.2,0.4,0.6],['','','','','','',''],fontsize=11,fontweight='bold')
	pl.xticks([0,50,100,150],['','','',''],fontsize=11,fontweight=1.8)
	pl.twiny()
	pl.xlim(0,radius)
	pl.xticks([0,20,40],['0','20','40'],fontsize=10.5,fontweight=1.8)    
	gr=13
	pl.subplot(2,5,5)
	radius = 2*np.pi*gr*(np.abs(az_min-az_max)/360.)
	thea,w = cal_warp(gr)
	thea,zw1,zw2,zw3,zw5=cal_warpc(gr)
	ind = np.where((rr>=(gr-0.5)) & (rr<(gr+0.5)))[0]
	zcen_3,zrms_3,vcen_3,vrms_3,rcen_3,rrms_3,azcen_3,azrms_3,zerr_3,verr_3,rerr_3,azerr_3,zcerr_3,vcerr_3,rcerr_3,azcerr_3=cal_zcen_zrms(az[ind],zz[ind],vv[ind],rr[ind],mass[ind],bii)
	az_3, zz_3, mass_3, zzz_13, mm_13, rr_13 = [], [], [], [], [], []
	az_3 = az[ind]
	mass_3 = mass[ind]
	zz_3 = zz[ind]
	rr_13 = rr[ind]
	maxval = int(np.max(azcen_3))
	minval = int(np.min(azcen_3))
	thea_13 = np.linspace(minval,maxval,maxval*25+1,endpoint=True)#np.arange(maxval+1)
	f = sp_interp.interp1d(azcen_3, zcen_3, fill_value='extrapolate')
	zcen_13 = f(thea_13)
#    print(thea_13)
	for i in range(len(ind)):
		ind2 = np.where(np.abs(az_3[i]-thea_13)<=0.0222)[0]
		if (len(ind2>=1)):
			zzz_13.append(zz_3[i]-zcen_13[ind2[0]])
			mm_13.append(mass_3[i])
			wei_13 = mm_13/np.linalg.norm(mm_13)

	pl.scatter(az[ind],zz[ind],s=all_dist[ind],c=col,alpha=0.3,edgecolors=col1)
	pl.plot([170,-40],[0,0],'--')
	pl.errorbar(azcen_3,zcen_3,yerr=np.array(zrms_3)*2.355/2.,fmt='^',c=col_err,markersize=0.05,elinewidth = 1, capsize=1.1)
	pl.plot(thea_13,zcen_13,'-',color=col_cen)
	pl.text(25,-0.45,'%skpc' %(gr),fontsize=11,fontweight='bold',color='b')
	pl.xlim([-30,165])
	pl.ylim([-0.61,0.61])
	pl.yticks([-0.6,-0.4,-0.2,0,0.2,0.4,0.6],['','','','','','',''],fontsize=11,fontweight='bold')
	pl.xticks([0,50,100,150],['','','',''],fontsize=11,fontweight=1.8)
	pl.twiny()
	pl.xlim(0,radius)
	pl.xticks([0,20,40],['0','20','40 kpc'],fontsize=10.5,fontweight=1.8)    
	gr=14
	pl.subplot(2,5,6)
	radius = 2*np.pi*gr*(np.abs(az_min-az_max)/360.)
	thea,w = cal_warp(gr)
	thea,zw1,zw2,zw3,zw5=cal_warpc(gr)
	ind = np.where((rr>=(gr-0.5)) & (rr<(gr+0.5)))[0]
	zcen_4,zrms_4,vcen_4,vrms_4,rcen_4,rrms_4,azcen_4,azrms_4,zerr_4,verr_4,rerr_4,azerr_4,zcerr_4,vcerr_4,rcerr_4,azcerr_4=cal_zcen_zrms(az[ind],zz[ind],vv[ind],rr[ind],mass[ind],bii)
	az_4, zz_4, mass_4, zzz_14, mm_14, rr_14 = [], [], [], [], [], []
	az_4 = az[ind]
	mass_4 = mass[ind]
	zz_4 = zz[ind]
	rr_14 = rr[ind]
	maxval = int(np.max(azcen_4))
	minval = int(np.min(azcen_4))
	thea_14 = np.linspace(minval,maxval,maxval*25+1,endpoint=True)#np.arange(maxval+1)
	f = sp_interp.interp1d(azcen_4, zcen_4, fill_value='extrapolate')
	zcen_14 = f(thea_14)
#    print(thea_14)
	for i in range(len(ind)):
		ind2 = np.where(np.abs(az_4[i]-thea_14)<=0.0294)[0]
		if (len(ind2>=1)):
			zzz_14.append(zz_4[i]-zcen_14[ind2[0]])
			mm_14.append(mass_4[i])
			wei_14 = mm_14/np.linalg.norm(mm_14)

	pl.scatter(az[ind],zz[ind],s=all_dist[ind],c=col,alpha=0.3,edgecolors=col1)
	pl.plot([170,-40],[0,0],'--')
	pl.errorbar(azcen_4,zcen_4,yerr=np.array(zrms_4)*2.355/2.,fmt='^',c=col_err,markersize=0.05,elinewidth = 1, capsize=1.1)
	pl.plot(thea_14,zcen_14,'-',color=col_cen)
	pl.text(25,-0.45,'%skpc' %(gr),fontsize=11,fontweight='bold',color='b')
	pl.xlim([-30,165])
	pl.ylim([-0.61,0.61])
	pl.xticks([0,50,100,150],[r'0$^{\circ}$',r'50$^{\circ}$',r'100$^{\circ}$',r'150$^{\circ}$'],fontsize=11,fontweight=1.8)
	pl.yticks([-0.6,-0.4,-0.2,0,0.2,0.4,0.6],['-0.6','-0.4','-0.2','0','0.2','0.4','0.6'],fontsize=11,fontweight=1.8)
#    
#    pl.xticks([0,50,100,150],['0','50','100','150'],fontsize=11,fontweight=1.8)
	pl.xlabel('Galactocentric Azimuth',fontsize=13,fontweight='bold')
	pl.ylabel('Z (kpc)',fontsize=13,fontweight='bold')  
	pl.twiny()
	pl.xlim(0,radius)
	pl.xticks([0,20,40],['0','20','40'],fontsize=10.5,fontweight=1.8)
	gr=15
	pl.subplot(2,5,7)
	radius = 2*np.pi*gr*(np.abs(az_min-az_max)/360.)
	thea,w = cal_warp(gr)
	thea,zw1,zw2,zw3,zw5=cal_warpc(gr)
	ind = np.where((rr>=(gr-0.5)) & (rr<(gr+0.5)))[0]
	zcen_5,zrms_5,vcen_5,vrms_5,rcen_5,rrms_5,azcen_5,azrms_5,zerr_5,verr_5,rerr_5,azerr_5,zcerr_5,vcerr_5,rcerr_5,azcerr_5=cal_zcen_zrms(az[ind],zz[ind],vv[ind],rr[ind],mass[ind],bii)
	az_5, zz_5, mass_5, zzz_15, mm_15, rr_15 = [], [], [], [], [], []
	az_5 = az[ind]
	mass_5 = mass[ind]
	zz_5 = zz[ind]
	rr_15 = rr[ind]
	maxval = int(np.max(azcen_5))
	minval = int(np.min(azcen_5))
	thea_15 = np.linspace(minval,maxval,maxval*25+1,endpoint=True)#np.arange(maxval+1)
	f = sp_interp.interp1d(azcen_5, zcen_5, fill_value='extrapolate')
	zcen_15 = f(thea_15)
#    print(thea_15)
	for i in range(len(ind)):
		ind2 = np.where(np.abs(az_5[i]-thea_15)<=0.0294)[0]
		if (len(ind2>=1)):
			zzz_15.append(zz_5[i]-zcen_15[ind2[0]])
			mm_15.append(mass_5[i])
			wei_15 = mm_15/np.linalg.norm(mm_15)

	pl.scatter(az[ind],zz[ind],s=all_dist[ind],c=col,alpha=0.3,edgecolors=col1)
	pl.plot([170,-40],[0,0],'--')
	pl.errorbar(azcen_5,zcen_5,yerr=np.array(zrms_5)*2.355/2.,fmt='^',c=col_err,markersize=0.05,elinewidth = 1, capsize=1.1)
	pl.plot(thea_15,zcen_15,'-',color=col_cen)
	pl.text(25,-0.45,'%skpc' %(gr),fontsize=11,fontweight='bold',color='b')
	pl.xlim([-30,165])
	pl.ylim([-0.61,0.61])
	pl.yticks([-0.6,-0.4,-0.2,0,0.2,0.4,0.6],['','','','','','',''],fontsize=11,fontweight='bold')
	pl.xticks([0,50,100,150],['','','',''],fontsize=11,fontweight=1.8)
	pl.twiny()
	pl.xlim(0,radius)
	pl.xticks([0,20,40],['0','20','40'],fontsize=10.5,fontweight=1.8)    
	gr=16
	pl.subplot(2,5,8)
	radius = 2*np.pi*gr*(np.abs(az_min-az_max)/360.)
	thea,w = cal_warp(gr)
	thea,zw1,zw2,zw3,zw5=cal_warpc(gr)
	ind = np.where((rr>=(gr-0.5)) & (rr<(gr+0.5)))[0]
	zcen_6,zrms_6,vcen_6,vrms_6,rcen_6,rrms_6,azcen_6,azrms_6,zerr_6,verr_6,rerr_6,azerr_6,zcerr_6,vcerr_6,rcerr_6,azcerr_6=cal_zcen_zrms(az[ind],zz[ind],vv[ind],rr[ind],mass[ind],bii)
	az_6, zz_6, mass_6, zzz_16, mm_16, rr_16 = [], [], [], [], [], []
	az_6 = az[ind]
	mass_6 = mass[ind]
	zz_6 = zz[ind]
	rr_16= rr[ind]
	maxval = int(np.max(azcen_6))
	minval = int(np.min(azcen_6))
	thea_16 = np.linspace(minval,maxval,maxval*25+1,endpoint=True)#np.arange(maxval+1)
	f = sp_interp.interp1d(azcen_6, zcen_6, fill_value='extrapolate')
	zcen_16 = f(thea_16)
#    print(thea_16)
	for i in range(len(ind)):
		ind2 = np.where(np.abs(az_6[i]-thea_16)<=0.0319)[0]
		if (len(ind2>=1)):
			zzz_16.append(zz_6[i]-zcen_16[ind2[0]])
			mm_16.append(mass_6[i])
			wei_16 = mm_16/np.linalg.norm(mm_16)

	pl.scatter(az[ind],zz[ind],s=all_dist[ind],c=col,alpha=0.3,edgecolors=col1)
	pl.plot([170,-40],[0,0],'--')
	pl.errorbar(azcen_6,zcen_6,yerr=np.array(zrms_6)*2.355/2.,fmt='^',c=col_err,markersize=0.05,elinewidth = 1, capsize=1.1)
	pl.plot(thea_16,zcen_16,'-',color=col_cen)
	pl.text(25,-0.45,'%skpc' %(gr),fontsize=11,fontweight='bold',color='b')
	pl.xlim([-30,165])
	pl.ylim([-0.61,0.61])
	pl.yticks([-0.6,-0.4,-0.2,0,0.2,0.4,0.6],['','','','','','',''],fontsize=11,fontweight='bold')
	pl.xticks([0,50,100,150],['','','',''],fontsize=11,fontweight=1.8)
	pl.twiny()
	pl.xlim(0,radius)
	pl.xticks([0,20,40],['0','20','40'],fontsize=10.5,fontweight=1.8)
	##########
	gr=17.5
	pl.subplot(2,5,9)
	radius = 2*np.pi*gr*(np.abs(az_min-az_max)/360.)
	thea,w = cal_warp(gr)
	thea,zw1,zw2,zw3,zw5=cal_warpc(gr)
	ind = np.where((rr>=(gr-1.0)) & (rr<(gr+1.0)))[0]
	zcen_7,zrms_7,vcen_7,vrms_7,rcen_7,rrms_7,azcen_7,azrms_7,zerr_7,verr_7,rerr_7,azerr_7,zcerr_7,vcerr_7,rcerr_7,azcerr_7=cal_zcen_zrms(az[ind],zz[ind],vv[ind],rr[ind],mass[ind],bii)
	az_7, zz_7, mass_7, zzz_17, mm_17, rr_17 = [], [], [], [], [], []
	az_7 = az[ind]
	mass_7 = mass[ind]
	zz_7 = zz[ind]
	rr_17 = rr[ind]
	maxval = int(np.max(azcen_7))
	minval = int(np.min(azcen_7))
	thea_17 = np.linspace(minval,maxval,maxval*25+1,endpoint=True)#np.arange(maxval+1)
	f = sp_interp.interp1d(azcen_7, zcen_7, fill_value='extrapolate')
	zcen_17 = f(thea_17)
#    print(thea_17)
	for i in range(len(ind)):
		ind2 = np.where(np.abs(az_7[i]-thea_17)<=0.0319)[0]
		if (len(ind2>=1)):
			zzz_17.append(zz_7[i]-zcen_17[ind2[0]])
			mm_17.append(mass_7[i])
			wei_17 = mm_17/np.linalg.norm(mm_17)

	pl.scatter(az[ind],zz[ind],s=all_dist[ind],c=col,alpha=0.3,edgecolors=col1)
	pl.plot([170,-40],[0,0],'--')
	pl.errorbar(azcen_7,zcen_7,yerr=np.array(zrms_7)*2.355/2.,fmt='^',c=col_err,markersize=0.05,elinewidth = 1, capsize=1.1)
	pl.plot(thea_17,zcen_17,'-',color=col_cen)
	pl.text(25,-0.45,'%skpc' %(gr),fontsize=11,fontweight='bold',color='b')
	pl.xlim([-30,165])
	pl.ylim([-0.61,0.61])
	pl.yticks([-0.6,-0.4,-0.2,0,0.2,0.4,0.6],['','','','','','',''],fontsize=11,fontweight='bold')
	pl.xticks([0,50,100,150],['','','',''],fontsize=11,fontweight=1.8)
	pl.twiny()
	pl.xlim(0,radius)
	pl.xticks([0,20,40],['0','20','40'],fontsize=10.5,fontweight=1.8)
	gr=19.5
	pl.subplot(2,5,10)
	radius = 2*np.pi*gr*(np.abs(az_min-az_max)/360.)
	thea,w = cal_warp(gr)
	thea,zw1,zw2,zw3,zw5=cal_warpc(gr)
	ind = np.where((rr>=(gr-1.0)) & (rr<(gr+1.0)))[0]
	zcen_8,zrms_8,vcen_8,vrms_8,rcen_8,rrms_8,azcen_8,azrms_8,zerr_8,verr_8,rerr_8,azerr_8,zcerr_8,vcerr_8,rcerr_8,azcerr_8=cal_zcen_zrms(az[ind],zz[ind],vv[ind],rr[ind],mass[ind],bii)
	az_8, zz_8, mass_8, zzz_18, mm_18, rr_18 = [], [], [], [], [], []
	az_8 = az[ind]
	mass_8 = mass[ind]
	zz_8 = zz[ind]
	rr_18 = rr[ind]
	maxval = int(np.max(azcen_8))
	minval = int(np.min(azcen_8))
	thea_18 = np.linspace(minval,maxval,maxval*25+1,endpoint=True)#np.arange(maxval+1)
	f = sp_interp.interp1d(azcen_8, zcen_8, fill_value='extrapolate')
	zcen_18 = f(thea_18)
#    print(thea_18)
	for i in range(len(ind)):
		ind2 = np.where(np.abs(az_8[i]-thea_18)<=0.0319)[0]
		if (len(ind2>=1)):
			zzz_18.append(zz_8[i]-zcen_18[ind2[0]])
			mm_18.append(mass_8[i])
			wei_18 = mm_18/np.linalg.norm(mm_18)

	pl.scatter(az[ind],zz[ind],s=all_dist[ind],c=col,alpha=0.3,edgecolors=col1)
	pl.plot([170,-40],[0,0],'--')
	pl.errorbar(azcen_8,zcen_8,yerr=np.array(zrms_8)*2.355/2.,fmt='^',c=col_err,markersize=0.05,elinewidth = 1, capsize=1.1)
	pl.plot(thea_18,zcen_18,'-',color=col_cen)
	pl.text(25,-0.45,'%skpc' %(gr),fontsize=11,fontweight='bold',color='b')
	pl.xlim([-30,165])
	pl.ylim([-0.61,0.61])
	pl.yticks([-0.6,-0.4,-0.2,0,0.2,0.4,0.6],['','','','','','',''],fontsize=11,fontweight='bold')
	pl.xticks([0,50,100,150],['','','',''],fontsize=11,fontweight=1.8)
	pl.twiny()
	pl.xlim(0,radius)
	pl.xticks([0,20,40,60],['0','20','40','60 kpc'],fontsize=10.5,fontweight=1.8)     
	
	gr=19
	pl.subplot(3,4,11)
	thea,w = cal_warp(gr)
	thea,zw1,zw2,zw3,zw5=cal_warpc(gr)
	ind = np.where((rr>=(gr-0.5)) & (rr<(gr+0.5)))[0]
	zcen_19,zrms_19,vcen_19,vrms_19,rcen_19,rrms_19,azcen_19,azrms_19,zerr_19,verr_19,rerr_19,azerr_19,zcerr_19,vcerr_19,rcerr_19,azcerr_19=cal_zcen_zrms(az[ind],zz[ind],vv[ind],rr[ind],mass[ind],bii)
	az_19, zz_19, mass_19, zzz_19, mm_19, rr_19 = [], [], [], [], [], []
	az_19 = az[ind]
	mass_19 = mass[ind]
	zz_19 = zz[ind]
	rr_19 = rr[ind]
	maxval = int(np.max(azcen_19))
	minval = int(np.min(azcen_19))
	thea_19 = np.linspace(minval,maxval,maxval*25+1,endpoint=True)#np.arange(maxval+1)
	f = sp_interp.interp1d(azcen_19, zcen_19, fill_value='extrapolate')
	zcen_199 = f(thea_19)
#    print(thea_19)
	for i in range(len(ind)):
		ind2 = np.where(np.abs(az_19[i]-thea_19)<=0.0319)[0]
		if (len(ind2>=1)):
			zzz_19.append(zz_19[i]-zcen_199[ind2[0]])
			mm_19.append(mass_19[i])
			wei_19 = mm_19/np.linalg.norm(mm_19)
  
	pl.scatter(az[ind],zz[ind],s=all_dist[ind],c=col,alpha=0.3,edgecolors=col1)
	pl.text(25,-0.45,'%skpc' %(gr),fontsize=11,fontweight='bold',color='b')
	pl.xlim([165,-30])
	pl.ylim([-0.61,0.61])
	pl.yticks([-0.6,-0.4,-0.2,0,0.2,0.4,0.6],['','','','',''],fontsize=11,fontweight='bold')
	pl.xticks([150,100,50,0],['150','100','50','0'],fontsize=11,fontweight='bold')
	pl.plot([165,-30],[0,0],'--')
	#pl.errorbar(azcen_19,zcen_19,yerr=np.array(zrms_19)*2.355/2.,fmt='^',c=col_err,markersize=5.,elinewidth = 1, capsize=1.1)
	pl.plot(thea_19,zcen_199,'-',color=col_cen)
	gr=20
	pl.subplot(3,4,12)
	thea,w = cal_warp(gr)
	thea,zw1,zw2,zw3,zw5=cal_warpc(gr)
	ind = np.where((rr>=(gr-0.5)) & (rr<(gr+2.5)))[0]
	zcen_20,zrms_20,vcen_20,vrms_20,rcen_20,rrms_20,azcen_20,azrms_20,zerr_20,verr_20,rerr_20,azerr_20,zcerr_20,vcerr_20,rcerr_20,azcerr_20=cal_zcen_zrms(az[ind],zz[ind],vv[ind],rr[ind],mass[ind],bii)
	az_20, zz_20, mass_20, zzz_20, mm_20, rr_20 = [], [], [], [], [], []
	az_20 = az[ind]
	mass_20 = mass[ind]
	zz_20 = zz[ind]
	rr_20 = rr[ind]
	maxval = int(np.max(azcen_20))
	minval = int(np.min(azcen_20))
	thea_20 = np.linspace(minval,maxval,maxval*25+1,endpoint=True)#np.arange(maxval+1)
	f = sp_interp.interp1d(azcen_20, zcen_20, fill_value='extrapolate')
	zcen_200 = f(thea_20)
#    print(thea_20)
	for i in range(len(ind)):
		ind2 = np.where(np.abs(az_20[i]-thea_20)<=0.0319)[0]
		if (len(ind2>=1)):
			zzz_20.append(zz_20[i]-zcen_200[ind2[0]])
			mm_20.append(mass_20[i])
			wei_20 = mm_20/np.linalg.norm(mm_20)        
	pl.scatter(az[ind],zz[ind],s=all_dist[ind],c=col,alpha=0.3,edgecolors=col1)
	pl.plot([165,-30],[0,0],'--')
	#pl.errorbar(azcen_20,zcen_20,yerr=np.array(zrms_20)*2.355/2.,fmt='^',c=col_err,markersize=5.,elinewidth = 1, capsize=1.1)
	pl.text(25,-0.45,'%skpc' %(gr),fontsize=11,fontweight='bold',color='b')
	pl.xlim([165,-30])
	pl.ylim([-0.61,0.61])
	pl.yticks([-0.6,-0.4,-0.2,0,0.2,0.4,0.6],['','','','',''],fontsize=11,fontweight='bold')
	pl.xticks([150,100,50,0],['150','100','50','0'],fontsize=11,fontweight='bold')
	pl.plot(thea_20,zcen_200,'-',color=col_cen)
	
	#fig1.savefig('az_z_warp_residual.png',format='png',bbox_inches='tight')
	pl.show()


'''
