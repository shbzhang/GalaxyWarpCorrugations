### obsolete, replaced by warp1.py
### plot Az-Z of MCs with Xfactor correction at different Rgal

import numpy as np
import pylab as pl
import scipy.interpolate as sp_interp
from scipy.stats import norm
import scipy.optimize as sp_opt
from shared import *
import matplotlib as mpl
import matplotlib.pyplot as plt

if __name__ == '__main__':
	az_min, az_max = -30, 165

	out_data = np.loadtxt('out_para_xf.txt',comments='#')
	osc_data = np.loadtxt('osc2_para_xf.txt',comments='#')
	per_data = np.loadtxt('per_para_xf.txt',comments='#')
	out_bb = out_data[:,1]
	out_ll = out_data[:,0]
	osc_bb = osc_data[:,1]
	osc_ll = osc_data[:,0]
	per_bb = per_data[:,1]
	per_ll = per_data[:,0]


	data = np.vstack((osc_data, out_data, per_data)).T
	az = data[19]
	rr = data[11]
	zz = data[14]
	vv = data[2]
	mass = data[7]
	com = data[8]
	#all_mass = np.hstack((osc_data[:,7]*1,out_data[:,7]*1.,per_data[:,7]))    #####################################pl mass use different scale
	scatter_size, scatter_color = mass_distribute(mass)

	theta_axis = np.linspace(-30,170,200)

	fig1, ax = plt.subplots(nrows=2, ncols=5, sharex=True, sharey=True, figsize=[17,7])
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
		ax[i].scatter(az[idx], zz[idx], s=scatter_size[idx]*2, **ring_kws_mc)

		### plot binned average
		azcen, azrms, zcen, zrms, vcen, vrms, rcen, rrms = cal_zcen_zrms(az[idx], zz[idx], vv[idx], rr[idx], weights=mass[idx], binsize = 6, nbin=60)

		# i think there is no need to fill gap with interpolating, just plot the binned average
		'''
		max_azcen = int(np.nanmax(azcen))
		min_azcen = int(np.nanmin(azcen))
		f = sp_interp.interp1d(azcen, zcen, fill_value='extrapolate')
		bin_axis = np.linspace(min_azcen, max_azcen, max_azcen*25+1, endpoint=True)
		bin_val = f(bin_axis)
		'''
		ax[i].plot(azcen, zcen, '-', **ring_kws_bin)
		### fill gap
		idx = np.isfinite(zcen)
		ax[i].plot(azcen[idx], zcen[idx], '--', **ring_kws_bin)

		### plot warp models
		# HI
		if i>0:
			w = cal_warp(gr, theta_axis)
			ax[i].plot(theta_axis, w, **ring_kws_hi)
		# cepheids, co
		zw1, zw2, zw4, zw5=cal_warpc(gr, theta_axis)
		ax[i].plot(theta_axis, zw1, **ring_kws_ceph)
		ax[i].plot(theta_axis, zw4, **ring_kws_co1)
		ax[i].plot(theta_axis, zw5, **ring_kws_co2)

		### plot gr text
		text = '%i' % gr if gr%1==0 else '%.1f' % gr
		if i==0: text += ' kpc'
		ax[i].text(0.02, 0.95, text, transform=ax[i].transAxes, **ring_kws_text)

		### axes
		ax[i].set_xticks(np.arange(-50, 200, 50))
		ax[i].set_yticks(np.arange(-3, 3, 0.5))
		ax[i].minorticks_on()
		#ax[i].tick_params(axis='both', which='major', labelsize=10.5, width=1.8)
		ax[i].set_xlim(az_min, az_max)
		ax[i].set_ylim([-0.8, 1.4])
		ax[i].tick_params(right=True, direction='in')#, labelsize=1000/self.dpi)
		ax[i].tick_params(which='minor', right=True, direction='in')#, labelsize=1000/self.dpi)

		#ax[i].set_yticklabels(fontsize=10.5, fontweight=1.8)
		upper = ax[i].twiny()
		upper.set_xticks(np.arange(0, 100, 20))
		upper.set_xticks(np.arange(0, 100, 5), minor=True)
		upper.set_xlim(0, radius)
		upper.tick_params(direction='in')
		upper.tick_params(which='minor', direction='in')

		if i >= 5:
			ax[i].set_xticklabels(['%i' % i for i in np.arange(-50, 200, 50)])
		else:
			ax[i].set_xticklabels([])
		if i == 5:
			ax[i].set_xlabel('Galactocentric Azimuth (deg)', fontsize=13, fontweight='bold')
			ax[i].set_ylabel('Z (kpc)', fontsize=13, fontweight='bold')  

		### slightly shift upper ticklabel position to avoid overlap
		if i<=2: upper.set_xticklabels(['  0', '20', None, None, None]) # avoid overlap
		elif i==3: upper.set_xticklabels(['  0', '20', '40    ', None, None]) # avoid overlap
		elif i==4: upper.set_xticklabels(['  0', '20', '40 kpc', None, None]) # add kpc at the end
		elif i<=8: upper.set_xticklabels(['  0', '20', '40', None, None]) # add kpc at the end
		elif i==9: upper.set_xticklabels(['  0', '20', '40', '60 kpc', None]) # add kpc at the end

		if i==9: ax[i].legend()

	fig1.savefig('fig/az_z_warp_xf.png',format='png',bbox_inches='tight', dpi=400)
	plt.show()



'''
import numpy as np
import pylab as pl
import scipy.interpolate as sp_interp
from scipy.stats import norm
import scipy.optimize as sp_opt
def function(x,a,b):

		y =  np.log(a)-x/b

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
#	p3=[-2.97397887e-02,1.72161063e-02,-1.66075673e+01,7.93265614e-03,6.05519348e+01]
#	p4=[-0.07616907,  0.10783739,  7.79992154,  0.9348509,  -9.02005785] ##1st notcorrect
	p4=[-0.06976495,  0.13775231,  8.12779254,  0.78453885, -8.47393439]##1st xf
	p5=[-0.069110292892,0.10257483877,7.868629069,1.0,-8.5999172236, 0.21143504143391126, 16.42369138978075, 1.0, -98.64681680983567]##a0, a1, Rw1, bw1, PHIw1, a2, Rw2, bw2, PHIw2
	thea, zw1, zw2, zw4,zw5 =[], [], [], [], []
	thea = np.linspace(-30,170,200)
	zw1=p1[2]*(rg-p1[0])*np.sin(np.deg2rad(thea-p1[1]))
	zw2=p2[2]*((rg-p2[0])**p2[3])*np.sin(np.deg2rad(thea-p2[1]))
#	zw2=-p3[0]F+(rg - 8)**2*(p3[1]*np.sin(np.deg2rad(thea-p3[2]))+p3[3]*np.sin(2.0*(thea-p3[4])/180.*np.pi))  #mwisp
	zw4=p4[0]+p4[1] * (rg-p4[2])**p4[3] * np.sin((thea-p4[4])/180*np.pi) #mwisp first component
	zw5=np.zeros(thea.size)        #mwisp 2nd comp
	if (rg>p5[6]):
		zw5 = p5[0]+p5[1] * (rg-p5[2])**p5[3] * np.sin((thea-p5[4])/180*np.pi)+p5[5] * (rg-p5[6])**p5[7] * np.sin(2*(thea-p5[8])/180*np.pi)#
	elif (rg<p5[6]):
		zw5 = p5[0]+p5[1] * (rg-p5[2])**p5[3] * np.sin((thea-p5[4])/180*np.pi)  
	return(thea,zw1,zw2,zw4,zw5)
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
#            thea.append((thea2+bi/2.))
			z_err.append(temp_zrms/np.sqrt(len(np.transpose(ind))))
			v_err.append(temp_vrms/np.sqrt(len(np.transpose(ind))))
			r_err.append(temp_rrms/np.sqrt(len(np.transpose(ind))))
			az_err.append(temp_azrms/np.sqrt(len(np.transpose(ind))))
			zcen_err.append(np.abs(temp_zcen)/np.sqrt(len(np.transpose(ind))))
			vcen_err.append(temp_vcen/np.sqrt(len(np.transpose(ind))))
			rcen_err.append(temp_rcen/np.sqrt(len(np.transpose(ind))))
			azcen_err.append(temp_azcen/np.sqrt(len(np.transpose(ind))))
#            xcen.append(temp_rcen*np.sin(np.deg2rad(thea2+bi/2.)))
#            ycen.append(temp_rcen*np.cos(np.deg2rad(thea2+bi/2.)))
	return(zcen,zrms,vcen,vrms,rcen,rrms,azcen,azrms,z_err,v_err,r_err,az_err,zcen_err,vcen_err,rcen_err,azcen_err)
if __name__ == '__main__':
   
	fcolor = '#1f8097'#'#1f77b4'
	edgecolor = 'none'
	col_hi = '#df7f0e'
	sty_hi = ':'
	col_cepheid = '#651fff'
	sty_cepheid = ':'
	col_co = '#e33a23'
	sty_co = '--'
	col_cen = '#0a3452'
	col_text = '#1e3f66'
	#fcolor,edgecolors,col_hi,col_cepheid,col_co,col_cen = ['coral','none','teal','blue','crimson','black']
	out_data = np.loadtxt('out_para_xf.txt',comments='#')
	osc_data = np.loadtxt('osc2_para_xf.txt',comments='#')
	per_data = np.loadtxt('per_para_xf.txt',comments='#')
	out_bb = out_data[:,1]
	out_ll = out_data[:,0]
	osc_bb = osc_data[:,1]
	osc_ll = osc_data[:,0]
	per_bb = per_data[:,1]
	per_ll = per_data[:,0]
	az_min, az_max = -30, 165

	#az = np.loadtxt('az_out1.txt')   ##9 to 19kpc   
	az = np.hstack((osc_data[:,19],out_data[:,19],per_data[:,19]))
	rr = np.hstack((osc_data[:,11],out_data[:,11],per_data[:,11]))
	zz = np.hstack((osc_data[:,14],out_data[:,14],per_data[:,14]))
	vv = np.hstack((osc_data[:,2],out_data[:,2],per_data[:,2]))
	mass = np.hstack((osc_data[:,7],out_data[:,7],per_data[:,7]))
	com = np.hstack((osc_data[:,8],out_data[:,8],per_data[:,8]))
#    for j in range(len(zz)):
#       	if az[j]>0:
#    	     print('%.2f\t%.2f\t%.3f\t%.2f'%(rr[j],az[j],zz[j],mass[j]))
#        else:
#    	     print('%.2f\t%.2f\t%.3f\t%.2f'%(rr[j],az[j]+360.,zz[j],mass[j]))     ###python warp.py > cloud.razm
	all_mass = np.hstack((osc_data[:,7]*1,out_data[:,7]*1.,per_data[:,7]))    #####################################pl mass use different scale
	all_dist,all_colors = mass_distribute(all_mass)
	bii=6
	fig1, ax = pl.subplots(sharex=True, figsize=[17,7])
	pl.subplots_adjust(wspace=0,hspace=0.2)
	pl.subplot(2,5,1)
	gr=9
	radius = 2*np.pi*gr*(np.abs(az_min-az_max)/360.)
	thea,w = cal_warp(gr)
	thea,zw1,zw2,zw4,zw5=cal_warpc(gr)
	ind = np.where((rr>=(gr-0.5)) & (rr<(gr+0.5)))[0]
	zcen_9,zrms_9,vcen_9,vrms_9,rcen_9,rrms_9,azcen_9,azrms_9,zerr_9,verr_9,rerr_9,azerr_9,zcerr_9,vcerr_9,rcerr_9,azcerr_9=cal_zcen_zrms(az[ind],zz[ind],vv[ind],rr[ind],mass[ind],bii)
	az_9, zz_9, mass_9, zzz_9, mm_9, rr_9, com_9 = [], [], [], [], [], [], []
	az_9 = az[ind]
	mass_9 = mass[ind]
	zz_9 = zz[ind]
	rr_9 = rr[ind]
	com_9 = com[ind]
	maxval = int(np.max(azcen_9))
	minval = int(np.min(azcen_9))
	thea_9 = np.linspace(minval,maxval,maxval*25+1,endpoint=True)#np.arange(maxval+1)
	f = sp_interp.interp1d(azcen_9, zcen_9, fill_value='extrapolate')
	zcen_99 = f(thea_9)
	for i in range(len(ind)):
		ind2 = np.where(np.abs(az_9[i]-azcen_9)<=0.02383)[0]
		if (len(ind2>=1)):
			zzz_9.append(zz_9[i]-zcen_9[ind2[0]])
			mm_9.append(mass_9[i])
			wei_9 = mm_9/np.linalg.norm(mm_9)
#    for j in range(len(zz_9)):
#    	print('%.4f\t%.2f'%(zzz_9[j],mm_9[j]))        ##python warp.py > 9.zm
#    	if az_9[j]>0:
#    	     print('%.2f\t%.2f\t%.3f\t%.2f\t%.3f'%(rr_9[j],az_9[j],zz_9[j],mass_9[j],com_9[j]))
#    	else:
#    	     print('%.2f\t%.2f\t%.3f\t%.2f\t%.3f'%(rr_9[j],az_9[j]+360.,zz_9[j],mass_9[j],com_9[j]))     ###python warp.py 9.razmc
	pl.scatter(az[ind],zz[ind],s=all_dist[ind],c=fcolor,alpha=0.3,edgecolors=edgecolor)
#    pl.plot(thea,w,ls=sty_hi,lw=1.5,color=col_hi)
	pl.plot(thea,zw1,ls=sty_cepheid,lw=1.5,color=col_cepheid)
#    pl.plot(thea,zw5,'-',lw=1.5,color='red')
	pl.plot(thea,zw4,ls=sty_co,lw=1.5,color=col_co)
#    pl.plot([170,-40],[0,0],ls='dotted')
	pl.plot(thea_9,zcen_99,'-',color=col_cen, zorder=1)
	pl.text(-25,1.3,'%s kpc' %(gr),fontsize=11,fontweight='bold',color=col_text)
	pl.xlim([-30,165])
	pl.ylim([-0.9,1.55])
	pl.yticks([-0.5,0,0.5,1.0,1.5],['-0.5','0','0.5','1.0','1.5'],fontsize=11,fontweight=1.8)
	pl.xticks([0,50,100,150],['','','',''],fontsize=11,fontweight=1.8)
	pl.twiny()
	pl.xlim(0,radius)
	pl.xticks([0,20],['0','20'],fontsize=10.5,fontweight=1.8)
	################################
	pl.subplot(2,5,2)
	gr=10
	radius = 2*np.pi*gr*(np.abs(az_min-az_max)/360.)
	thea,w = cal_warp(gr)
#    thea,zw1,zw2,zw4,zw5=cal_warpc(gr)
	ind = np.where((rr>=(gr-0.5)) & (rr<(gr+0.5)))[0]
	zcen_0,zrms_0,vcen_0,vrms_0,rcen_0,rrms_0,azcen_0,azrms_0,zerr_0,verr_0,rerr_0,azerr_0,zcerr_0,vcerr_0,rcerr_0,azcerr_0=cal_zcen_zrms(az[ind],zz[ind],vv[ind],rr[ind],mass[ind],bii)
	az_0, zz_0, mass_0, zzz_10, mm_10, rr_10, com_10 = [], [], [], [], [], [], []
	az_0 = az[ind]
	mass_0 = mass[ind]
	zz_0 = zz[ind]
	rr_10 =rr[ind]
	com_10 = com[ind]
	maxval = int(np.max(azcen_0))
	minval = int(np.min(azcen_0))
	thea_10 = np.linspace(minval, maxval, maxval*25+1, endpoint=True)
	f = sp_interp.interp1d(azcen_0, zcen_0, fill_value='extrapolate')
	zcen_10 = f(thea_10)
	for i in range(len(ind)):
		ind2 = np.where(np.abs(az_0[i]-thea_10)<=0.02383)[0]
		if (len(ind2>=1)):
			zzz_10.append(zz_0[i]-zcen_10[ind2[0]])
			mm_10.append(mass_0[i])
			wei_10 = mm_10/np.linalg.norm(mm_10)
#    for j in range(len(zz_0)):
#    	print('%.4f\t%.2f'%(zzz_10[j],mm_10[j]))        ##python warp.py > 10.zm
#    	if az_0[j]>0:
#    	     print('%.2f\t%.2f\t%.3f\t%.2f\t%.3f'%(rr_10[j],az_0[j],zz_0[j],mass_0[j],com_10[j]))
#    	else:
#    	     print('%.2f\t%.2f\t%.3f\t%.2f\t%.3f'%(rr_10[j],az_0[j]+360.,zz_0[j],mass_0[j],com_10[j]))     ###python warp.py > 10.razmc
	pl.scatter(az[ind],zz[ind],s=all_dist[ind],c=fcolor,alpha=0.3,edgecolors=edgecolor)
	pl.plot(thea,w,ls=sty_hi,lw=1.5,color=col_hi)
	pl.plot(thea,zw1,ls=sty_cepheid,lw=1.5,color=col_cepheid)
#    pl.plot(thea,zw5,'-',lw=1.5,color='red')
	pl.plot(thea,zw4,ls=sty_co,lw=1.5,color=col_co)
#    pl.plot([170,-40],[0,0],ls='dotted')
	pl.plot(thea_10,zcen_10,'-',color=col_cen, zorder=1)
	pl.text(-25,1.3,'%s' %(gr),fontsize=11,fontweight='bold',color=col_text)
	pl.xlim([-30,165])
	pl.ylim([-0.9,1.55])
	pl.yticks([-0.5,0,0.5,1.0,1.5],['','','','',''],fontsize=11,fontweight=1.8)
	pl.xticks([0,50,100,150],['','','',''],fontsize=11,fontweight=1.8)
	pl.twiny()
	pl.xlim(0,radius)
	pl.xticks([0,20],['0','20'],fontsize=10.5,fontweight=1.8)    
	gr=11
	pl.subplot(2,5,3)
	radius = 2*np.pi*gr*(np.abs(az_min-az_max)/360.)
	thea,w = cal_warp(gr)
	thea,zw1,zw2,zw4,zw5=cal_warpc(gr)
	ind = np.where((rr>=(gr-0.5)) & (rr<(gr+0.5)))[0]
	zcen_1,zrms_1,vcen_1,vrms_1,rcen_1,rrms_1,azcen_1,azrms_1,zerr_1,verr_1,rerr_1,azerr_1,zcerr_1,vcerr_1,rcerr_1,azcerr_1=cal_zcen_zrms(az[ind],zz[ind],vv[ind],rr[ind],mass[ind],bii)
	az_1, zz_1, mass_1, zzz_11, mm_11, rr_11, com_11 = [], [], [], [], [], [], []
	az_1 = az[ind]
	mass_1 = mass[ind]
	zz_1 = zz[ind]
	rr_11 =rr[ind]
	com_11 = com[ind]
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
#    for j in range(len(zz_1)):
#    	print('%.4f\t%.2f'%(zzz_11[j],mm_11[j]))        ##python warp.py > 11.zm
#    	if az_1[j]>0:
#    	     print('%.2f\t%.2f\t%.3f\t%.2f\t%.3f'%(rr_11[j],az_1[j],zz_1[j],mass_1[j],com_11[j]))
#    	else:
#    	     print('%.2f\t%.2f\t%.3f\t%.2f\t%.3f'%(rr_11[j],az_1[j]+360.,zz_1[j],mass_1[j],com_11[j]))     ###python warp.py 11.razmc
	pl.scatter(az[ind],zz[ind],s=all_dist[ind],c=fcolor,alpha=0.3,edgecolors=edgecolor)
	pl.plot(thea,w,ls=sty_hi,lw=1.5,color=col_hi)
	pl.plot(thea,zw1,ls=sty_cepheid,lw=1.5,color=col_cepheid)
#    pl.plot(thea,zw5,'-',lw=1.5,color='red')
	pl.plot(thea,zw4,ls=sty_co,lw=1.5,color=col_co)
	pl.plot(thea_11,zcen_11,'-',color=col_cen, zorder=1)
	pl.text(-25,1.3,'%s' %(gr),fontsize=11,fontweight='bold',color=col_text)
	pl.xlim([-30,165])
	pl.ylim([-0.9,1.55])
	pl.yticks([-0.5,0,0.5,1.0,1.5],['','','','',''],fontsize=11,fontweight=1.8)
	pl.xticks([0,50,100,150],['','','',''],fontsize=11,fontweight=1.8)
	pl.twiny()
	pl.xlim(0,radius)
	pl.xticks([0,20],['0','20'],fontsize=10.5,fontweight=1.8)    
	gr=12
	pl.subplot(2,5,4)
	radius = 2*np.pi*gr*(np.abs(az_min-az_max)/360.)
	thea,w = cal_warp(gr)
	thea,zw1,zw2,zw4,zw5=cal_warpc(gr)
	ind = np.where((rr>=(gr-0.5)) & (rr<(gr+0.5)))[0]
	zcen_2,zrms_2,vcen_2,vrms_2,rcen_2,rrms_2,azcen_2,azrms_2,zerr_2,verr_2,rerr_2,azerr_2,zcerr_2,vcerr_2,rcerr_2,az_cerr_2=cal_zcen_zrms(az[ind],zz[ind],vv[ind],rr[ind],mass[ind],bii)
	az_2, zz_2, mass_2, zzz_12, mm_12, rr_12, com_12 = [], [], [], [], [], [], []
	az_2 = az[ind]
	mass_2 = mass[ind]
	zz_2 = zz[ind]
	rr_12 = rr[ind]
	com_12 = com[ind]
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
#    for j in range(len(zz_2)):
#    	print('%.4f\t%.2f'%(zzz_12[j],mm_12[j]))        ##python warp.py > 12.zm
#    	if az_2[j]>0:
#    	     print('%.2f\t%.2f\t%.3f\t%.2f\t%.3f'%(rr_12[j],az_2[j],zz_2[j],mass_2[j],com_12[j]))
#    	else:
#    	     print('%.2f\t%.2f\t%.3f\t%.2f\t%.3f'%(rr_12[j],az_2[j]+360.,zz_2[j],mass_2[j],com_12[j]))     ###python warp.py 12.razmc
	pl.scatter(az[ind],zz[ind],s=all_dist[ind],c=fcolor,alpha=0.3,edgecolors=edgecolor)
	pl.plot(thea,w,ls=sty_hi,lw=1.5,color=col_hi)
	pl.plot(thea,zw1,ls=sty_cepheid,lw=1.5,color=col_cepheid)
#    pl.plot(thea,zw5,'-',lw=1.5,color='red')
	pl.plot(thea,zw4,ls=sty_co,lw=1.5,color=col_co)
	pl.plot(thea_12,zcen_12,'-',color=col_cen, zorder=1)
	pl.text(-25,1.3,'%s' %(gr),fontsize=11,fontweight='bold',color=col_text)
	pl.xlim([-30,165])
	pl.ylim([-0.9,1.55])
	pl.yticks([-0.5,0,0.5,1.0,1.5],['','','','',''],fontsize=11,fontweight=1.8)
	pl.xticks([0,50,100,150],['','','',''],fontsize=11,fontweight=1.8)
	pl.twiny()
	pl.xlim(0,radius)
	pl.xticks([0,20,40],['0','20','40'],fontsize=10.5,fontweight=1.8)    
	gr=13
	pl.subplot(2,5,5)
	radius = 2*np.pi*gr*(np.abs(az_min-az_max)/360.)
	thea,w = cal_warp(gr)
	thea,zw1,zw2,zw4,zw5=cal_warpc(gr)
	ind = np.where((rr>=(gr-0.5)) & (rr<(gr+0.5)))[0]
	zcen_3,zrms_3,vcen_3,vrms_3,rcen_3,rrms_3,azcen_3,azrms_3,zerr_3,verr_3,rerr_3,azerr_3,zcerr_3,vcerr_3,rcerr_3,azcerr_3=cal_zcen_zrms(az[ind],zz[ind],vv[ind],rr[ind],mass[ind],bii)
	az_3, zz_3, mass_3, zzz_13, mm_13, rr_13, com_13 = [], [], [], [], [], [], []
	az_3 = az[ind]
	mass_3 = mass[ind]
	zz_3 = zz[ind]
	rr_13 = rr[ind]
	com_13 = com[ind]
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
#    for j in range(len(zz_3)):
#    	print('%.4f\t%.2f'%(zzz_13[j],mm_13[j]))        ##python warp.py > 13.zm
#    	if az_3[j]>0:
#    	     print('%.2f\t%.2f\t%.3f\t%.2f\t%.3f'%(rr_13[j],az_3[j],zz_3[j],mass_3[j],com_13[j]))
#    	else:
#    	     print('%.2f\t%.2f\t%.3f\t%.2f\t%.3f'%(rr_13[j],az_3[j]+360.,zz_3[j],mass_3[j],com_13[j]))     ###python warp.py 13.razmc
	pl.scatter(az[ind],zz[ind],s=all_dist[ind],c=fcolor,alpha=0.3,edgecolors=edgecolor)
	pl.plot(thea,w,ls=sty_hi,lw=1.5,color=col_hi)
	pl.plot(thea,zw1,ls=sty_cepheid,lw=1.5,color=col_cepheid)
#    pl.plot(thea,zw5,'-',lw=1.5,color='red')
	pl.plot(thea,zw4,ls=sty_co,lw=1.5,color=col_co)
	pl.plot(thea_13,zcen_13,'-',color=col_cen, zorder=1)
	pl.text(-25,1.3,'%s' %(gr),fontsize=11,fontweight='bold',color=col_text)
	pl.xlim([-30,165])
	pl.ylim([-0.9,1.55])
	pl.yticks([-0.5,0,0.5,1.0,1.5],['','','','',''],fontsize=11,fontweight=1.8)
	pl.xticks([0,50,100,150],['','','',''],fontsize=11,fontweight=1.8)
	pl.twiny()
	pl.xlim(0,radius)
	pl.xticks([0,20,40],['0','20','40 kpc'],fontsize=10.5,fontweight=1.8)    
	gr=14
	pl.subplot(2,5,6)
	radius = 2*np.pi*gr*(np.abs(az_min-az_max)/360.)
	thea,w = cal_warp(gr)
	thea,zw1,zw2,zw4,zw5=cal_warpc(gr)
	ind = np.where((rr>=(gr-0.5)) & (rr<(gr+0.5)))[0]
	zcen_4,zrms_4,vcen_4,vrms_4,rcen_4,rrms_4,azcen_4,azrms_4,zerr_4,verr_4,rerr_4,azerr_4,zcerr_4,vcerr_4,rcerr_4,azcerr_4=cal_zcen_zrms(az[ind],zz[ind],vv[ind],rr[ind],mass[ind],bii)
	az_4, zz_4, mass_4, zzz_14, mm_14, rr_14, com_14 = [], [], [], [], [], [], []
	az_4 = az[ind]
	mass_4 = mass[ind]
	zz_4 = zz[ind]
	rr_14 = rr[ind]
	com_14 = com[ind]
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
#    for j in range(len(zz_4)):
#    	print('%.4f\t%.2f'%(zzz_14[j],mm_14[j]))        ##python warp.py > 14.zm
#    	if az_4[j]>0:
#    	     print('%.2f\t%.2f\t%.3f\t%.2f\t%.3f'%(rr_14[j],az_4[j],zz_4[j],mass_4[j],com_14[j]))
#    	else:
#    	     print('%.2f\t%.2f\t%.3f\t%.2f\t%.3f'%(rr_14[j],az_4[j]+360.,zz_4[j],mass_4[j],com_14[j]))     ###python warp.py 14.razmc
	pl.scatter(az[ind],zz[ind],s=all_dist[ind],c=fcolor,alpha=0.3,edgecolors=edgecolor)
	pl.plot(thea,w,ls=sty_hi,lw=1.5,color=col_hi)
	pl.plot(thea,zw1,ls=sty_cepheid,lw=1.5,color=col_cepheid)
#    pl.plot(thea,zw5,'-',lw=1.5,color='red')
	pl.plot(thea,zw4,ls=sty_co,lw=1.5,color=col_co)
	pl.plot(thea_14,zcen_14,'-',color=col_cen, zorder=1)
	pl.text(-25,1.3,'%s' %(gr),fontsize=11,fontweight='bold',color=col_text)
	pl.xlim([-30,165])
	pl.ylim([-0.9,1.55])
	pl.yticks([-0.5,0,0.5,1.0,1.5],['-0.5','0','0.5','1.0','1.5'],fontsize=11,fontweight=1.8)
	pl.xticks([0,50,100,150],[r'0$^{\circ}$',r'50$^{\circ}$',r'100$^{\circ}$',r'150$^{\circ}$'],fontsize=11,fontweight=1.8)
	pl.xlabel('Galactocentric Azimuth',fontsize=13,fontweight='bold')
	pl.ylabel('Z (kpc)',fontsize=13,fontweight='bold')  
	pl.twiny()
	pl.xlim(0,radius)
	pl.xticks([0,20,40],['0','20','40'],fontsize=10.5,fontweight=1.8)
	gr=15
	pl.subplot(2,5,7)
	radius = 2*np.pi*gr*(np.abs(az_min-az_max)/360.)
	thea,w = cal_warp(gr)
	thea,zw1,zw2,zw4,zw5=cal_warpc(gr)
	ind = np.where((rr>=(gr-0.5)) & (rr<(gr+0.5)))[0]
	zcen_5,zrms_5,vcen_5,vrms_5,rcen_5,rrms_5,azcen_5,azrms_5,zerr_5,verr_5,rerr_5,azerr_5,zcerr_5,vcerr_5,rcerr_5,azcerr_5=cal_zcen_zrms(az[ind],zz[ind],vv[ind],rr[ind],mass[ind],bii)
	az_5, zz_5, mass_5, zzz_15, mm_15, rr_15, com_15 = [], [], [], [], [], [], []
	az_5 = az[ind]
	mass_5 = mass[ind]
	zz_5 = zz[ind]
	rr_15 = rr[ind]
	com_15 = com[ind]
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
#    for j in range(len(zz_5)):
#    	print('%.4f\t%.2f'%(zzz_15[j],mm_15[j]))        ##python warp.py > 15.zm
#    	if az_5[j]>0:
#    	     print('%.2f\t%.2f\t%.3f\t%.2f\t%.3f'%(rr_15[j],az_5[j],zz_5[j],mass_5[j],com_15[j]))
#    	else:
#    	     print('%.2f\t%.2f\t%.3f\t%.2f\t%.3f'%(rr_15[j],az_5[j]+360.,zz_5[j],mass_5[j],com_15[j]))     ###python warp.py > 15.razmc
	pl.scatter(az[ind],zz[ind],s=all_dist[ind],c=fcolor,alpha=0.3,edgecolors=edgecolor)
	pl.plot(thea,w,ls=sty_hi,lw=1.5,color=col_hi)
	pl.plot(thea,zw1,ls=sty_cepheid,lw=1.5,color=col_cepheid)
#    pl.plot(thea,zw5,'-',lw=1.5,color='red')
	pl.plot(thea,zw4,ls=sty_co,lw=1.5,color=col_co)
	pl.plot(thea_15,zcen_15,'-',color=col_cen, zorder=1)
	pl.text(-25,1.3,'%s' %(gr),fontsize=11,fontweight='bold',color=col_text)
	pl.xlim([-30,165])
	pl.ylim([-0.9,1.55])
	pl.yticks([-0.5,0,0.5,1.0,1.5],['','','','',''],fontsize=11,fontweight=1.8)
	pl.xticks([0,50,100,150],['','','',''],fontsize=11,fontweight=1.8)
	pl.twiny()
	pl.xlim(0,radius)
	pl.xticks([0,20,40],['0','20','40'],fontsize=10.5,fontweight=1.8)    
	gr=16
	pl.subplot(2,5,8)
	radius = 2*np.pi*gr*(np.abs(az_min-az_max)/360.)
	thea,w = cal_warp(gr)
	thea,zw1,zw2,zw4,zw5=cal_warpc(gr)
	ind = np.where((rr>=(gr-0.5)) & (rr<(gr+0.5)))[0]
	zcen_6,zrms_6,vcen_6,vrms_6,rcen_6,rrms_6,azcen_6,azrms_6,zerr_6,verr_6,rerr_6,azerr_6,zcerr_6,vcerr_6,rcerr_6,azcerr_6=cal_zcen_zrms(az[ind],zz[ind],vv[ind],rr[ind],mass[ind],bii)
	az_6, zz_6, mass_6, zzz_16, mm_16, rr_16, com_16 = [], [], [], [], [], [], []
	az_6 = az[ind]
	mass_6 = mass[ind]
	zz_6 = zz[ind]
	rr_16= rr[ind]
	com_16 = com[ind]
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
#    for j in range(len(zz_6)):
#    	print('%.4f\t%.2f'%(zzz_16[j],mm_16[j]))        ##python warp.py > 16.zm
#    	if az_6[j]>0:
#    	     print('%.2f\t%.2f\t%.3f\t%.2f\t%.3f'%(rr_16[j],az_6[j],zz_6[j],mass_6[j],com_16[j]))
#    	else:
#    	     print('%.2f\t%.2f\t%.3f\t%.2f\t%.3f'%(rr_16[j],az_6[j]+360.,zz_6[j],mass_6[j],com_16[j]))     ###python warp.py 16.razmc
	pl.scatter(az[ind],zz[ind],s=all_dist[ind],c=fcolor,alpha=0.3,edgecolors=edgecolor)
	pl.plot(thea,w,ls=sty_hi,lw=1.5,color=col_hi)
	pl.plot(thea,zw1,ls=sty_cepheid,lw=1.5,color=col_cepheid)
#    pl.plot(thea,zw5,'-',lw=1.5,color='red')
	pl.plot(thea,zw4,ls=sty_co,lw=1.5,color=col_co)
	pl.plot(thea_16,zcen_16,'-',color=col_cen, zorder=1)
	pl.text(-25,1.3,'%s' %(gr),fontsize=11,fontweight='bold',color=col_text)
	pl.xlim([-30,165])
	pl.ylim([-0.9,1.55])
	pl.yticks([-0.5,0,0.5,1.0,1.5],['','','','',''],fontsize=11,fontweight=1.8)
	pl.xticks([0,50,100,150],['','','',''],fontsize=11,fontweight=1.8)
	pl.twiny()
	pl.xlim(0,radius)
	pl.xticks([0,20,40],['0','20','40'],fontsize=10.5,fontweight=1.8)
	##########
	gr=17.5
	pl.subplot(2,5,9)
	radius = 2*np.pi*gr*(np.abs(az_min-az_max)/360.)
	thea,w = cal_warp(gr)
	thea,zw1,zw2,zw4,zw5=cal_warpc(gr)
	ind = np.where((rr>=(gr-1)) & (rr<(gr+1)))[0]
#    ind = np.where(rr>=(gr-0.5))[0]
	zcen_7,zrms_7,vcen_7,vrms_7,rcen_7,rrms_7,azcen_7,azrms_7,zerr_7,verr_7,rerr_7,azerr_7,zcerr_7,vcerr_7,rcerr_7,azcerr_7=cal_zcen_zrms(az[ind],zz[ind],vv[ind],rr[ind],mass[ind],bii)
	az_7, zz_7, mass_7, zzz_17, mm_17, rr_17, com_17 = [], [], [], [], [], [], []
	az_7 = az[ind]
	mass_7 = mass[ind]
	zz_7 = zz[ind]
	rr_17 = rr[ind]
	com_17 = com[ind]
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
#    for j in range(len(zz_7)):
##    	print('%.4f\t%.2f'%(zzz_17[j],mm_17[j]))        ##python warp.py > 17.zm    
#    	if az_7[j]>0:
#    	     print('%.2f\t%.2f\t%.3f\t%.2f\t%.3f'%(rr_17[j],az_7[j],zz_7[j],mass_7[j],com_17[j]))
#    	else:
#    	     print('%.2f\t%.2f\t%.3f\t%.2f\t%.3f'%(rr_17[j],az_7[j]+360.,zz_7[j],mass_7[j],com_17[j]))     ###python warp.py > 17.razmc
	pl.scatter(az[ind],zz[ind],s=all_dist[ind],c=fcolor,alpha=0.3,edgecolors=edgecolor)
	pl.plot(thea,w,ls=sty_hi,lw=1.5,color=col_hi)
	pl.plot(thea,zw1,ls=sty_cepheid,lw=1.5,color=col_cepheid)
#    pl.plot(thea,zw5,'-',lw=1.5,color='red')
	pl.plot(thea,zw4,ls=sty_co,lw=1.5,color=col_co)
	pl.plot(thea_17,zcen_17,'-',color=col_cen, zorder=1)
	pl.text(-25,1.3,'%s' %(gr),fontsize=11,fontweight='bold',color=col_text)
	pl.xlim([-30,165])
	pl.ylim([-0.9,1.55])
	pl.yticks([-0.5,0,0.5,1.0,1.5],['','','','',''],fontsize=11,fontweight=1.8)
	pl.xticks(fontsize=10.5,fontweight=1.8)
	pl.xticks([0,50,100,150],['','','',''],fontsize=11,fontweight=1.8)
	pl.twiny()
	pl.xlim(0,radius)
	pl.xticks([0,20,40],['0','20','40'],fontsize=10.5,fontweight=1.8)
	gr=19.5
	pl.subplot(2,5,10)
	radius = 2*np.pi*gr*(np.abs(az_min-az_max)/360.)
	thea,w = cal_warp(gr)
	thea,zw1,zw2,zw4,zw5=cal_warpc(gr)
	ind = np.where((rr>=(gr-1)) & (rr<(gr+1)))[0]
	zcen_8,zrms_8,vcen_8,vrms_8,rcen_8,rrms_8,azcen_8,azrms_8,zerr_8,verr_8,rerr_8,azerr_8,zcerr_8,vcerr_8,rcerr_8,azcerr_8=cal_zcen_zrms(az[ind],zz[ind],vv[ind],rr[ind],mass[ind],bii)
	az_8, zz_8, mass_8, zzz_18, mm_18, rr_18, com_18 = [], [], [], [], [], [], []
	az_8 = az[ind]
	mass_8 = mass[ind]
	zz_8 = zz[ind]
	rr_18 = rr[ind]
	com_18 = com[ind]
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
#    for j in range(len(zz_8)):
#    	print('%.4f\t%.2f'%(zzz_18[j],mm_18[j]))        ##python warp.py > 18.zm    
#    	if az_8[j]>0:
#    	     print('%.2f\t%.2f\t%.3f\t%.2f\t%.3f'%(rr_18[j],az_8[j],zz_8[j],mass_8[j],com_18[j]))
#    	else:
#    	     print('%.2f\t%.2f\t%.3f\t%.2f\t%.3f'%(rr_18[j],az_8[j]+360.,zz_8[j],mass_8[j],com_18[j]))     ###python warp.py 18.razmc
	pl.scatter(az[ind],zz[ind],s=all_dist[ind],c=fcolor,alpha=0.3,edgecolors=edgecolor)
	pl.plot(thea,w,ls=sty_hi,lw=1.5,color=col_hi,label='HI')
	pl.plot(thea,zw1,ls=sty_cepheid,lw=1.5,color=col_cepheid,label='Cepheids')
#    pl.plot(thea,zw5,'-',lw=1.5,color='red',label='CO 2comp')
	pl.plot(thea,zw4,ls=sty_co,lw=1.5,color=col_co,label='CO 1comp xf')
	pl.plot(thea_18,zcen_18,'-',color=col_cen, zorder=1)
	pl.text(-25,1.3,'%s' %(gr),fontsize=11,fontweight='bold',color=col_text)
	pl.xlim([-30,165])
	pl.ylim([-0.9,1.55])
	pl.yticks([-0.5,0,0.5,1.0,1.5],['','','','',''],fontsize=11,fontweight=1.8)
	pl.xticks([0,50,100,150],['','','',''],fontsize=11,fontweight=1.8)
	pl.legend() 
	pl.twiny()
	pl.xlim(0,radius)
	pl.xticks([0,20,40,60],['0','20','40','60 kpc'],fontsize=10.5,fontweight=1.8)     
	fig1.savefig('az_z_warp_xf.png',format='png',bbox_inches='tight')
	pl.show()
'''