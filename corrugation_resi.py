### plot Az-(Z-warp model) of MCs in different Rgal

import numpy as np
import pylab as pl
import scipy.interpolate as sp_interp
from scipy.stats import norm
import scipy.optimize as sp_opt
from shared import *
import matplotlib as mpl
import matplotlib.pyplot as plt


comp = 2

if __name__ == '__main__':
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
	#mass = (np.vstack((osc_data*3, out_data*1.5, per_data)).T)[7]
	com = data[8]

	print(az.min(), az.max())
	scatter_size, scatter_color = mass_distribute(mass)
	gr_axis = np.linspace(7.5, 22, 200)

	fig1, ax = plt.subplots(nrows=3, ncols=5, sharex=True, sharey=True, figsize=[17,10])
	plt.subplots_adjust(wspace=0, hspace=0)
	ax = ax.ravel()

	theta_sep = [155, 143, 131, 119,107, 95, 83, 71, 59, 47, 38, 29, 20, 11, -11, -17, -27]
	sep = -1
	for i in range(len(theta_sep)-2):
		### skip 11 ~ -11
		if i == 13: sep += 2
		else: sep += 1
		theta1, theta2 = theta_sep[sep:sep+2]
		theta = (theta1+theta2)/2

		### filter theta
		idx = (az >= theta2) & (az < theta1)

		### plot cloud
		ax[i].scatter(rr[idx], zz[idx], s=scatter_size[idx]*2, **rad_kws_mc)

		### plot binned average and errorbar
		rcen, rrms, zcen, zrms = cal_zcen_zrms(rr[idx], zz[idx], weights=mass[idx], binsize=-1, nbin=30, bin0=8.15 if i<5 else 8)
		ax[i].errorbar(rcen, zcen, yerr=zrms*2.355/2., **rad_kws_err)
		#ax[i].plot(rcen, zcen, '-', color=darker_hex(col_mc), lw=2, alpha=0.8, zorder=10)

		### plot h line
		if comp==1:
			rad_kws_co1['zorder']=0
			ax[i].plot([0, 30], [0, 0], **rad_kws_co1)
		else:
			rad_kws_co2['zorder']=0
			ax[i].plot([0, 30], [0, 0], **rad_kws_co2)

		### plot text
		rad_kws_text['fontsize']=14
		text = r'$\mathbf{\phi_{%i}}$' % (i+1)#=[%s$^{\circ}$, %s$^{\circ}$]' % (i+1, theta1, theta2)
		ax[i].text(0.02, 0.98, text, transform=ax[i].transAxes, **rad_kws_text)

		### ticks
		ax[i].set_xticks(np.arange(0, 30, 4))
		ax[i].set_yticks(np.arange(-3, 3, 0.5))
		ax[i].minorticks_on()
		ax[i].tick_params(axis='x', which='major', labelsize=12)#, direction='in')
		ax[i].tick_params(axis='y', which='major', labelsize=12)#, direction='in')
		ax[i].tick_params(top=True, right=True, direction='in')#, labelsize=1000/self.dpi)
		ax[i].tick_params(which='minor', top=True, right=True, direction='in')#, labelsize=1000/self.dpi)


		if i == 10:
			### sin component
			raxis = np.linspace(0,30,200)
			phiaxis = np.repeat(34, 200)
			zaxis = function_warp((raxis, phiaxis), p=[0,0,0,0,0], sin=p_sin1comp if comp==1 else p_sin2comp)
			ax[i].plot(raxis, zaxis, **rad_kws_sin)
			### arm labels
			ax[i].plot([9.6, 12.6, 17.1], [-0.4, -0.2, -0.25/comp], linestyle='None', color='grey', marker=r'$\uparrow$', ms=8, zorder=30)
			for x,y,t in zip([9.6, 12.6, 17.1], [-0.4, -0.2, -0.25/comp], ['   Perseus', '   Outer', 'OSC']):
				ax[i].text(x,y-0.1,t, ha='center', va='top', fontsize=12)
		if i == 10:
			### label
			ax[i].set_xlabel('R (kpc)', fontsize=13, fontweight='bold')
			ax[i].set_ylabel('Residuals in Z (kpc)', fontsize=13, fontweight='bold')
		if i>=10:
			xtk = ax[i].get_xticks().astype(str)
			xtk[2] = '  '+xtk[2]	#shift '8' a little right
			ax[i].set_xticklabels(xtk) # add kpc at the end
		#if i == 14: ax[i].legend()

		ax[i].set_xlim(8, 21.5)
		ax[i].set_ylim(-0.8, 0.8)

	ax[0].text(-0.25, 0.9, 'b', color='black', font=dict(size=36, family="Arial Black"), transform=ax[0].transAxes)

	plt.savefig('fig/r_z_corrugation_resi_%icomp.png' % comp,format='png',bbox_inches='tight', dpi=400)
	plt.show()

'''
import numpy as np
import pylab as pl
import scipy.interpolate as sp_interp
from scipy.stats import norm
import scipy.optimize as sp_opt
def fp2p(free_params):
	#free params -> mixed params with free and fixed
	
	return params[r'$a_0$'][0], params[r'$a_1$'][0], params[r'$R_{w0}$'][0], params[r'$b_{w0}$'], params[r'$\phi_{w0}$'], \
		free_params
	
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
def cal_flaring(rg):
	b0=0.15
	r0=9.8
	br=b0*np.exp((rg-8.15)/r0)
	return(br)
def cal_warpc(thea):
	p1=[9.26,17.4,0.148]
	p2=[7.72,17.5,0.06,1.33]
	p3=[-2.97397887e-02,1.72161063e-02,-1.66075673e+01,7.93265614e-03,6.05519348e+01]
	p4=[-0.07623866250085598, 0.10786410974035177, 7.799981504630573, 0.9349929777917091, -9.018240871527818] ##1st notcorrect
	#p5=[0.1356796,7.902843844516651,2.6634521838919905,0.16057916214431697,30.674799398760037,8.759276457793465,0.010257305565027436] #a3,Rsin, Period0, Period1, PHIsin, PHIsinwid
	p5=[0.13407291765265944, 8.035984791800805, 2.4779148488734806, 0.17532946235760719, 31.1201876289705, 8.219785507458754]
	p6=[0.11612725147813213, 7.988470812105598, 2.2615529572158253, 0.18997888415358558, 23.75, 5.65]
	p6=[0.11530288482575485, 7.915141586910033, 2.33913108001671, 0.18318029951336756, 23.75, 4.65]
	#p5=[0.1268542768826838, 6.017271758222874, 4.817207909592683, 0.13247450817047815, 144.5, 8.65] #a3, Rsin, Period0, Period1, PHIsin, PHIsinwid
	p6=[0.155, 7.75, 3.0, 0.138, 32.1, 8.28, 0.0145] #a3, Rsin, Period0, Period1, PHIsin, PHIsinwid, a4
	p6=[1.39573585e-01, 7.74237334, 2.86256551, 1.57132833e-01,3.21995180e+01, 9.11163275, 1.06050455e-02]
	p6= [1.77364170e-01,8.11971549,2.28645802,1.91232869e-01,2.89850422e+01,3.99999462,-1.33220676e-03] #a3, Rsin, Period0, Period1, PHIsin, PHIsinwid, a4
	rg, zw1, zw2, zw3, w= [], [], [], [], []
	rg = np.linspace(9,22,100)
	zw1=p1[2]*(rg-p1[0])*np.sin(np.deg2rad(thea-p1[1]))        #chen2019
	zw2=p2[2]*((rg-p2[0])**p2[3])*np.sin(np.deg2rad(thea-p2[1])) #chen2019
	zw2=-p3[0]+(rg - 8)**2*(p3[1]*np.sin(np.deg2rad(thea-p3[2]))+p3[3]*np.sin(2.0*(thea-p3[4])/180.*np.pi))  #mwisp
	zw3=p4[0]+p4[1] * (rg-p4[2])**p4[3] * np.sin((thea-p4[4])/180*np.pi) #mwisp first component
	zw5=np.zeros(rg.size)
	index = rg>p5[1]
	zw5 += p5[0]*np.sin( (rg[index] - p5[1]) / (p5[2] + p5[3]*(rg[index]-p5[1])) * 2*np.pi ) * np.exp(-(thea-p5[4])**2/ 2 / p5[5]**2)
	zw6=np.zeros(rg.size)
	index = rg>p6[1]
	#zw6 += p6[0]*np.sin( (rg[index] - p6[1]) / (p6[2] + p6[3]*(rg[index]-p6[1])) * 2*np.pi ) * np.exp(-(thea-p6[4])**2/ 2 / p6[5]**2)
	zw6 += ( p6[6]*(rg[index] - p6[1]) + p6[0]*np.sin( (rg[index] - p6[1]) / (p6[2] + p6[3]*(rg[index]-p6[1])) * 2*np.pi ) ) * np.exp(-(thea-p6[4])**2/ 2 / p6[5]**2)
	w1p=[9,197,10,-3.1] # [k0,k1,rk1,k2]                      ##levenine 2006
	w0p=[-66,150,15,-0.47] # [k0,k1,rk1,k2]
	w2p=[-70,171,15,-5.3] # [k0,k1,rk1,k2]
	for i in range(len(rg)):                                       
		w0=w0p[0]+w0p[1]*(rg[i]-w0p[2])+w0p[3]*(rg[i]-w0p[2])**2
		w1=w1p[0]+w1p[1]*(rg[i]-w1p[2])+w1p[3]*(rg[i]-w1p[2])**2
		w2=w2p[0]+w2p[1]*(rg[i]-w2p[2])+w2p[3]*(rg[i]-w2p[2])**2
		if ((w2 >= 150) & (rg[i] > 15)):
			temp_w= (w0+w1*np.sin(np.deg2rad(thea))+w2*np.sin(np.deg2rad(2.*thea)))/1000.
			w.append(temp_w)
		elif ((w2 < 150) & (rg[i] > 15)):
			temp_w= (w0+w1*np.sin(np.deg2rad(thea)))/1000.
			w.append(temp_w)
		else:
			temp_w= (w1*np.sin(np.deg2rad(thea)))/1000.
			w.append(temp_w)
	return(rg,zw1,zw2,zw3,zw5,zw6,w)
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
def cal_zcen_zrms(az,zz,vv,rr,mass,r1,r2,thea1,thea2):
    zcen,zrms,vcen,vrms,rcen,rrms,r,z_err,v_err,r_err,zcen_err,vcen_err,rcen_err,xcen,ycen = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    ind = np.where((az>thea2) & (az<thea1) & (rr>=r1) & (rr<r2))
    if (len(np.transpose(ind))>2):
    	zcen,zrms = weighted_avg_and_std(zz[ind],mass[ind])
    	vcen,vrms = weighted_avg_and_std(vv[ind],mass[ind])
    	rcen,rrms = weighted_avg_and_std(rr[ind],mass[ind])
    	r = (r1+r2)/2.
    	thea = (thea1+thea2)/2.
    	z_err = zrms/np.sqrt(len(np.transpose(ind)))
    	v_err = vrms/np.sqrt(len(np.transpose(ind)))
    	r_err = rrms/np.sqrt(len(np.transpose(ind)))
    	zcen_err = zcen/np.sqrt(len(np.transpose(ind)))
    	vcen_err = vcen/np.sqrt(len(np.transpose(ind)))
    	rcen_err = rcen/np.sqrt(len(np.transpose(ind)))
    	xcen = rcen*np.sin(np.deg2rad(thea))
    	ycen = rcen*np.cos(np.deg2rad(thea))
    return(zcen,zrms,vcen,vrms,rcen,rrms,r,z_err,v_err,r_err,zcen_err,vcen_err,rcen_err,xcen,ycen,ind)
if __name__ == '__main__':

	out_data = np.loadtxt('out_para.txt',comments='#')
	osc_data = np.loadtxt('osc2_para.txt',comments='#')
	per_data = np.loadtxt('per_para.txt',comments='#')
	residual = np.loadtxt('residual_1com.razm')
	col,col1,col_hi,col_cepheid,col_co,col_cen,col_err = ['teal','none','olive','blue','crimson','black','#36454F']

	out_mass = out_data[:,7]#/out_data[:,8]
	out_rr = out_data[:,11]
	out_zz = out_data[:,14]#*1000
	out_az = out_data[:,19]
	out_dd = out_data[:,4]
	out_vv = out_data[:,2]
	out_bb = out_data[:,1]
	out_ll = out_data[:,0]
	osc_mass = osc_data[:,7]#/osc_data[:,8]
	osc_rr = osc_data[:,11]
	osc_zz = osc_data[:,14]#*1000.
	osc_az = osc_data[:,19]
	osc_dd = osc_data[:,4]
	osc_vv = osc_data[:,2]
	osc_bb = osc_data[:,1]
	osc_ll = osc_data[:,0]
	per_mass = per_data[:,7]#/per_data[:,8]
	per_rr = per_data[:,11]
	per_zz = per_data[:,14]#*1000.
	per_az = per_data[:,19]
	per_dd = per_data[:,4]
	per_vv = per_data[:,2]
	per_bb = per_data[:,1]
	per_ll = per_data[:,0]
	best=[-0.08225466502024195,0.10192102165531088, 7.779011564372266, 1.0, -11.000005387066379, -0.11935558287557957, 13.791351404542997, 1.0, -23.76326592401999, 0.0, 8.46, 1.228, 0.476, 29.83, 10.0, 0.0007]

#	az = np.loadtxt('az_out2.txt')   ##9 to 19kpc
	#az = np.hstack((osc_az,out_az,per_az))
	l = np.hstack((osc_ll,out_ll,per_ll))
	vv = np.hstack((osc_vv,out_vv,per_vv))
	#rr = np.hstack((osc_rr,out_rr,per_rr))
	
	rr = residual[:,0]
	az = residual[:,1]
	zz = residual[:,2] #zz-function((rr, az), best) ##residual[:,2]
	mass = residual[:,3]
	mas = np.hstack((osc_mass,out_mass,per_mass))
	ind = np.where((mas > 0) & ~((l > 160) & (l <= 165) | (l >= 195) & (l < 200)))[0]
	dd = np.hstack((osc_dd,out_dd,per_dd))
	vv = vv[ind]	
	all_mass = np.hstack((osc_mass*1,out_mass*1.,per_mass))
	all_mass = all_mass[ind]    #####################################
	all_dist,all_colors = mass_distribute(all_mass)
	rho_1,rherr_1,rcen_1,zcen_1,zrms_1,ind1 = [],[],[],[],[],[]
	fig1, ax = pl.subplots(sharex=True, figsize=[17,10.2])
#	fig1.tight_layout()
	pl.subplots_adjust(wspace=0,hspace=0)
	pl.subplot(3,5,1)
	#155_145
	thea=[160,148]
	rg,zw1,zw2,zw3,zw5,zw6,w=cal_warpc((thea[0]+thea[1])/2.)
	for i in range(11):
		zcen,zrms,vcen,vrms,rcen,rrms,radius,zerr,verr,rerr,zcerr,vcerr,rcerr,xcen,ycen,ind=cal_zcen_zrms(az,zz,vv,rr,mass,8.15+i,8.15+(i+1),thea[0],thea[1])
		if (len(np.transpose(ind))>0):
			pl.scatter(rr[ind],zz[ind],s=all_dist[ind],alpha=0.35,color=col,edgecolors=col1)
		if np.isnan(zrms) == False:
			rho = 18.525*vrms**2/(zrms)**2
			rherr = np.sqrt(4.*(verr/vrms)**2+4.*(zerr/zrms)**2)
			rho_1.append(rho)
			rherr_1.append(rherr)
			rcen_1.append(rcen)
			zcen_1.append(zcen)
			zrms_1.append(zrms)
			pl.errorbar(rcen,zcen,yerr=zrms*2.355/2.,fmt='^',c=col_err,markersize=5.,elinewidth = 1, capsize=2)
	pl.text(8.5,0.8,'[%s,%s]' %(thea[0],thea[1]),fontsize=11,fontweight=1.8)
	pl.xlim([8,21.5])
	pl.ylim([-0.9,0.9])
#	pl.plot(rg,w,'--',color='black')
#	pl.plot(rg,zw1,'-',color='blue')
#	pl.plot(rg,zw2,'--',color='red')
#	pl.plot(rg,zw5,'--',color='green')
	pl.plot([8,21.5],[0,0],ls='dotted',color='red',lw=1)
	pl.xticks([8,10,12,14,16,18,20],['','','','','','',''],fontsize=11,fontweight=1.8)
	pl.yticks([-0.8,-0.4,0,0.4,0.8],['','','','',''],fontsize=11,fontweight=1.8)
	#145_135
	thea=[148,136]
	rg,zw1,zw2,zw3,zw5,zw6,w=cal_warpc((thea[0]+thea[1])/2.)
	rho_2,rherr_2,rcen_2,zcen_2,zrms_2,ind2 = [],[],[],[],[],[]
	pl.subplot(3,5,2)
	for i in range(6):
		zcen,zrms,vcen,vrms,rcen,rrms,radius,zerr,verr,rerr,zcerr,vcerr,rcerr,xcen,ycen,ind=cal_zcen_zrms(az,zz,vv,rr,mass,8.15+i,8.15+(i+1),thea[0],thea[1])
		if (len(np.transpose(ind))>0):
			pl.scatter(rr[ind],zz[ind],s=all_dist[ind],alpha=0.35,color=col,edgecolors=col1)
		if np.isnan(zrms) == False:
			rho = 18.525*vrms**2/(zrms)**2
			rherr = np.sqrt(4.*(verr/vrms)**2+4.*(zerr/zrms)**2)
			rho_2.append(rho)
			rherr_2.append(rherr)
			rcen_2.append(rcen)
			zcen_2.append(zcen)
			zrms_2.append(zrms)
			pl.errorbar(rcen,zcen,yerr=zrms*2.355/2.,fmt='^',c=col_err,markersize=5.,elinewidth = 1, capsize=2)
	pl.text(8.5,0.8,'[%s,%s]' %(thea[0],thea[1]),fontsize=11,fontweight=1.8)
	pl.xlim([8,21.5])
	pl.ylim([-0.9,0.9])
#	pl.plot(rg,w,'--',color='black')
#	pl.plot(rg,zw1,'-',color='blue')
#	pl.plot(rg,zw2,'--',color='red')
#	pl.plot(rg,zw5,'--',color='green')
	pl.plot([8,21.5],[0,0],ls='dotted',color='red',lw=1)
	pl.xticks([8,10,12,14,16,18,20],['8','10','12','14','16','18','20'],fontsize=11,fontweight=1.8)
	pl.yticks([-0.8,-0.4,0,0.4,0.8],['','','','',''],fontsize=11,fontweight=1.8)
	#135_125
	thea=[136,124]
	rg,zw1,zw2,zw3,zw5,zw6,w=cal_warpc((thea[0]+thea[1])/2.)
	rho_3,rherr_3,rcen_3,zcen_3,zrms_3,ind3 = [],[],[],[],[],[]
	pl.subplot(3,5,3)
	for i in range(11):
		zcen,zrms,vcen,vrms,rcen,rrms,radius,zerr,verr,rerr,zcerr,vcerr,rcerr,xcen,ycen,ind=cal_zcen_zrms(az,zz,vv,rr,mass,8.15+i,8.15+(i+1),thea[0],thea[1])
		if (len(np.transpose(ind))>0):
			pl.scatter(rr[ind],zz[ind],s=all_dist[ind],alpha=0.35,color=col,edgecolors=col1)
		if np.isnan(zrms) == False:
			rho = 18.525*vrms**2/(zrms)**2
			rherr = np.sqrt(4.*(verr/vrms)**2+4.*(zerr/zrms)**2)
			rho_3.append(rho)
			rherr_3.append(rherr)
			rcen_3.append(rcen)
			zcen_3.append(zcen)
			zrms_3.append(zrms)
			pl.errorbar(rcen,zcen,yerr=zrms*2.355/2.,fmt='^',c=col_err,markersize=5.,elinewidth = 1, capsize=2)
	pl.text(8.5,0.8,'[%s,%s]' %(thea[0],thea[1]),fontsize=11,fontweight=1.8)
	pl.xlim([8,21.5])
	pl.ylim([-0.9,0.9])
#	pl.plot(rg,w,'--',color='black')
#	pl.plot(rg,zw1,'-',color='blue')
#	pl.plot(rg,zw2,'--',color='red')
#	pl.plot(rg,zw5,'--',color='green')
	pl.plot([8,21.5],[0,0],ls='dotted',color='red',lw=1)
	pl.xticks([8,10,12,14,16,18,20],['8','10','12','14','16','18','20'],fontsize=11,fontweight=1.8)
	pl.yticks([-0.8,-0.4,0,0.4,0.8],['','','','',''],fontsize=11,fontweight=1.8)
	#125_115
	thea=[124,112]
	rg,zw1,zw2,zw3,zw5,zw6,w=cal_warpc((thea[0]+thea[1])/2.)
	rho_4,rherr_4,rcen_4,zcen_4,zrms_4,ind4 = [],[],[],[],[],[]
	pl.subplot(3,5,4)
	for i in range(11):
		zcen,zrms,vcen,vrms,rcen,rrms,radius,zerr,verr,rerr,zcerr,vcerr,rcerr,xcen,ycen,ind=cal_zcen_zrms(az,zz,vv,rr,mass,8.15+i,8.15+(i+1),thea[0],thea[1])
		if (len(np.transpose(ind))>0):
			pl.scatter(rr[ind],zz[ind],s=all_dist[ind],alpha=0.35,color=col,edgecolors=col1)
		if np.isnan(zrms) == False:
			rho = 18.525*vrms**2/(zrms)**2
			rherr = np.sqrt(4.*(verr/vrms)**2+4.*(zerr/zrms)**2)
			rho_4.append(rho)
			rherr_4.append(rherr)
			rcen_4.append(rcen)
			zcen_4.append(zcen)
			zrms_4.append(zrms)
			pl.errorbar(rcen,zcen,yerr=zrms*2.355/2.,fmt='^',c=col_err,markersize=5.,elinewidth = 1, capsize=2)
	pl.text(8.5,0.8,'[%s,%s]' %(thea[0],thea[1]),fontsize=11,fontweight=1.8)
	pl.xlim([8,21.5])
	pl.ylim([-0.9,0.9])
#	pl.plot(rg,w,'--',color='black')
#	pl.plot(rg,zw1,'-',color='blue')
#	pl.plot(rg,zw2,'--',color='red')
#	pl.plot(rg,zw5,'--',color='green')
	pl.plot([8,21.5],[0,0],ls='dotted',color='red',lw=1)
	pl.xticks([8,10,12,14,16,18,20],['8','10','12','14','16','18','20'],fontsize=11,fontweight=1.8)
	pl.yticks([-0.8,-0.4,0,0.4,0.8],['','','','',''],fontsize=11,fontweight=1.8)
	#115_105
	thea=[112,100]
	rg,zw1,zw2,zw3,zw5,zw6,w=cal_warpc((thea[0]+thea[1])/2.)
	rho_5,rherr_5,rcen_5,zcen_5,zrms_5,ind5 = [],[],[],[],[],[]
	pl.subplot(3,5,5)
	for i in range(11):
		zcen,zrms,vcen,vrms,rcen,rrms,radius,zerr,verr,rerr,zcerr,vcerr,rcerr,xcen,ycen,ind=cal_zcen_zrms(az,zz,vv,rr,mass,8.15+i,8.15+(i+1),thea[0],thea[1])
		if (len(np.transpose(ind))>0):
			pl.scatter(rr[ind],zz[ind],s=all_dist[ind],alpha=0.35,color=col,edgecolors=col1)
		if np.isnan(zrms) == False:
			rho = 18.525*vrms**2/(zrms)**2
			rherr = np.sqrt(4.*(verr/vrms)**2+4.*(zerr/zrms)**2)
			rho_5.append(rho)
			rherr_5.append(rherr)
			rcen_5.append(rcen)
			zcen_5.append(zcen)
			zrms_5.append(zrms)
			pl.errorbar(rcen,zcen,yerr=zrms*2.355/2.,fmt='^',c=col_err,markersize=5.,elinewidth = 1, capsize=2)
	pl.text(8.5,0.8,'[%s,%s]' %(thea[0],thea[1]),fontsize=11,fontweight=1.8)
	pl.xlim([8,21.5])
	pl.ylim([-0.9,0.9])
#	pl.plot(rg,w,'--',color='black')
#	pl.plot(rg,zw1,'-',color='blue')
#	pl.plot(rg,zw2,'--',color='red')
#	pl.plot(rg,zw5,'--',color='green')
	pl.plot([8,21.5],[0,0],ls='dotted',color='red',lw=1)
	pl.xticks([8,10,12,14,16,18,20],['','','','','','',''],fontsize=11,fontweight=1.8)
	pl.yticks([-0.8,-0.4,0,0.4,0.8],['','','','',''],fontsize=11,fontweight=1.8)	
	#105_95
	thea=[100,88]
	rg,zw1,zw2,zw3,zw5,zw6,w=cal_warpc((thea[0]+thea[1])/2.)
	rho_6,rherr_6,rcen_6,zcen_6,zrms_6,ind6 = [],[],[],[],[],[]
	pl.subplot(3,5,6)
	for i in range(11):
		zcen,zrms,vcen,vrms,rcen,rrms,radius,zerr,verr,rerr,zcerr,vcerr,rcerr,xcen,ycen,ind=cal_zcen_zrms(az,zz,vv,rr,mass,8.+i,8.+(i+1),thea[0],thea[1])
		if (len(np.transpose(ind))>0):
			pl.scatter(rr[ind],zz[ind],s=all_dist[ind],alpha=0.35,color=col,edgecolors=col1)
		if np.isnan(zrms) == False:
			rho = 18.525*vrms**2/(zrms)**2
			rherr = np.sqrt(4.*(verr/vrms)**2+4.*(zerr/zrms)**2)
			rho_6.append(rho)
			rherr_6.append(rherr)
			rcen_6.append(rcen)
			zcen_6.append(zcen)
			zrms_6.append(zrms)
			pl.errorbar(rcen,zcen,yerr=zrms*2.355/2.,fmt='^',c=col_err,markersize=5.,elinewidth = 1, capsize=2)
	pl.text(8.5,0.8,'[%s,%s]' %(thea[0],thea[1]),fontsize=11,fontweight=1.8)
	pl.xlim([8,21.5])
	pl.ylim([-0.9,0.9])
#	pl.plot(rg,w,'--',color='black')
#	pl.plot(rg,zw1,'-',color='blue')
#	pl.plot(rg,zw2,'--',color='red')
#	pl.plot(rg,zw5,'--',color='green')
	pl.plot([8,21.5],[0,0],ls='dotted',color='red',lw=1)
	pl.xticks([8,10,12,14,16,18,20],['','','','','','',''],fontsize=11,fontweight=1.8)
	pl.yticks([-0.8,-0.4,0,0.4,0.8],['','','','',''],fontsize=11,fontweight=1.8)	
	#95_85
	thea=[88,76]
	rg,zw1,zw2,zw3,zw5,zw6,w=cal_warpc((thea[0]+thea[1])/2.)
	rho_7,rherr_7,rcen_7,zcen_7,zrms_7,ind7 = [],[],[],[],[],[]
	pl.subplot(3,5,7)
	for i in range(11):
		zcen,zrms,vcen,vrms,rcen,rrms,radius,zerr,verr,rerr,zcerr,vcerr,rcerr,xcen,ycen,ind=cal_zcen_zrms(az,zz,vv,rr,mass,8.+i,8.+(i+1),thea[0],thea[1])
		if (len(np.transpose(ind))>0):
			pl.scatter(rr[ind],zz[ind],s=all_dist[ind],alpha=0.35,color=col,edgecolors=col1)
		if np.isnan(zrms) == False:
			rho = 18.525*vrms**2/(zrms)**2
			rherr = np.sqrt(4.*(verr/vrms)**2+4.*(zerr/zrms)**2)
			rho_7.append(rho)
			rherr_7.append(rherr)
			rcen_7.append(rcen)
			zcen_7.append(zcen)
			zrms_7.append(zrms)
			pl.errorbar(rcen,zcen,yerr=zrms*2.355/2.,fmt='^',c=col_err,markersize=5.,elinewidth = 1, capsize=2)
	pl.text(8.5,0.8,'[%s,%s]' %(thea[0],thea[1]),fontsize=11,fontweight=1.8)
	pl.xlim([8,21.5])
	pl.ylim([-0.9,0.9])
#	pl.plot(rg,w,'--',color='black')
#	pl.plot(rg,zw1,'-',color='blue')
#	pl.plot(rg,zw2,'--',color='red')
#	pl.plot(rg,zw5,'--',color='green')
	pl.plot([8,21.5],[0,0],ls='dotted',color='red',lw=1)
	pl.xticks([8,10,12,14,16,18,20],['8','10','12','14','16','18','20'],fontsize=11,fontweight=1.8)
	pl.yticks([-0.8,-0.4,0,0.4,0.8],['','','','',''],fontsize=11,fontweight=1.8)
	#85_75
	thea=[76,64]
	rg,zw1,zw2,zw3,zw5,zw6,w=cal_warpc((thea[0]+thea[1])/2.)
	rho_8,rherr_8,rcen_8,zcen_8,zrms_8,ind8 = [],[],[],[],[],[]
	pl.subplot(3,5,8)
	for i in range(14):
		zcen,zrms,vcen,vrms,rcen,rrms,radius,zerr,verr,rerr,zcerr,vcerr,rcerr,xcen,ycen,ind=cal_zcen_zrms(az,zz,vv,rr,mass,8.+i,8.+(i+1),thea[0],thea[1])
		if (len(np.transpose(ind))>0):
			pl.scatter(rr[ind],zz[ind],s=all_dist[ind],alpha=0.35,color=col,edgecolors=col1)
		if np.isnan(zrms) == False:
			rho = 18.525*vrms**2/(zrms)**2
			rherr = np.sqrt(4.*(verr/vrms)**2+4.*(zerr/zrms)**2)
			rho_8.append(rho)
			rherr_8.append(rherr)
			rcen_8.append(rcen)
			zcen_8.append(zcen)
			zrms_8.append(zrms)
			pl.errorbar(rcen,zcen,yerr=zrms*2.355/2.,fmt='^',c=col_err,markersize=5.,elinewidth = 1, capsize=2)
	pl.text(8.5,0.8,'[%s,%s]' %(thea[0],thea[1]),fontsize=11,fontweight=1.8)
	pl.xlim([8,21.5])
	pl.ylim([-0.9,0.9])
#	pl.plot(rg,w,'--',color='black')
#	pl.plot(rg,zw1,'-',color='blue')
#	pl.plot(rg,zw2,'--',color='red')
	#pl.plot(rg,zw5,'--',color='green')
	pl.plot([8,21.5],[0,0],ls='dotted',color='red',lw=1)
	pl.xticks([8,10,12,14,16,18,20],['','','','','','',''],fontsize=11,fontweight=1.8)
	pl.yticks([-0.8,-0.4,0,0.4,0.8],['','','','',''],fontsize=11,fontweight=1.8)		
	#75_65
	thea=[64,52]
	rg,zw1,zw2,zw3,zw5,zw6,w=cal_warpc((thea[0]+thea[1])/2.)
	rho_9,rherr_9,rcen_9,zcen_9,zrms_9,ind9 = [],[],[],[],[],[]
	pl.subplot(3,5,9)
	for i in range(14):
		zcen,zrms,vcen,vrms,rcen,rrms,radius,zerr,verr,rerr,zcerr,vcerr,rcerr,xcen,ycen,ind=cal_zcen_zrms(az,zz,vv,rr,mass,8.+i,8.+(i+1),thea[0],thea[1])
		if (len(np.transpose(ind))>0):
			pl.scatter(rr[ind],zz[ind],s=all_dist[ind],alpha=0.35,color=col,edgecolors=col1)
		if np.isnan(zrms) == False:
			rho = 18.525*vrms**2/(zrms)**2
			rherr = np.sqrt(4.*(verr/vrms)**2+4.*(zerr/zrms)**2)
			rho_9.append(rho)
			rherr_9.append(rherr)
			rcen_9.append(rcen)
			zcen_9.append(zcen)
			zrms_9.append(zrms)
			pl.errorbar(rcen,zcen,yerr=zrms*2.355/2.,fmt='^',c=col_err,markersize=5.,elinewidth = 1, capsize=2)
	pl.text(8.5,0.8,'[%s,%s]' %(thea[0],thea[1]),fontsize=11,fontweight=1.8)
	pl.xlim([8,21.5])
	pl.ylim([-0.9,0.9])
#	pl.plot(rg,w,'--',color='black')
#	pl.plot(rg,zw1,'-',color='blue')
#	pl.plot(rg,zw2,'--',color='red')
#	pl.plot(rg,zw5,'--',color='green')
	pl.plot([8,21.5],[0,0],ls='dotted',color='red',lw=1)
	pl.xticks([8,10,12,14,16,18,20],['','','','','','',''],fontsize=11,fontweight=1.8)
	pl.yticks([-0.8,-0.4,0,0.4,0.8],['','','','',''],fontsize=11,fontweight=1.8)		
	#65_55
	thea=[52,40]
	rg,zw1,zw2,zw3,zw5,zw6,w=cal_warpc((thea[0]+thea[1])/2.)
	rho_10,rherr_10,rcen_10,zcen_10,zrms_10,ind10 = [],[],[],[],[],[]
	pl.subplot(3,5,10)
	for i in range(14):
		zcen,zrms,vcen,vrms,rcen,rrms,radius,zerr,verr,rerr,zcerr,vcerr,rcerr,xcen,ycen,ind=cal_zcen_zrms(az,zz,vv,rr,mass,8.+i,8.+(i+1),thea[0],thea[1])
		if (len(np.transpose(ind))>0):
			pl.scatter(rr[ind],zz[ind],s=all_dist[ind],alpha=0.35,color=col,edgecolors=col1)
		if np.isnan(zrms) == False:
			rho = 18.525*vrms**2/(zrms)**2
			rherr = np.sqrt(4.*(verr/vrms)**2+4.*(zerr/zrms)**2)
			rho_10.append(rho)
			rherr_10.append(rherr)
			rcen_10.append(rcen)
			zcen_10.append(zcen)
			zrms_10.append(zrms)
			pl.errorbar(rcen,zcen,yerr=zrms*2.355/2.,fmt='^',c=col_err,markersize=5.,elinewidth = 1, capsize=2)
	pl.text(8.5,0.8,'[%s,%s]' %(thea[0],thea[1]),fontsize=11,fontweight=1.8)
	pl.xlim([8,21.5])
	pl.ylim([-0.9,0.9])
#	pl.plot(rg,w,'--',color='black')
#	pl.plot(rg,zw1,'-',color='blue')
#	pl.plot(rg,zw2,'--',color='red')
#	pl.plot(rg,zw5,'--',color='green')
	pl.plot([8,21.5],[0,0],ls='dotted',color='red',lw=1)
	pl.xticks([8,10,12,14,16,18,20],['8','10','12','14','16','18','20'],fontsize=11,fontweight=1.8)
	pl.yticks([-0.8,-0.4,0,0.4,0.8],['','','','',''],fontsize=11,fontweight=1.8)
	#55_45
	thea=[40,28]
	rg,zw1,zw2,zw3,zw5,zw6,w=cal_warpc((thea[0]+thea[1])/2.)
	rho_11,rherr_11,rcen_11,zcen_11,zrms_11,ind11 = [],[],[],[],[],[]
	pl.subplot(3,5,11)
	for i in range(14):
		zcen,zrms,vcen,vrms,rcen,rrms,radius,zerr,verr,rerr,zcerr,vcerr,rcerr,xcen,ycen,ind=cal_zcen_zrms(az,zz,vv,rr,mass,8.+i,8.+(i+1),thea[0],thea[1])
		if (len(np.transpose(ind))>0):
			pl.scatter(rr[ind],zz[ind],s=all_dist[ind],alpha=0.35,color=col,edgecolors=col1)
		if np.isnan(zrms) == False:
			rho = 18.525*vrms**2/(zrms)**2
			rherr = np.sqrt(4.*(verr/vrms)**2+4.*(zerr/zrms)**2)
			rho_11.append(rho)
			rherr_11.append(rherr)
			rcen_11.append(rcen)
			zcen_11.append(zcen)
			zrms_11.append(zrms)
			pl.errorbar(rcen,zcen,yerr=zrms*2.355/2.,fmt='^',c=col_err,markersize=5.,elinewidth = 1, capsize=2)
	pl.text(8.5,0.8,'[%s,%s]' %(thea[0],thea[1]),fontsize=11,fontweight=1.8)
	pl.xlim([8,21.5])
	pl.ylim([-0.9,0.9])
#	pl.plot(rg,w,'--',color='black')
#	pl.plot(rg,zw1,'-',color='blue')
#	pl.plot(rg,zw2,'--',color='red')
	pl.plot(rg,zw6,ls='-',color='orange',lw=1.5)
	pl.plot([8,21.5],[0,0],ls='dotted',color='red',lw=1)
	pl.xticks([8,10,12,14,16,18,20],['8','10','12','14','16','18','20'],fontsize=11,fontweight=1.8)
	pl.yticks([-0.8,-0.4,0,0.4,0.8],['-0.8','-0.4','0','0.4','0.8'],fontsize=11,fontweight=1.8)
	pl.xlabel('R (kpc)',fontsize=13,fontweight='bold')
	pl.ylabel('Z (kpc)',fontsize=13,fontweight='bold')
	#45_35
	thea=[28,16]
	rg,zw1,zw2,zw3,zw5,zw6,w=cal_warpc((thea[0]+thea[1])/2.)
	rho_12,rherr_12,rcen_12,zcen_12,zrms_12,ind12 = [],[],[],[],[],[]
	pl.subplot(3,5,12)
	for i in range(14):
		zcen,zrms,vcen,vrms,rcen,rrms,radius,zerr,verr,rerr,zcerr,vcerr,rcerr,xcen,ycen,ind=cal_zcen_zrms(az,zz,vv,rr,mass,8.+i,8.+(i+1),thea[0],thea[1])
		if (len(np.transpose(ind))>0):
			pl.scatter(rr[ind],zz[ind],s=all_dist[ind],alpha=0.35,color=col,edgecolors=col1)
		if np.isnan(zrms) == False:
			rho = 18.525*vrms**2/(zrms)**2
			rherr = np.sqrt(4.*(verr/vrms)**2+4.*(zerr/zrms)**2)
			rho_12.append(rho)
			rherr_12.append(rherr)
			rcen_12.append(rcen)
			zcen_12.append(zcen)
			zrms_12.append(zrms)
			pl.errorbar(rcen,zcen,yerr=zrms*2.355/2.,fmt='^',c=col_err,markersize=5.,elinewidth = 1, capsize=2)
	pl.text(8.5,0.8,'[%s,%s]' %(thea[0],thea[1]),fontsize=11,fontweight=1.8)
	pl.xlim([8,21.5])
	pl.ylim([-0.9,0.9])
#	pl.plot(rg,w,'--',color='black')
#	pl.plot(rg,zw1,'-',color='blue')
#	pl.plot(rg,zw2,'--',color='red')
#	pl.plot(rg,zw5,'-',color='green')
	pl.plot([8,21.5],[0,0],ls='dotted',color='red',lw=1)
	pl.xticks([8,10,12,14,16,18,20],['8','10','12','14','16','18','20'],fontsize=11,fontweight=1.8)
	pl.yticks([-0.8,-0.4,0,0.4,0.8],['','','','',''],fontsize=11,fontweight=1.8)
	#35_25
	thea=[16,4]
	rg,zw1,zw2,zw3,zw5,zw6,w=cal_warpc((thea[0]+thea[1])/2.)
	rho_13,rherr_13,rcen_13,zcen_13,zrms_13,ind13 = [],[],[],[],[],[]
	pl.subplot(3,5,13)
	for i in range(14):
		zcen,zrms,vcen,vrms,rcen,rrms,radius,zerr,verr,rerr,zcerr,vcerr,rcerr,xcen,ycen,ind=cal_zcen_zrms(az,zz,vv,rr,mass,8.+i,8.+(i+1),thea[0],thea[1])
		if (len(np.transpose(ind))>0):
			pl.scatter(rr[ind],zz[ind],s=all_dist[ind],alpha=0.35,color=col,edgecolors=col1)
		if np.isnan(zrms) == False:
			rho = 18.525*vrms**2/(zrms)**2
			rherr = np.sqrt(4.*(verr/vrms)**2+4.*(zerr/zrms)**2)
			rho_13.append(rho)
			rherr_13.append(rherr)
			rcen_13.append(rcen)
			zcen_13.append(zcen)
			zrms_13.append(zrms)
			pl.errorbar(rcen,zcen,yerr=zrms*2.355/2.,fmt='^',c=col_err,markersize=5.,elinewidth = 1, capsize=2)
	pl.text(8.5,0.8,'[%s,%s]' %(thea[0],thea[1]),fontsize=11,fontweight=1.8)
	pl.xlim([8,21.5])
	pl.ylim([-0.9,0.9])
#	pl.plot(rg,w,'--',color='black')
#	pl.plot(rg,zw1,'-',color='blue')
#	pl.plot(rg,zw2,'--',color='red')
#	pl.plot(rg,zw6,'--',color='green')
	pl.plot([8,21.5],[0,0],ls='dotted',color='red',lw=1)
	pl.xticks([8,10,12,14,16,18,20],['8','10','12','14','16','18','20'],fontsize=11,fontweight=1.8)
	pl.yticks([-0.8,-0.4,0,0.4,0.8],['','','','',''],fontsize=11,fontweight=1.8)	

	#25_15
	thea=[-4,-16]
	rg,zw1,zw2,zw3,zw5,zw6,w=cal_warpc((thea[0]+thea[1])/2.)
	rho_14,rherr_14,rcen_14,zcen_14,zrms_14,ind14 = [],[],[],[],[],[]
	pl.subplot(3,5,14)
	for i in range(14):
		zcen,zrms,vcen,vrms,rcen,rrms,radius,zerr,verr,rerr,zcerr,vcerr,rcerr,xcen,ycen,ind=cal_zcen_zrms(az,zz,vv,rr,mass,8.+i,8.+(i+1),thea[0],thea[1])
		if (len(np.transpose(ind))>0):
			pl.scatter(rr[ind],zz[ind],s=all_dist[ind],alpha=0.35,color=col,edgecolors=col1)
		if np.isnan(zrms) == False:
			rho = 18.525*vrms**2/(zrms)**2
			rherr = np.sqrt(4.*(verr/vrms)**2+4.*(zerr/zrms)**2)
			rho_14.append(rho)
			rherr_14.append(rherr)
			rcen_14.append(rcen)
			zcen_14.append(zcen)
			zrms_14.append(zrms)
			pl.errorbar(rcen,zcen,yerr=zrms*2.355/2.,fmt='^',c=col_err,markersize=5.,elinewidth = 1, capsize=2)
	pl.text(8.5,0.8,'[%s,%s]' %(thea[0],thea[1]),fontsize=11,fontweight=1.8)
	pl.xlim([8,21.5])
	pl.ylim([-0.9,0.9])
#	pl.plot(rg,w,'--',color='black')
#	pl.plot(rg,zw1,'-',color='blue')
#	pl.plot(rg,zw2,'--',color='red')
#	pl.plot(rg,zw5,'--',color='green')
	pl.plot([8,21.5],[0,0],ls='dotted',color='red',lw=1)
	pl.xticks([8,10,12,14,16,18,20],['8','10','12','14','16','18','20'],fontsize=11,fontweight=1.8)
	pl.yticks([-0.8,-0.4,0,0.4,0.8],['','','','',''],fontsize=11,fontweight=1.8)
	
	#15_5
	thea=[-16,-28]
	rg,zw1,zw2,zw3,zw5,zw6,w=cal_warpc((thea[0]+thea[1])/2.)
	rho_15,rherr_15,rcen_15,zcen_15,zrms_15,ind15 = [],[],[],[],[],[]
	pl.subplot(3,5,15)
	for i in range(14):
		zcen,zrms,vcen,vrms,rcen,rrms,radius,zerr,verr,rerr,zcerr,vcerr,rcerr,xcen,ycen,ind=cal_zcen_zrms(az,zz,vv,rr,mass,8.+i,8.+(i+1),thea[0],thea[1])
		if (len(np.transpose(ind))>0):
			pl.scatter(rr[ind],zz[ind],s=all_dist[ind],alpha=0.35,color=col,edgecolors=col1)
		if np.isnan(zrms) == False:
			rho = 18.525*vrms**2/(zrms)**2
			rherr = np.sqrt(4.*(verr/vrms)**2+4.*(zerr/zrms)**2)
			rho_15.append(rho)
			rherr_15.append(rherr)
			rcen_15.append(rcen)
			zcen_15.append(zcen)
			zrms_15.append(zrms)
			pl.errorbar(rcen,zcen,yerr=zrms*2.355/2.,fmt='^',c=col_err,markersize=5.,elinewidth = 1, capsize=2)
	pl.text(8.5,0.8,'[%s,%s]' %(thea[0],thea[1]),fontsize=11,fontweight=1.8)
	pl.xlim([8,21.5])
	pl.ylim([-0.9,0.9])
#	pl.plot(rg,w,'--',color='black')
#	pl.plot(rg,zw1,'-',color='blue')
#	pl.plot(rg,zw2,'--',color='red')
#	pl.plot(rg,zw5,'--',color='green')
	pl.plot([8,21.5],[0,0],ls='dotted',color='red',lw=1)
	pl.xticks([8,10,12,14,16,18,20],['8','10','12','14','16','18','20'],fontsize=11,fontweight=1.8)
	pl.yticks([-0.8,-0.4,0,0.4,0.8],['','','','',''],fontsize=11,fontweight=1.8)
	#-6_-16
	
	thea=[-16,-26]
	rg,zw1,zw2,zw3,zw5,zw6,w=cal_warpc((thea[0]+thea[1])/2.)
	rho_16,rherr_16,rcen_16,zcen_16,zrms_16,ind16 = [],[],[],[],[],[]
	pl.subplot(3,5,15)
	for i in range(14):
		zcen,zrms,vcen,vrms,rcen,rrms,radius,zerr,verr,rerr,zcerr,vcerr,rcerr,xcen,ycen,ind=cal_zcen_zrms(az,zz,vv,rr,mass,8.+i,8.+(i+1),thea[0],thea[1])
		if (len(np.transpose(ind))>0):
			pl.scatter(rr[ind],zz[ind],s=all_dist[ind],alpha=0.5)
		if np.isnan(zrms) == False:
			rho = 18.525*vrms**2/(zrms)**2
			rherr = np.sqrt(4.*(verr/vrms)**2+4.*(zerr/zrms)**2)
			rho_16.append(rho)
			rherr_16.append(rherr)
			rcen_16.append(rcen)
			zcen_16.append(zcen)
			zrms_16.append(zrms)
			pl.errorbar(rcen,zcen,yerr=zrms*2.355/2.,fmt='^',c=col_err,markersize=5.,elinewidth = 1, capsize=2)
	pl.text(8.5,0.8,'[%s,%s]' %(thea[0],thea[1]),fontsize=11,fontweight=1.8)
	pl.xlim([8,21.5])
	pl.ylim([-0.9,0.9])
#	pl.plot(rg,w,'--',color='black')
#	pl.plot(rg,zw1,'-',color='blue')
#	pl.plot(rg,zw2,'--',color='red')
#	pl.plot(rg,zw5,'--',color='green')
	pl.plot([8,21.5],[0,0],'--',color='red')
	pl.xticks([8,10,12,14,16,18,20],['8','10','12','14','16','18','20'],fontsize=11,fontweight=1.8)
	pl.yticks([-0.8,-0.4,0,0.4,0.8],['','','','',''],fontsize=11,fontweight=1.8)
    
	fig1.savefig('/home/yansun/HDwmisp/survey/reduced/r_z_corrugation_residual.png',format='png',bbox_inches='tight')
	fig2, ax = pl.subplots(nrows=1,ncols=1,sharex=True, figsize=[6,6])
	pl.plot(rcen_1[:-1],2.355*np.array(zrms_1[:-1]),'o')
	pl.plot(rcen_2[:-1],2.355*np.array(zrms_2[:-1]),'o')
	pl.plot(rcen_3[:-1],2.355*np.array(zrms_3[:-1]),'o')
	pl.plot(rcen_4[:-1],2.355*np.array(zrms_4[:-1]),'o')
	pl.plot(rcen_5[1:-1],2.355*np.array(zrms_5[1:-1]),'o')
	pl.plot(rcen_6[1:-1],2.355*np.array(zrms_6[1:-1]),'o')
	pl.plot(rcen_7[1:-1],2.355*np.array(zrms_7[1:-1]),'o')
	pl.plot(rcen_8[:-3],2.355*np.array(zrms_8[:-3]),'o')
	pl.plot(rcen_9[:-1],2.355*np.array(zrms_9[:-1]),'o')
	pl.plot(rcen_10[:-1],2.355*np.array(zrms_10[:-1]),'o')
	pl.plot(rcen_11[:-2],2.355*np.array(zrms_11[:-2]),'o')
	pl.plot(rcen_12[:-3],2.355*np.array(zrms_12[:-3]),'o')
	pl.plot(rcen_13[:-3],2.355*np.array(zrms_13[:-3]),'o')
	pl.plot(rcen_14[:-3],2.355*np.array(zrms_14[:-3]),'o')
	pl.plot(rcen_15[:-3],2.355*np.array(zrms_15[:-3]),'o')
#	pl.plot(rcen_16[:-3],2.355*np.array(zrms_16[:-3]),'o')
	rg1 = np.linspace(8,20,100)
	br1 = cal_flaring(rg1)
	pl.plot(rg1,br1,'--')
pl.xlim([8,21])
	
pl.show()
'''