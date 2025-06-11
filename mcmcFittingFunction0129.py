import numpy as np 
import matplotlib.pyplot as plt
import emcee
import corner
import os

np.random.seed(42)

#comp1/comp2 * dR * sin
component = 2
excluded = True
sin = True 	#set to True to fit sin component and plot 3D residual
print('Running with %i component(s) and %s l in [195, 200] and %s SIN component' % (component, 'excluding' if excluded else 'including', 'with' if sin else 'without'))

if component == 1: path = 'oneComp'
elif component ==2: path = 'twoComp'
else: pass
if excluded: path += 'Exc'

#cat = 'residual.razm'
#cat = 'cloud.razm'
#cat = 'catalog/osc_sy.txt'
#cat = 'cloud_notcorrect.razm'
'''
Set parameters
parameter should be
(inital_value, 'free'/'fixed', value_range)
'''
params = {}
params[r'$a_0$'] = (-0.0, 'fixed', [-0.1, +0.1])	#a_0 in kpc [ 0.09414761  8.57497436  1.04215655 -1.06578696]

###first warp compo
params[r'$a_1$'] = (0.0, 'free' if not sin else 'fixed', [-0.501, 0.50])	#a_1 in kpc
params[r'$R_{w1}$'] = (8.56849936, 'free' if not sin else 'fixed', [7, 11])		#R_w1 in kpc
params[r'$b_{w1}$'] = (1, 'free' if (component==1) & (not sin) else 'fixed', [0.1, 2.5])	#b_w1 index 0.9435
params[r'$\phi_{w1}$'] = (-0.7049789, 'free' if not sin else 'fixed', [-90, 90])	#phi_w1 in deg

###second warp component
params[r'$a_2$'] = (0.0, 'free' if (component==2) & (not sin) else 'fixed', [-0.50, 0.50])	#a_2 in kpc  -0.21
params[r'$R_{w2}$'] = (12.7166892, 'free' if (component==2) & (not sin) else 'fixed', [10, 17])		#R_w2 in kpc
params[r'$b_{w2}$'] = (2.03007431, 'free' if (component==2) & (not sin) else 'fixed', [0.1, 2.5])	#b_w2 index
params[r'$\phi_{w2}$'] = (-20.4479290, 'free' if (component==2) & (not sin) else 'fixed', [-90, 90])	#phi_w2 in deg

###sin component
params[r'$a_3$'] = (0.0, 'free' if sin else 'fixed', [0, 0.50]) #a_3 in kpc
params[r'$R_{c}$'] = (7.54, 'free' if sin else 'fixed', [6, 9])  #start of sin component in kpc
params[r'$P_0$'] = (1.9, 'free' if sin else 'fixed', [1.1, 5.0])  #period of sin component in kpc
params[r'$P_1$'] = (0.0, 'free' if sin else 'fixed', [-0.2, 0.3])  #period increasement of sin component in kpc
params[r'$\phi_{sin}$'] = (33.5, 'fixed', [28, 50])  #phi center in deg
params[r'$\phi_{sin,width}$'] = (4,'fixed', [1, 5]) #phi width in deg
params[r'$a_4$'] = (0.00, 'free' if sin else 'fixed', [-0.2, 0.1])  #baseline

#[1.39573585e-01, 7.74237334, 2.86256551, 1.57132833e-01,3.21995180e+01, 9.11163275, 1.06050455e-02]
#params[r'$a_5$'] = (0.0102, 'free', [0.001, 6])		#period of sin component in kpc
## [-0.07616907  0.10783739  7.79992154  0.9348509  -9.02005785]   constant xf
##sin [1.36492974e-01 8.47536353e+00 1.22706616e+00 4.68903037e-01 2.94463588e+01 9.84374395e+00 7.01067072e-04]
## [-0.069754    0.13771357  8.12746402  0.78466442 -8.47120866]   ##xfactor
###MCMC setting  [ -0.08204349   0.10188      7.77984319 -10.99553919  -0.11923537 13.79306634 -23.78646842]

nwalkers = 32 	#how many thread
burnin = 1000 	#number of iteration to reach local minimum (usually 20~30% of niter)
niter = 3000	#number of iteration

'''
one component: fix a0=0, free (a1, Rw1, bw1, phiw1)
two component: fix a0=0, bw1=1, free (a1, Rw1, phiw1, a2, Rw2, bw2, phiw2)
phiw at diff R: free (phiW1) only
sin component: fix phisin=33.5, phisinwidth=4, free (a3, Rsin, P0, P1, a4)
'''


###get free parameters
initial = []
free_params_name = []
free_params_limit = []
ndim = 0
for k, v in params.items():
	if v[1] == 'free':
		initial.append(v[0])
		free_params_name.append(k)
		free_params_limit.append(v[2])
		ndim += 1

best = initial


def fp2p(free_params):
	#free params -> mixed params with free and fixed
	'''
	return params[r'$a_0$'][0], params[r'$a_1$'][0], params[r'$R_{w0}$'][0], params[r'$b_{w0}$'], params[r'$\phi_{w0}$'], \
		free_params
	'''
	i = 0
	p = []
	for v in params.values():
		if v[1]=='free':
			p.append(free_params[i])
			i += 1
		else:
			p.append(v[0])
	return p


def function(x, free_params, warp=not sin, sin=sin):
	a0, a1, Rw1, bw1, PHIw1, a2, Rw2, bw2, PHIw2, a3, Rsin, Period0, Period1, PHIsin, PHIsinwid, a4 = fp2p(free_params)
	#a0, a1, Rw1, bw1, PHIw1, a2, bw2, PHIw2, a3, Rsin, Period0, Period1, PHIsin, PHIsinwid, a4 = fp2p(free_params)
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

def lnlike(free_params, data):
	#likelihood, or probability
	R, PHI, Z, mass = data
	LnProb = -0.5 * np.sum( (Z - function((R, PHI), free_params))**2 * mass )
	return LnProb
def lnlike_normalized(free_params, data):
	#likelihood, or probability
	R, PHI, Z, mass = data
	LnProb = -0.5 * np.sum( (Z - function((R, PHI), free_params))**2 * norm_mass )
	return LnProb
def lnprior(free_params):
	#criteria of fitting
	for p, lim in zip(free_params, free_params_limit):
		if not (lim[0] < p < lim[1]): return -np.inf
	return 0.0

def lnprob(free_params, data):
	lp = lnprior(free_params) #call lnprior
	if not np.isfinite(lp): #check if lp is infinite:
		return -np.inf
	return lp + lnlike(free_params, data) #recall if lp not -inf, its 0, so this just returns likelihood


def main(*data, initial=initial, nwalkers=nwalkers, burnin=burnin, niter=niter, ndim=ndim, lnprob=lnprob):
	#guess initial values
	ndim = len(initial)
	p0 = [np.array(initial) + 1e-7 * np.random.randn(ndim) for i in range(nwalkers)]

	#mcmc sampler
	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=data, \
		moves=emcee.moves.StretchMove(0.5))
	
	print("Initial = ", initial)
	print("Running burn-in...")
	p0, _, _ = sampler.run_mcmc(p0, burnin)
	print("After burn-in:\n", p0)
	sampler.reset()

	print("Running production...")
	pos, prob, state = sampler.run_mcmc(p0, niter)
	
	return sampler, pos, prob, state


def _excludedXYZ(x, y, z):
	### remove excluded surface
	r, phi = XY2RPHI(x, y)
	from warp_out_b import gc2g
	l, b, d = gc2g(r, phi, z)
	idx = (l>165) & (l<199) & (x**2+y**2>16**2)
	x[idx] = np.nan
	y[idx] = np.nan
	z[idx] = np.nan
	return x, y, z


def convolveZ(x, y, z, w, kernelsize, sampleX, sampleY):
	# get average z at (sampleX, sampleY) by convolving points (x,y,z,w) along z direction
	distance2 = ((sampleX[...,np.newaxis] - x)**2 + (sampleY[...,np.newaxis] - y)**2)
	sigma2 = kernelsize**2 / (8*np.log(2))
	kernelweight = np.exp(- distance2 / 2 / sigma2)

	kernelweight[distance2>kernelsize**2*4] = 0
	#w=1
	sampleZ = np.sum(z * w * kernelweight, axis=-1) / np.sum(w * kernelweight, axis=-1)

	### remove excluded l range
	sampleR, samplePHI = XY2RPHI(sampleX, sampleY)
	from warp_out_b import gc2g
	l, b, d = gc2g(sampleR, samplePHI, sampleZ)

	idx = ( (l>165) & (l<199) & (sampleR>16) ) | ( (l>165) & (l<195) )
	sampleZ[idx] = np.nan

	return sampleZ


def XY2RPHI(x, y):
	r = np.sqrt(x**2 + y**2)
	phi = np.arctan2(x, y)/np.pi*180
	return r, phi

def RPHI2XY(r, phi):
	x = r*np.sin(phi/180*np.pi)
	y = r*np.cos(phi/180*np.pi)
	return x, y

def aicbic(data, free_params, lnlike=lnlike_normalized):
	k = len(free_params)
	lnlk = lnlike(free_params, data)
	n = len(data[0])
	aic = 2 * k - 2 * lnlk
	bic = k * np.log(n) - 2 * lnlk
	print('AIC=%f, BIC=%f' % (aic, bic))

if __name__ == '__main__':
	if 1:
		###load observation
		l,b,v,lwidth,distance,err_dist,size_pc,mas,\
		complete,surface_density,viral_para,gal_ridus,\
		xx,yy,zz,index,l_rms,b_rms,v_rms,beta,T_int,area,tpeak,tex,nh2,pint,fill,lp,bp,vp = \
		np.hstack((np.loadtxt('osc2_para.txt').T, np.loadtxt('out_para.txt').T, np.loadtxt('per_para.txt').T))
		#ind = ~((l >= 15) & (l < 20) | (l > 160) & (l <= 165) | (l >= 195) & (l < 200))
		#ind = (mas > 0) & ~((l > 160) & (l <= 165) | (l >= 195) & (l < 200)) & (gal_ridus >=16) & (gal_ridus < 19)

		'''
		import pandas as pd
		df = pd.read_csv('clumps.csv')
		xx = df['xx'].to_numpy()
		yy = df['yy'].to_numpy()
		zz = df['zz'].to_numpy()
		l = df['l'].to_numpy()
		b = df['b'].to_numpy()
		mas = df['mass'].to_numpy()
		'''

		if excluded:
			#ind = (mas > 0) & ~((l > 160) & (l <= 165) | (l >= 195) & (l < 200))
			#ind = (mas > 0) & ((l <= 160) | (l >= 200))
			### exclude distant but massive clouds within 195~200
			ind = (mas > 0) & ~((l>195) & (l<199) & (xx**2+yy**2>16**2))# & (mas>1e4))
		else:
			ind = (mas > 0)
		print('Total:', ind.size, 'After excluding:', ind.sum())
		X = xx[ind]
		Y = yy[ind]
		Z = zz[ind]
		mass = mas[ind] #/complete[ind]) #*gal_ridus[ind]
		norm_mass = mass/np.sum(mass)

		R, PHI = XY2RPHI(X, Y)

		#idx = (l>180) & (l<200)
		#plt.scatter(X, Y, c=Z, s=mass/1e3, cmap='RdYlBu_r', vmin=-0.5, vmax=0.5)
		#plt.show()

		#data = data.T[mass>1e3].T
		'''
		data = np.loadtxt(cat).T
		R, PHI, Z, mass = data

		X, Y = RPHI2XY(R, PHI)
		'''

		### filter R range
		if 0:
			i = 8
			R_sep = [9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
			R1, R2 = R_sep[i:i+2]
			print(R1, R2)
			idx = (R>=R1) & (R<R2)
			R = R[idx]
			PHI = PHI[idx]
			Z = Z[idx]
			mass = mass[idx]
			norm_mass = mass/np.sum(mass)
			suffix = '_%i_%i' % (R1, R2)
		elif sin:
			### sin component
			suffix = '_sin'
		else:
			suffix = ''


		data = np.vstack([R, PHI, Z, mass])

		### load residual data
		if sin:
			data = np.load('residual_%icomp.npy' % (component))
			R, PHI, Z, mass = data
			X, Y = RPHI2XY(R, PHI)

			#plt.scatter(X, Y, c=Z, cmap='RdYlBu_r', vmin=-0.3, vmax=0.3)
			#plt.show()

		print('Data in shape:', data.shape)
	else:
		###simulate
		ns = 1000
		R = np.random.normal(9, 10, ns)
		R[R<0] = -R[R<0]
		PHI = np.random.rand(ns) * 150 + 30
		mass = np.random.rand(ns) + 0.1

		###a0, a, Rw, PHIw, b
		Z = function((R, PHI), initial)
		Z += np.random.normal(0, 5, ns)

		print(np.isnan(Z).sum())

		data = np.vstack([R, PHI, Z, mass])


	if 0:
		###run MCMC
		sampler, pos, prob, state = main(data)
		#samples = sampler.flatchain
		steps = sampler.flatchain
		probs = sampler.flatlnprobability
		np.save(os.path.join(path, 'steps%s.npy' % suffix), steps)
		np.save(os.path.join(path, 'probs%s.npy' % suffix), probs)
	else:
		steps = np.load(os.path.join(path, 'steps%s.npy' % suffix))
		probs = np.load(os.path.join(path, 'probs%s.npy' % suffix))

	### get median / best / best med
	if 1:
		med  = np.median(steps, axis=0)
		best = steps[np.argmax(probs)]
		bestmed = np.mean(steps[probs > np.percentile(probs, 95)], axis=0)
		print('The best fitting result is: ', best)
		print('The median fitting result is: ', med)
		print('The median of the best is: ', bestmed)



	### AIC / BIC
	if 1:
		# Number of data points and free parameters
		n = len(data[0])  # number of data points
		k = len(best)     # number of free parameters

		# Log-likelihood for the best-fitting parameters
		log_likelihood = lnlike_normalized(best, data)

		# Calculate AIC and BIC
		aic = 2 * k - 2 * log_likelihood
		bic = k * np.log(n) - 2 * log_likelihood


		print('AIC for the best-fitting model:', aic)
		print('BIC for the best-fitting model:', bic)
		aicbic(data, bestmed)


	### corner plot
	if 0:
		def corners(steps, probs, bins=10, labels=[None]*ndim, point=best, top=90):
			def lim(s):
				ma = np.nanmax(s)
				mi = np.nanmin(s)
				return mi-(ma-mi)*0.1, ma+(ma-mi)*0.1
			steps = steps[np.argsort(probs)]
			probs = np.sort(probs)
			def scientificNotationTitle(v, l, u):
				e = np.floor(np.log10(np.abs(v)))
				if e == 1: return '$%+.2f_{%+.2f}^{%+.2f}$' % (v, l, u)
				elif e == 0: return '$%+.3f_{%+.3f}^{%+.3f}$' % (v, l, u)
				else: return '$%+.3f_{%+.3f}^{%+.3f}$ x$10^{%i}$' % (v/10**e, l/10**e, u/10**e, e)
			fig, ax = plt.subplots(nrows=ndim, ncols=ndim, figsize=(9,9))
			plt.subplots_adjust(bottom=0.09, top=0.91, left=0.09, right=0.91, wspace=0.05, hspace=0.05)
			### only plot the top% steps
			idx = probs > np.percentile(probs, top)
			for i in range(ndim):
				### hist at [i,i]
				h,xe,_ = ax[i,i].hist(steps[:,i][idx], color='grey', bins=bins, density=True)
				### title
				lower = np.quantile(steps[:,i], 0.16)-point[i]
				upper = np.quantile(steps[:,i], 0.84)-point[i]
				ax[i,i].text(0, 1.02, '%s=%s' % (labels[i], scientificNotationTitle(point[i], lower, upper)), \
					ha='left', va='bottom', fontsize=10, transform=ax[i,i].transAxes)
				### marker
				peak = float(h[np.argwhere(point[i]>xe)[-1]])
				ax[i,i].plot([point[i],], peak*1.1, 'v', color='#e33a23')
				ax[i,i].plot([point[i]+lower]*2, [0, peak], '--', color='k', lw=1, alpha=0.5)
				ax[i,i].plot([point[i]+upper]*2, [0, peak], '--', color='k', lw=1, alpha=0.5)

				### hide some xticks and all yticks
				if i<ndim-1: ax[i,i].xaxis.set_ticklabels([])
				ax[i,i].set_yticklabels([])

				### ticks
				ax[i,i].minorticks_on()
				ax[i,i].tick_params(top=True, left=False, direction='in', labelrotation=45, length=5)
				ax[i,i].tick_params(which='minor', top=True, left=False, direction='in', length=2)
				#ax[i,i].tick_params(labelrotation=45)
				ax[i,i].set_xlim(lim(steps[:,i][idx]))

				### plot hist2d
				for j in range(i):
					### hist
					h,xe,ye,_ = ax[i,j].hist2d(steps[:,j][idx], steps[:,i][idx], bins=bins, cmap='gray_r')
					#ax[i,j].scatter(steps[:,j][idx], steps[:,i][idx], c=probs[idx], s=0.2, cmap='RdYlBu_r')
					#x,y = (xe[1:]+xe[:-1])/2, (ye[1:]+ye[:-1])/2
					#ax[i,j].contour(x, y, h.T, colors='gray',levels=[h.max()*0.5, h.max()*0.75])
					#ax[i,j].scatter(steps[:,j], steps[:,i], c=probs, s=0.2)#, vmin=np.percentile(probs,0.8))
					#ax[i,j].plot(best[j],best[i], 'gx')
					### marker
					ax[i,j].plot(point[j],point[i], '+', color='#e33a23', ms=8, mew=2)
					
					### hide some
					if i<ndim-1: ax[i,j].xaxis.set_ticklabels([])
					if j>0: ax[i,j].yaxis.set_ticklabels([])

					### ticks
					ax[i,j].minorticks_on()
					ax[i,j].set_yticks(ax[i,i].get_xticks())
					ax[i,j].tick_params(top=True, right=True, direction='in', labelrotation=45, length=5)
					ax[i,j].tick_params(which='minor', top=True, right=True, direction='in', length=2)
					ax[i,j].set_xlim(lim(steps[:,j][idx]))
					ax[i,j].set_ylim(lim(steps[:,i][idx]))
					### hide upper right
					ax[j,i].axes.set_axis_off()

			if sin:
				if component==1:
					ax[0,0].text(-0.25, 1, 'a', ha='left', va='top', color='black', font=dict(size=18, family="Arial Black"), transform=ax[0,0].transAxes)
				else:
					ax[0,0].text(-0.25, 1, 'b', ha='left', va='top', color='black', font=dict(size=18, family="Arial Black"), transform=ax[0,0].transAxes)

		#plt.plot(probs, '.')
		#plt.show()
		corners(steps[2000:], probs[2000:], bins=31, labels=free_params_name, point=bestmed, top=40)
		plt.savefig('fig/corner_%icomp%s.png' % (component, '_sin' if sin else ''), format='png', bbox_inches='tight', dpi=400)
		#corner.corner(steps, show_titles=True, plot_contours=False, plot_datapoints=True, quantiles=[0.16, 0.5, 0.84], range=None)#[(v-0.5,v+0.5) for v in best])
		plt.show()


	### export residual
	if not sin:
		data[2] -= function(data[:2], bestmed, warp=True, sin=False)
		print(data[:2], function(data[:2], bestmed, warp=True, sin=False))
		np.save('residual_%icomp.npy' % component, data)
		print('Export to residual_%icomp.npy' % component)


	### plot face-on gplane
	if 0:
		def gplane(fig, ax, x, y, c, s, im=False, **kws):
			if im:
				h, xe, ye = np.histogram2d(x, y, bins=[np.arange(-10,18,0.6), np.arange(-15,24,0.6)], weights=c)
				n, xe, ye = np.histogram2d(x, y, bins=[np.arange(-10,18,0.6), np.arange(-15,24,0.6)])
				im = ax.imshow(h.T/n.T, origin='lower', extent=(xe[0], xe[-1], ye[0], ye[-1]), **kws)
			else:
				im = ax.scatter(x, y, c=c, s=s, **kws)
			ax.set_aspect('equal')
			ax.plot(0,0,'ko')
			xlim = ax.get_xlim()
			ylim = ax.get_ylim()
			for ray in range(0, 360, 10):
				ax.plot([0,np.sin(ray/180*np.pi)*30],[0,np.cos(ray/180*np.pi)*30], '--', color='gray', linewidth=1)
			for rad in range(2, 30, 2):
				ax.plot(rad*np.sin(np.linspace(0, 2*np.pi, 360)), rad*np.cos(np.linspace(0, 2*np.pi, 360)), '--', color='gray', linewidth=1)
			fig.colorbar(im, ax=ax)
			ax.set_xlim(xlim)
			ax.set_ylim(ylim)

		fig, ax = plt.subplots(ncols=2, nrows=2, sharex=True, sharey=True, figsize=(14, 9))
		kws = dict(cmap='RdBu', vmin=-0.10, vmax=0.10)
		gplane(fig, ax[0,0], X, Y, c=Z, s=mass/mass.max()*200, **kws)
		gplane(fig, ax[0,1], X, Y, c=Z-function((R, PHI), best), s=mass/mass.max()*200, **kws)

		kws = dict(cmap='RdBu', vmin=-0.010, vmax=0.010)
		gplane(fig, ax[1,0], X, Y, c=Z, s=mass/mass.max()*200, im=True, **kws)
		gplane(fig, ax[1,1], X, Y, c=Z-function((R, PHI), best), s=mass/mass.max()*200, im=True, **kws)



	###prepare functions for 3D visualization
	if True:
		def _point(x, y, z):
			xyz = np.vstack((x,y,z)).T
			point = pv.PointSet(xyz)
			return point

		def _surface(gx, gy, gz, c=None):
			surface = pv.StructuredGrid(gx, gy, gz)
			return surface

		def _hist2dData2Grid(step=0.5):
			### 
			bins = [np.arange(-10, 18, step), np.arange(-15, 24, step)]
			gridZM, xe, ye = np.histogram2d(X, Y, bins=bins, weights=Z*weight)
			gridM, xe, ye = np.histogram2d(X, Y, bins=bins, weights=weight)
			gridX, gridY = np.meshgrid((xe[1:]+xe[:-1])/2, (ye[1:]+ye[:-1])/2)
			gridM[gridM==0] = np.nan
			gridZ = (gridZM / gridM).T
			return gridX, gridY, gridZ

		def _convolveData2Grid(step=None, polar=True, residual=False, **kws):
			'''
			step: grid step for (dx, dy) or (dr, dphi)
			polar: grid in xy or polar
			residual: return data or residual
			'''
			if polar:
				if step is None: step=(.25, 2)	#
				sampleR, samplePHI = np.meshgrid(np.arange(1, 22.1, step[0]), np.arange(0, 360.1, step[1]))
				sampleX, sampleY = RPHI2XY(sampleR, samplePHI)
			else:
				if step is None: step=(.25, .25)
				sampleX, sampleY = np.meshgrid(np.arange(-10, 18, step[0]), np.arange(-15, 24, step[1]))
			
			if not residual:
				sampleZ = convolveZ(X, Y, Z, mass, 0.75, sampleX, sampleY)
			else:
				sampleZ = convolveZ(X, Y, Z-function((R, PHI), bestmed, **kws), mass, 0.75, sampleX, sampleY)

			return sampleX, sampleY, sampleZ

		def _modelGrid(step=(1,15), **kws):
			gridR, gridPHI = np.meshgrid(np.arange(8, 22.1, step[0]), np.arange(0, 360.1, step[1]))
			gridZ = function((gridR, gridPHI), bestmed, **kws)
			gridX, gridY = RPHI2XY(gridR, gridPHI)
			return gridX, gridY, gridZ

		def _smoothModelGrid(plotter, step=(1,15), **kws):
			for phi in np.arange(0, 360, step[1]):
				lineR = np.arange(8, 22.1, 0.2)
				linePHI = np.repeat(phi, lineR.size)
				lineX, lineY = RPHI2XY(lineR, linePHI)
				lineZ = function((lineR, linePHI), bestmed, **kws)
				xyz = np.vstack((lineX, lineY, lineZ))
				plotter.add_lines(xyz, color='k', width=2, connected=True)

			for r in np.arange(8, 23.1, 1):
				linePHI = np.linspace(0, 360, int(360*r)//40)
				lineR = np.repeat(r, linePHI.size)
				lineX, lineY = RPHI2XY(lineR, linePHI)
				lineZ = function((lineR, linePHI), bestmed, **kws)
				xyz = np.vstack((lineX, lineY, lineZ)).T
				plotter.add_lines(xyz, color='k', width=2, connected=True)
	

		def _residualGrid(**kws):
			gridX, gridY, gridZ = _dataConvolve(**kws)
			gridR, gridPHI = XY2RPHI(gridX, gridY)
			gridResidual = gridZ - function((gridR, gridPHI), bestmed, sin=False)
			return gridX, gridY, gridResidual


	### 3d visualization with pyvista
	if 0:
		import pyvista as pv

		### plot axes
		p = pv.Plotter(window_size=(1600,1600))#, lighting=None)

		if 1:
			### data + model

			### DATA points
			#p.add_points(_point(X, Y, Z), style='points_gaussian', scalars=np.log10(mass), cmap='gray', clim=[0,5], opacity=0.2, point_size=5, show_scalar_bar=False)
			#p.add_points(_point(X, Y, Z), scalars=Z, style='points_gaussian', render_points_as_spheres=False, emissive=False, cmap='RdYlBu_r', clim=[-1,1], opacity=0.3, point_size=5, show_scalar_bar=False)

			### DATA regrid surface
			### use a small step for surface to make it smooth, eg. (0.25, 2)
			### use a larger step for wireframe to make it clear to see, eg. (0.5, 4)
			#gridX, gridY, gridZ = _convolveData2Grid(step=(0.75, 5), polar=True, warp=False, residual=False, sin=False)
			#gridX, gridY, gridZ = _dataConvolve(step=(0.5, 4), polar=True)
			#p.add_mesh(_surface(gridX, gridY, gridZ), style='points_gaussian', color='k', clim=[-1,1], opacity=1.0, point_size=14, render_points_as_spheres=False)
			
			#p.add_points(_point(gridX.ravel(), gridY.ravel(), gridZ.ravel()), style='points_gaussian', scalars=gridZ.ravel(), cmap='RdYlBu_r', clim=[-1,1], opacity=1.0, point_size=3, render_points_as_spheres=False)
			#p.add_mesh(_surface(gridX, gridY, gridZ), style='surface', scalars=gridZ.T, cmap='RdYlBu_r', clim=[-1,1], opacity=0.9, point_size=5, show_scalar_bar=False)

			### MODEL wireframe
			#gridX, gridY, gridZ = _modelGrid(sin=False)
			grid = _surface(gridX, gridY, gridZ)
			grid.texture_map_to_plane(inplace=True)

			p.add_mesh(grid, style='surface', opacity=1.0, line_width=1, texture = pv.read_texture('ssc2008-10a_ext.jpg'))
			#p.add_mesh(_surface(gridX, gridY, gridZ), style='surface', scalars=gridZ.T, cmap='RdYlBu_r', clim=[-1,1], point_size=10)

			### smoooth MODEL wireframe
			#_smoothModelGrid(p, sin=False)

			p.set_scale(xscale=1, yscale=1, zscale=3)	#scale z by 10
			p.set_position([-120, 70, 12])	#position of eye
			p.set_focus([0, 0, 0])	#position to look at
			p.set_viewup([0, 0, 1])	#z is up

			p.set_background('white')
			p.show_bounds(bounds=[-25,25,-25,25,-2,2], grid=True, bold=True, \
				font_size=10, color='gray', location='outer', padding=0, \
				xtitle='X (kpc)', \
				ytitle='Y (kpc)', \
				ztitle='Z (kpc)')

		else:
			### residual

			#for i, phi in enumerate([155.0, 149.75, 144.5, 139.25, 134.0, 128.75, 123.5, 118.25, 113.0, 107.75, 102.5, 97.25, 92.0, 86.75, 81.5, 76.25, 71.0, 65.75, 60.5, 55.25, 50.0, 44.75, 39.5, 34.25, 29.0, 23.75, 18.5, 13.25, 8.0, 0.0, -6.0, -11.0, -16.0, -21.0, -26.0]):
			for i, phi in enumerate([149.75, 139.25, 128.75, 118.25, 107.75, 97.25, 86.75, 76.25, 65.75, 55.25, 44.75, 34.25, 23.75, 13.25, -11.0, -21.0]):
				lineR = np.arange(8, 22.1, 0.2)
				linePHI = np.repeat(phi, lineR.size)
				lineX, lineY = RPHI2XY(lineR, linePHI)
				lineResidual = convolveZ(X, Y, Z, mass, 0.5, lineX, lineY) - function((lineR, linePHI), best, sin=False)
				xyr = np.vstack((lineX, lineY, lineResidual)).T
				p.add_lines(xyr, color='k', width=2)
				#p.add_mesh(pv.StructuredGrid(lineX, lineY, lineResidual), style='wireframe', scalars=lineResidual, color='RdYlBu_r', clim=[-.05,.05], opacity=1.0, line_width=5 if i%2==1 else 0)

			for r in np.arange(8, 22.1, 1):
				linePHI = np.arange(-30, 180, 1)
				lineR = np.repeat(r, linePHI.size)
				lineX, lineY = RPHI2XY(lineR, linePHI)
				lineResidual = convolveZ(X, Y, Z, mass, 0.5, lineX, lineY) - function((lineR, linePHI), best, sin=False)
				xyr = np.vstack((lineX, lineY, lineResidual)).T
				p.add_lines(xyr, color='k', width=2)
				#p.add_mesh(pv.StructuredGrid(lineX, lineY, lineResidual), style='wireframe', scalars=lineResidual, color='RdYlBu_r', clim=[-.05,.05], opacity=1.0, line_width=5 if r%1==0 else 0)
			
			gridX, gridY, gridZ = _residualGrid(step=(0.2, 1), polar=True)
			#p.add_mesh(_surface(gridX, gridY, gridZ), style='wireframe', scalars=gridZ.T, cmap='RdYlBu_r', clim=[-.1,.1], opacity=1.0, line_width=3)
			p.add_mesh(_surface(gridX, gridY, gridZ), style='surface', scalars=gridZ.T, cmap='RdYlBu_r', clim=[-.1,.1], opacity=0.7)
		
			#gridX, gridY, gridZ = _modelGrid(warp=False)
			#p.add_mesh(_surface(gridX, gridY, gridZ), style='wireframe', color='k', opacity=1.0, line_width=2)
			#p.add_mesh(_surface(gridX, gridY, gridZ), style='surface', scalars=gridZ.T, cmap='RdYlBu_r', clim=[-1,1], point_size=10)

			p.set_scale(xscale=1, yscale=1, zscale=5)	#scale z by 10
			p.set_position([-60, 30, 15])	#position of eye
			p.set_focus([4, 3, 0])	#position to look at
			p.set_viewup([0, 0, 1])	#z is up

			p.set_background('white')
			p.show_bounds(bounds=[-10,20,-15,20,-1,1], grid=True, bold=True, \
				font_size=10, color='gray', location='outer', padding=0, \
				xlabel='X (kpc)', \
				ylabel='Y (kpc)', \
				zlabel='Z - model (kpc)')
		
		#p.show_bounds(all_edges=False)
		p.add_points(pv.PointSet([[-24,-24,-2],[-24,-24, 2]]))
		'''
		###xy axis
		gx, gy = np.meshgrid(np.arange(-24,26,8), np.arange(-24,26,8))
		gz = np.zeros(gx.shape)-2
		axis = pv.StructuredGrid(gx, gy, gz)
		p.add_mesh(axis, style='wireframe', color='lightgrey', line_width=2)

		###xz axis
		gx, gz = np.meshgrid(np.arange(-24,26,8), np.arange(-2,3,1))
		gy = np.zeros(gx.shape)-24
		axis = pv.StructuredGrid(gx, gy, gz)
		p.add_mesh(axis, style='wireframe', color='lightgrey', line_width=2)

		###yz axis
		gy, gz = np.meshgrid(np.arange(-24,26,8), np.arange(-2,3,1))
		gx = np.zeros(gx.shape)+24
		axis = pv.StructuredGrid(gx, gy, gz)
		p.add_mesh(axis, style='wireframe', color='lightgrey', line_width=2)
		'''

		p.show()
		#p.save_graphic('datamodel.eps')


	### 3D visualization with plotly
	if not sin:
		#best = initial
		import plotly.graph_objs as go
		### Data
		#points = go.Scatter3d(x=X, y=Y, z=Z - function((R,PHI), bestmed, sin=False), mode='markers', \
		#	marker=dict(size=3, color=Z, colorscale='RdYlBu_r', cmin=-1, cmax=1, opacity=0.2, colorbar=dict(thickness=20)))
		#points = go.Scatter3d(x=X, y=Y, z=Z, mode='markers', \
		#	marker=dict(size=10, color=np.log10(mass), colorscale='gray', opacity=0.3))#, colorbar=dict(thickness=20)

		
		### median
		gridX, gridY, gridZ = _convolveData2Grid(step=(0.75, 5), polar=True, warp=False, residual=False, sin=False)
		gridR, gridPHI = XY2RPHI(gridX, gridY)

		#gridZ = gridZ - function((gridR, gridPHI), bestmed, warp1=True, sin=False)
		#gridR, gridPHI, gridZ = _convolveData2Grid(step=(0.5, 0.5), polar=False)
		points = go.Scatter3d(x=gridX.ravel(), y=gridY.ravel(), z=gridZ.ravel(), mode='markers', \
			marker=dict(size=5, color=gridZ.ravel(), colorscale='RdYlBu_r', cmin=-0.5, cmax=0.5, opacity=1.0, \
			colorbar=dict(thickness=20, title='Z (kpc)', x=0.24, y=0.12, len=0.3, orientation='h')))


		### wireframe
		wires = []
		line_marker = dict(color='#000000', width=1)#, colorscale='RdYlBu_r', cmin=-0.5, cmax=0.5)
		### R axis
		for phi in list(range(323, 155, -12))+[155, 143, 131, 119, 107, 95, 83, 71, 59, 47, 38, 29, 20, 11, -11, -17, -27]:
			#for phi in np.arange(0, 360, 15):
			lineR = np.arange(8.5, 20.6, 0.5)
			linePHI = np.repeat(phi, lineR.size)
			lineX, lineY = RPHI2XY(lineR, linePHI)
			lineZ = function((lineR, linePHI), bestmed, warp=True, sin=True)
			#line_marker['color'] = lineZ
			wires.append(go.Scatter3d(x=lineX, y=lineY, z=lineZ, mode='lines', line=line_marker))

		### phi axis
		for r in [8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 18.5, 20.5]:
			linePHI = np.linspace(0, 360, int(360*r)//40)
			lineR = np.repeat(r, linePHI.size)
			lineX, lineY = RPHI2XY(lineR, linePHI)
			lineZ = function((lineR, linePHI), bestmed, warp=True, sin=True)
			#line_marker['color'] = lineZ
			wires.append(go.Scatter3d(x=lineX, y=lineY, z=lineZ, mode='lines', line=line_marker))

		#sectorPosPhi = np.arange(165, -35, -15)-7.5
		sectorPosPhi = np.array([(p1+p2)/2 for p1,p2 in zip([155, 143, 131, 119, 107, 95, 83, 71, 59, 47, 38, 29, 20, -11, -17], [143, 131, 119, 107, 95, 83, 71, 59, 47, 38, 29, 20, 11, -17, -27])])
		sectorPosR = np.repeat(22.5, len(sectorPosPhi))
		sectorPosZ = function((sectorPosR, sectorPosPhi), bestmed)
		sectorPosX, sectorPosY = RPHI2XY(sectorPosR, sectorPosPhi)
		sectorText = []
		#for i in [0, 6, 12, 14]:
		for i in range(len(sectorPosPhi)):
			sectorText.append(dict(x=sectorPosX[i], y=sectorPosY[i], z=sectorPosZ[i], text='$\phi_{%i}$' % (i+1), showarrow=False, \
				xanchor = 'center', xshift=5, yanchor='middle', font=dict(color='#444444', size=20)))

		sectorPosR = np.array([(p1+p2)/2 for p1,p2 in zip([8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 18.5], [9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 18.5, 20.5])])
		sectorPosPhi = np.repeat(-42, len(sectorPosR))
		sectorPosZ = function((sectorPosR, sectorPosPhi), bestmed)
		sectorPosX, sectorPosY = RPHI2XY(sectorPosR, sectorPosPhi)
		for i in [0, 9]:#range(len(sectorPosR)):
			sectorText.append(dict(x=sectorPosX[i], y=sectorPosY[i], z=sectorPosZ[i], text='$R_{%i}$' % (i+1), showarrow=False, \
				xanchor = 'center', xshift=5, yanchor='middle', font=dict(color='#444444', size=15)))

		sun = go.Scatter3d(x=[0], y=[8.15], z=[0], mode='markers', \
				marker=dict(color='#ff6600', size=6))
		sunText = dict(x=0, y=8.15, z=0, text='Sun', showarrow=False,
			xanchor='left', xshift=5, yanchor='bottom', font=dict(color='#ff6600', size=16))

		gc = go.Scatter3d(x=[0], y=[0], z=[0], mode='markers', \
			marker=dict(color='#555555', size=6))
		gcText = dict(x=0, y=0, z=0, text='GC', showarrow=False,
			xanchor='left', xshift=5, yanchor='bottom', font=dict(color='#666666', size=16))


		eyeD = 1.6
		eyeAz = -35
		eyeEl = 25
		eyePos = dict(x = eyeD*np.cos(eyeEl*np.pi/180)*np.sin(eyeAz*np.pi/180),
			y = eyeD*np.cos(eyeEl*np.pi/180)*np.cos(eyeAz*np.pi/180),
			z = eyeD*np.sin(eyeEl*np.pi/180))

		
		### R vs phi_w
		if component==1:
			_r, _phi, _z = [], [], []
			for i in range(9,18):
				suffix = '_%i_%i' % (i, i+1)
				steps = np.load(os.path.join(path, 'steps%s.npy' % suffix))
				probs = np.load(os.path.join(path, 'probs%s.npy' % suffix))
				bestmed = np.mean(steps[probs > np.percentile(probs, 95)], axis=0)
				_r.append(i+0.5)
				_phi.append(bestmed[-1])
				_z.append(0)
			_x, _y = RPHI2XY(np.array(_r), np.array(_phi))
			_rphiw = go.Scatter3d(x=_x, y=_y, z=_z, mode='lines+markers',\
				marker=dict(symbol='square', size=2.5, color='darkgreen', opacity=1), line=dict(color='darkgreen', width=3))
			_r = [8.5, 20.5]
			_phi = [bestmed[-1]]*2
			_x, _y = RPHI2XY(np.array(_r), np.array(_phi))
			_phiw_1comp = go.Scatter3d(x=_x, y=_y, z=[0, 0], mode='lines',\
				marker=dict(symbol='square', size=1, color='darkgreen', opacity=1), line=dict(dash='longdash', color='darkgreen', width=3))
		

		### layout format
		layout = go.Layout(
			scene=dict(
				xaxis_title='X (kpc)',
				xaxis=dict(
					range=[-24.01, 24.01], tickvals=np.arange(-24,25,8),
					ticks='outside', tickwidth=3,
					gridcolor='#dddddd',
					#zerolinecolor='#000000',
					#zerolinewidth=3,
					showbackground=False,
					backgroundcolor='white',
					tickfont=dict(size=12)
				),
				xaxis_showspikes=True,
				yaxis_title='Y (kpc)',
				yaxis=dict(
					range=[-24.01, 24.01], tickvals=np.arange(-24,25,8),
					ticks='outside', tickwidth=3,
					gridcolor='#dddddd',
					#zerolinecolor='#000000',
					#zerolinewidth=3,
					showbackground=False,
					backgroundcolor='white',
					tickfont=dict(size=12)
				),
				zaxis_title='Z (kpc)',
				zaxis=dict(
					range=[-2, 2.21], tickvals=np.arange(-1, 3, 1),
					ticks='outside', tickwidth=3,
					gridcolor='#dddddd',
					zerolinecolor='#000000',
					zerolinewidth=3,
					showbackground=False,
					backgroundcolor='white',
					tickfont=dict(size=12)
				),
				annotations=[
					sunText,
					gcText,
					*sectorText,
					]
			),
			scene_camera = dict(
				up=dict(x=0, y=0, z=1),
				center=dict(x=0.11, y=0, z=-0.3),
				eye=eyePos
				#eye=dict(x=-0.4, y=-1.6, z=1.2)
				),
			showlegend=False,
			width=1200,
			height=900,
			margin=dict(r=0, l=0, b=0, t=0),
			scene_aspectmode='manual',
			scene_aspectratio=dict(x=1, y=1, z=.45),
		)


		# Create a Figure object and add the 3d objects to it
		fig = go.Figure(data=[points, *wires, sun, gc, _rphiw, _phiw_1comp], layout=layout)

		### panel ID
		fig.add_annotation(x=0.05, xref='paper', y=0.95, yref='paper', text='a', showarrow=False, font=dict(color='black', size=70, family="Arial Black"))

		### Show the plot
		if 0: fig.show()
		else: fig.write_image("fig/median_model_%icomp.png" % component, scale=3)


	### 3D visualization of residual with plotly
	if sin:
		import plotly.graph_objs as go
		
		
		### residual
		resX, resY, resZ = _convolveData2Grid(step=(0.25, 1), polar=True, residual=False)
		resX, resY, resZ = _excludedXYZ(resX, resY, resZ)
		residual = go.Surface(x=resX, y=resY, z=resZ, surfacecolor=resZ, colorscale='RdYlBu_r', cmin=-0.2, cmax=0.2, opacity=1,
			colorbar = dict(thickness=20, title='Residual (kpc)', x=0.7, y=0.05, len=0.3, orientation='h'))
		#residual = go.Scatter3d(x=resX.ravel(), y=resY.ravel(), z=resZ.ravel(), mode='markers', \
		#	marker=dict(size=5, color=resZ.ravel(), colorscale='RdYlBu_r', cmin=-0.2, cmax=0.2, opacity=1.0, \
		#	colorbar=dict(thickness=20, title='Z (kpc)', x=0.24, y=0.7, len=0.05, len=0.3, orientation='h')))
		

		### residual wireframe
		wires = []
		line_marker = dict(color='#223322', width=2)
		#for phi in np.arange(0, 360, 15):
		for phi in list(range(323, 155, -12))+[155, 143, 131, 119, 107, 95, 83, 71, 59, 47, 38, 29, 20, 11, -11, -17, -27]:
			lineR = np.arange(8, 22.1, 0.25)
			linePHI = np.repeat(phi, lineR.size)
			lineX, lineY = RPHI2XY(lineR, linePHI)
			lineZ = convolveZ(X, Y, Z, mass, 0.75, lineX, lineY)
			if phi in [29,38]: line_marker['width']=6
			else: line_marker['width']=2
			wires.append(go.Scatter3d(x=lineX, y=lineY, z=lineZ+0.01, mode='lines', line=line_marker))

		for r in np.arange(8, 22.1, 1):
			#linePHI = np.linspace(0, 360, int(360*r)//40)
			linePHI = np.arange(0, 360.1, 1)
			lineR = np.repeat(r, linePHI.size)
			lineX, lineY = RPHI2XY(lineR, linePHI)
			lineZ = convolveZ(X, Y, Z, mass, 0.75, lineX, lineY)
			wires.append(go.Scatter3d(x=lineX, y=lineY, z=lineZ+0.01, mode='lines', line=line_marker))

		
		#sectorPosPhi = np.arange(165, -35, -15)-7.5
		sectorPosPhi = np.array([(p1+p2)/2 for p1,p2 in zip([155, 143, 131, 119, 107, 95, 83, 71, 59, 47, 38, 29, 20, -11, -17], [143, 131, 119, 107, 95, 83, 71, 59, 47, 38, 29, 20, 11, -17, -27])])
		sectorPosR = np.linspace(20, 24, sectorPosPhi.size)
		sectorPosZ = np.zeros(sectorPosPhi.size)
		sectorPosX, sectorPosY = RPHI2XY(sectorPosR, sectorPosPhi)
		sectorText = []
		for i in range(len(sectorPosPhi)):
			#for i in [0, 10, 14]:#
			sectorText.append(dict(x=sectorPosX[i], y=sectorPosY[i], z=sectorPosZ[i], text='$\phi_{%i}$' % (i+1), showarrow=False, \
				xanchor = 'center', xshift=5, yanchor='middle', font=dict(color='#000000', size=30)))
		
		### bar ellipse
		pa = 30.5 /180*np.pi
		a = np.linspace(0,360,100) /180*np.pi
		x = 1.5*np.sin(a)
		y = 4.5*np.cos(a)
		z = np.zeros(x.shape)
		xrot = x * np.cos(pa) + y * np.sin(pa)
		yrot = - x * np.sin(pa) + y * np.cos(pa)

		bar = go.Mesh3d(x=xrot, y=yrot, z=z, color='gray', opacity=0.8, delaunayaxis='z')
		barText = dict(x=xrot.min(), y=yrot.min(), z=0, text='Bar', showarrow=False,
			xanchor='left', xshift=5, yanchor='bottom', font=dict(color='#666666', size=16))

		sun = go.Scatter3d(x=[0], y=[8.15], z=[0], mode='markers', \
				marker=dict(color='#ff6600', size=6))
		sunText = dict(x=0, y=8.15, z=0, text='Sun', showarrow=False,
			xanchor='left', xshift=5, yanchor='bottom', font=dict(color='#ff6600', size=16))

		gc = go.Scatter3d(x=[0], y=[0], z=[0], mode='markers', \
			marker=dict(color='#555555', size=6))
		gcText = dict(x=0, y=0, z=0, text='GC', showarrow=False,
			xanchor='left', xshift=5, yanchor='bottom', font=dict(color='#666666', size=16))

		### arrow
		arrowR = 25
		a = np.linspace(24.5, 53, 200) /180*np.pi
		arrowX = np.concatenate( (arrowR * np.sin(a), [(arrowR+0.5) * np.sin(a[-20])]) )
		arrowY = np.concatenate( (arrowR * np.cos(a), [(arrowR+0.5) * np.cos(a[-20])]) )
		arrowZ = np.zeros(arrowX.shape)
		arrow = go.Scatter3d(x=arrowX, y=arrowY, z=arrowZ, mode='lines', line=dict(color='black', width=4))

		eyeD = 1.5
		eyeAz = -67
		eyeEl = 45
		eyePos = dict(x = eyeD*np.cos(eyeEl*np.pi/180)*np.sin(eyeAz*np.pi/180),
			y = eyeD*np.cos(eyeEl*np.pi/180)*np.cos(eyeAz*np.pi/180),
			z = eyeD*np.sin(eyeEl*np.pi/180))

		### layout format
		layout = go.Layout(
			scene=dict(
				xaxis_title='X (kpc)',
				xaxis=dict(
					range=[-10.01, 25.01], tickvals=np.arange(-10,25.1,5),
					ticks='outside', tickwidth=3,
					gridcolor='#dddddd',
					#zerolinecolor='#000000',
					#zerolinewidth=3,
					showbackground=False,
					backgroundcolor='white',
					tickfont=dict(size=12)
				),
				xaxis_showspikes=True,
				yaxis_title='Y (kpc)',
				yaxis=dict(
					range=[-20.01, 25.01], tickvals=np.arange(-20,25.1,5),
					ticks='outside', tickwidth=3,
					gridcolor='#dddddd',
					#zerolinecolor='#000000',
					#zerolinewidth=3,
					showbackground=False,
					backgroundcolor='white',
					tickfont=dict(size=12)
				),
				zaxis_title='Residuals in Z (kpc)',#𝛿
				zaxis=dict(
					range=[-1, 1.01], tickvals=np.arange(-0.5, 2, 0.5),
					ticks='outside', tickwidth=3,
					gridcolor='#dddddd',
					zerolinecolor='#000000',
					zerolinewidth=3,
					showbackground=False,
					backgroundcolor='white',
					tickfont=dict(size=12)
				),
				annotations=[
					sunText,
					gcText,
					barText,
					*sectorText
					]
			),
			scene_camera = dict(
				up=dict(x=0, y=0, z=1),
				center=dict(x=0.10, y=-0.06, z=-0.32),
				eye=eyePos
				#eye=dict(x=-0.4, y=-1.6, z=1.2)
				),
			showlegend=False,
			width=1200,
			height=900,
			margin=dict(r=0, l=0, b=0, t=0),
			scene_aspectmode='manual',
			scene_aspectratio=dict(x=40/48, y=1, z=.3),
		)

		# Create a Figure object and add the 3d objects to it
		fig = go.Figure(data=[residual, *wires, bar, sun, gc, arrow], layout=layout)#

		### panel ID
		fig.add_annotation(x=0.05, xref='paper', y=0.95, yref='paper', text='a', showarrow=False, font=dict(color='black', size=60, family="Arial Black"))

		### Show the plot
		if 0: fig.show()
		else: fig.write_image("fig/residual_%icomp.png" % component, scale=3)


		if 0:
			### export different angle to pictures
			d = 1.5
			for az in np.arange(0, 360, 3):
				for el in np.arange(20, 70, 110):
					fig.update_layout(
						scene_camera=dict(eye=dict(
							x=d*np.cos(el/180*np.pi)*np.sin(az/180*np.pi),
							y=d*np.cos(el/180*np.pi)*np.cos(az/180*np.pi),
							z=d*np.sin(el/180*np.pi))
						)
					)

					# Save the figure to an image file using the write_image() function
					fig.write_image('image_%i_%i.png' % (el, az), width=600, height=400)

		if 0:
			from PIL import Image, ImageDraw, ImageFont
			### put picture of different angle together
			mosaic = Image.new("RGB", (600*24, 400*6))
			for i,az in enumerate(np.arange(0, 360, 15)):
				for j,el in enumerate(np.arange(10, 70, 10)):
					panel = Image.open('image_%i_%i.png' % (el, az))
					draw = ImageDraw.Draw(panel)
					#print(help(draw.text))
					draw.text((50, 100), '#%02i%02i' % (i,j), (255,0,0), font=ImageFont.truetype('Arial', 36))
					mosaic.paste(panel, (i*600, j*400))
			mosaic.save('mosaic.png')
			mosaic.close()

		if 0:
			###export a gif
			from PIL import Image
			frames = []
			for az in range(0, 360, 3):
				frames.append(Image.open('image_%i_%i.png' % (20, az)))
			frames[0].save('image.gif', format='GIF', append_images=frames[1:], save_all=True, duration=50, loop=True)



	### 3D visualization with mplot3d, a fake 3D plot
	if 0:
		from mpl_toolkits.mplot3d import axes3d

		# Plot
		fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(12,6))

		gridX, gridY, gridZ = _modelGrid()
		ax.plot_wireframe(gridX, gridY, gridZ, rcount=33, ccount=15, color='k', linewidth=0.5, zorder=2)#, rstride=1, cstride=1)
		#ax.plot_surface(gridX, gridY, gridZ, zorder=2)#, rstride=1, cstride=1)

		gridX, gridY, gridZ = _dataGrid()
		valid = np.isfinite(gridZ)
		ax.scatter(gridX[valid], gridY[valid], gridZ[valid], c=gridZ[valid], vmin=-1, vmax=1, cmap='RdYlBu_r', alpha=1, edgecolor='none', zorder=1)
		#gridX, gridY, gridResidual = _residualGrid()
		#above = gridResidual>=0
		#ax.scatter(gridX[above], gridY[above], gridZ[above], c=gridZ[above], vmin=-1, vmax=1, cmap='RdYlBu_r', alpha=1, edgecolor='none', zorder=1)
		#below = gridResidual<0
		#ax.scatter(gridX[below], gridY[below], gridZ[below], c=gridZ[below], vmin=-1, vmax=1, cmap='RdYlBu_r', alpha=0.5, edgecolor='none', zorder=3)

		ax.set(box_aspect=(1,1,.4), fc='none', 
			xticks=np.arange(-24,25,8), xlabel='X (kpc)', xlim=[-24,24],
			yticks=np.arange(-24,25,8), ylabel='Y (kpc)', ylim=[-24,24],
			zticks=np.arange(-2,2.1,1), zlabel='Z (kpc)', zlim=[-2,2])
		ax.view_init(25,155)
		ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
		ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
		ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
		plt.show()



	###R-Z map along different PHI directions
	if 0:
		step=15
		fig, axes = plt.subplots(ncols=4, nrows=4, sharex=True, sharey=True, figsize=(14,10))
		fig.subplots_adjust(hspace=0.05, wspace=0.05)
		axes = axes.ravel()
		aran = [[155. , 144.5], [144.5, 134. ], [134. , 123.5], [123.5, 113. ], \
			[113. , 102.5], [102.5,  92. ], [ 92. ,  81.5], [ 81.5,  71. ], \
			[ 71. ,  60.5], [ 60.5,  50. ], [ 50. ,  39.5], [ 39.5,  29. ], \
			[ 29. ,  18.5], [ 18.5,   8. ], [ -6. , -16. ], [-16. , -26. ]]#[9:15]
		#aran = np.vstack((np.linspace(0,360,17)[:-1], np.linspace(0,360,17)[1:])).T
		PHI[PHI > 200] -= 360
		for i,ran in enumerate(aran):
			ax = axes[i]
			index = (PHI<=ran[0]) & (PHI>ran[1])
			#print(index.sum())
			ax.scatter(R[index], Z[index], s=mass[index]/5000)

			#print(function((np.linspace(8,22,500), np.repeat(np.mean(ran),500)), best))
			dots = 500
			axis_r = np.linspace(7, 22, dots)
			axis_phi = np.repeat(np.mean(ran), dots)
			axis_z = function((axis_r, axis_phi), best)
			ax.plot(axis_r, axis_z, '--', color='red')

			ax.scatter(R[index], Z[index]-function((R[index], PHI[index]), best)-1.5, s=mass[index]/5000)
			ax.plot(axis_r, np.zeros(axis_r.size)-1.5, '--', color='grey')
			
			ax.text(0,1,'%.1f~%.1f' % tuple(ran), transform = ax.transAxes, ha='left', va='top')
			ax.set_xlim(7, 22)
			ax.set_ylim(-2.5, 1.5)
		plt.show()


	###PHI-Z map at different R distances.
	if 0:
		step=1
		fig, axes = plt.subplots(ncols=3, nrows=4, sharex=True, sharey=True, figsize=(14,10))
		fig.subplots_adjust(hspace=0.05, wspace=0.05)
		axes = axes.ravel()
		aran = [[8.5, 9.5], [9.5, 10.5], [10.5, 11.5], [11.5, 12.5], \
			[12.5, 13.5], [13.5, 14.5], [14.5, 15.5], [15.5, 16.5], \
			[16.5, 17.5], [17.5, 18.5], [18.5, 19.5], [19.5, 20.5]]
		#aran = np.vstack((np.linspace(0,360,17)[:-1], np.linspace(0,360,17)[1:])).T
		PHI[PHI > 200] -= 360
		for i,ran in enumerate(aran):
			ax = axes[i]
			index = (R>=ran[0]) & (R<ran[1])
			ax.scatter(PHI[index], Z[index], s=mass[index]/5000)

			#print(function((np.linspace(8,22,500), np.repeat(np.mean(ran),500)), best))
			dots = 500
			axis_r = np.repeat(np.mean(ran), dots)
			axis_phi = np.linspace(180, -30, dots)
			axis_z = function((axis_r, axis_phi), best)
			ax.plot(axis_phi, axis_z, '--', color='red')

			#ax.scatter(PHI[index], Z[index]-function((R[index], PHI[index]), best)-1.5, s=mass[index]/5000)
			#ax.plot(axis_phi, np.zeros(axis_r.size)-1.5, '--', color='grey')
			
			ax.text(0, 1, '%i kpc' % (np.mean(ran)), transform = ax.transAxes, ha='left', va='top')
			ax.set_xlim(180, -30)
			ax.set_ylim(-1, 1.5)


	if 0:
		f = open('residual.razm', 'w')
		C=Z-function((R, PHI), best)
		for i in range(len(R)):
			f.write('%.2f\t%.2f\t%.3f\t%.2f\n' % (R[i],PHI[i],C[i],mass[i]))      ### > r 
		f.close()
	
	plt.show()

