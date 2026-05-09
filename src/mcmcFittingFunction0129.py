import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import colormaps
import emcee, multiprocessing
#import corner
import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import circmean, circstd
from shared import textwidth, gc2g, g2gc

#np.random.seed(42)

#comp1/comp2 * dR * sin
component = 2
sn = range(0,600)#None#
excluded = True # whether to exclude clouds near 180deg
sin = False 	# set to False to fit WARP, plot WARP corner, and plot warp 3D model
			# set to True to fit SIN component, plot SIN corner, and plot 3D corrugation after subtracting warp
radwave = True
radnum = 0
ringwave = False
print('Running with %i component(s) and %s l in [195, 200] and in %s component' % (component, 'excluding' if excluded else 'including', 'corrugation' if sin else 'warp'))

if component == 1: path = 'oneComp'
elif component ==2: path = 'twoComp'
else: raise ValueError('No such component')
if excluded: path += 'Exc'

os.makedirs(path, exist_ok=True)

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
params[r'$a_1$'] = [0.09363, 'free' if not sin else 'fixed', [-0.501, 0.50]]	#a_1 in kpc
params[r'$R_{w1}$'] = [8.56849936, 'free' if not sin else 'fixed', [7, 11]]		#R_w1 in kpc
params[r'$b_{w1}$'] = [1, 'free' if (component==1) & (not sin) else 'fixed', [0.1, 2.5]]	#b_w1 index 0.9435
params[r'$\phi_{w1}$'] = (-0.7049789, 'free' if not sin else 'fixed', [-90, 90])	#phi_w1 in deg

###second warp component
params[r'$a_2$'] = (0.0, 'free' if (component==2) & (not sin) else 'fixed', [-0.50, 0.50])	#a_2 in kpc  -0.21
params[r'$R_{w2}$'] = (12.7166892, 'free' if (component==2) & (not sin) else 'fixed', [10, 17])		#R_w2 in kpc
params[r'$b_{w2}$'] = (2.03007431, 'free' if (component==2) & (not sin) else 'fixed', [0.1, 2.5])	#b_w2 index
params[r'$\phi_{w2}$'] = (7, 'free' if (component==2) & (not sin) else 'fixed', [-90, 90])	#phi_w2 in deg

###sin component
if sin and radwave:
	# Radial direction (Radcliffe wave form)
	if radnum == 0:
		#phi=31
		params[r'$A_{rad}$'] = (0.23, 'free' if sin else 'fixed', [0.01, 0.5])	# amplitude in kpc
		params[r'$R_{rad}$'] = (11, 'free' if sin else 'fixed', [9, 15])  #start of sin component in kpc
		params[r'$\sigma_{R_{rad}}$'] = (3, 'free' if sin else 'fixed', [1, 9.5]) # amplitude decay width in 1/kpc2 ([1, 6.5] for comp1)
		params[r'$P_{rad}$'] = (4., 'free' if sin else 'fixed', [3, 7.5]) # period in kpc
		params[r'$lg\gamma_{rad}$'] = (2, 'fixed', [1.4, 10]) # period decay x dmax
		params[r'$phase_{rad}$'] = (0.0, 'fixed', [-np.pi/2, np.pi]) # phase
		params[r'$A_{rad,0}$'] = (0.0, 'free' if sin else 'fixed', [-0.2, 0.5]) # phase
		# Circumferential width
		params[r'$\phi_{rad}$'] = (33.5, 'free', [30, 38])  #phi center in deg
		params[r'$\sigma_{\phi_{rad}}$'] = (4, 'free', [1, 10]) #phi width in deg
	elif radnum == 1:
		#phi=-15
		params[r'$A_{rad}$'] = (-0.2, 'free' if sin else 'fixed', [-0.5, -0.01])	# amplitude in kpc
		params[r'$R_{rad}$'] = (14, 'free' if sin else 'fixed', [8, 17])  #start of sin component in kpc
		params[r'$\sigma_{R_{rad}}$'] = (4, 'fixed' if sin else 'fixed', [1, 9.5]) # amplitude decay width in 1/kpc2 ([1, 6.5] for comp1)
		params[r'$P_{rad}$'] = (4., 'free' if sin else 'fixed', [3, 7.5]) # period in kpc
		params[r'$lg\gamma_{rad}$'] = (2, 'fixed', [1.4, 10]) # period decay x dmax
		params[r'$phase_{rad}$'] = (0.0, 'fixed', [-np.pi/2, np.pi]) # phase
		params[r'$A_{rad,0}$'] = (0.0, 'fixed' if sin else 'fixed', [-0.2, 0.5]) # phase
		# Circumferential width
		params[r'$\phi_{rad}$'] = (-18.0, 'free', [-40, -10])  #phi center in deg
		params[r'$\sigma_{\phi_{rad}}$'] = (9.5, 'fixed', [1, 10]) #phi width in deg
	elif radnum == 2:
		#phi=62
		params[r'$A_{rad}$'] = (-0.2, 'free' if sin else 'fixed', [-0.5, -0.01])	# amplitude in kpc
		params[r'$R_{rad}$'] = (11, 'free' if sin else 'fixed', [9.5, 12])  #start of sin component in kpc
		params[r'$\sigma_{R_{rad}}$'] = (4, 'fixed' if sin else 'fixed', [1, 8.5]) # amplitude decay width in 1/kpc2 ([1, 6.5] for comp1)
		params[r'$P_{rad}$'] = (4., 'free' if sin else 'fixed', [2, 7]) # period in kpc
		params[r'$lg\gamma_{rad}$'] = (2, 'fixed', [1.4, 10]) # period decay x dmax
		params[r'$phase_{rad}$'] = (0.0, 'fixed', [-np.pi/2, np.pi]) # phase
		params[r'$A_{rad,0}$'] = (0.0, 'fixed' if sin else 'fixed', [-0.2, 0.5]) # phase
		# Circumferential width
		params[r'$\phi_{rad}$'] = (65, 'free', [45, 85])  #phi center in deg
		params[r'$\sigma_{\phi_{rad}}$'] = (9.5, 'fixed', [1, 10]) #phi width in deg
	elif radnum == 3:
		#phi=106, fail
		params[r'$A_{rad}$'] = (-0.2, 'free' if sin else 'fixed', [-0.5, -0.01])	# amplitude in kpc
		params[r'$R_{rad}$'] = (10, 'free' if sin else 'fixed', [9.5, 12])  #start of sin component in kpc
		params[r'$\sigma_{R_{rad}}$'] = (4, 'fixed' if sin else 'fixed', [1, 6.5]) # amplitude decay width in 1/kpc2 ([1, 6.5] for comp1)
		params[r'$P_{rad}$'] = (4., 'free' if sin else 'fixed', [2, 5]) # period in kpc
		params[r'$lg\gamma_{rad}$'] = (2, 'fixed', [1.4, 10]) # period decay x dmax
		params[r'$phase_{rad}$'] = (0.0, 'fixed', [-np.pi/2, np.pi]) # phase
		params[r'$A_{rad,0}$'] = (0.0, 'fixed' if sin else 'fixed', [-0.2, 0.5]) # phase
		# Circumferential width
		params[r'$\phi_{rad}$'] = (100, 'free', [90, 110])  #phi center in deg
		params[r'$\sigma_{\phi_{rad}}$'] = (9.5, 'fixed', [1, 11]) #phi width in deg
	elif radnum == 4:
		#phi=106, fail
		params[r'$A_{rad}$'] = (0.2, 'free' if sin else 'fixed', [0.01, 0.6])	# amplitude in kpc
		params[r'$R_{rad}$'] = (12, 'free' if sin else 'fixed', [9.5, 15])  #start of sin component in kpc
		params[r'$\sigma_{R_{rad}}$'] = (2, 'fixed' if sin else 'fixed', [1, 6.5]) # amplitude decay width in 1/kpc2 ([1, 6.5] for comp1)
		params[r'$P_{rad}$'] = (4., 'free' if sin else 'fixed', [1, 9]) # period in kpc
		params[r'$lg\gamma_{rad}$'] = (2, 'fixed', [1.4, 10]) # period decay x dmax
		params[r'$phase_{rad}$'] = (0.0, 'fixed', [-np.pi/2, np.pi]) # phase
		params[r'$A_{rad,0}$'] = (0.0, 'fixed' if sin else 'fixed', [-0.2, 0.5]) # phase
		# Circumferential width
		params[r'$\phi_{rad}$'] = (140, 'free', [130, 155])  #phi center in deg
		params[r'$\sigma_{\phi_{rad}}$'] = (9.5, 'fixed', [1, 11]) #phi width in deg

elif sin and ringwave:
	# Circumferential direction
	params[r'$A_{az}$'] = (0.2, 'free' if sin else 'fixed', [0.08, 0.5])	# amplitude in kpc
	params[r'$\phi_{az}$'] = (20., 'free' if sin else 'fixed', [0, 40]) # period in kpc
	params[r'$\sigma_{\phi_{az}}$'] = (100., 'fixed' if sin else 'fixed', [1, 200]) # period in kpc
	params[r'$P_{az}$'] = (52., 'free' if sin else 'fixed', [40, 90]) # period in kpc
	params[r'$lg\gamma_{az}$'] = (2.8, 'fixed' if sin else 'fixed', [1, 4]) # period in kpc
	# Radial width
	params[r'$R_{az}$'] = (12.7, 'free' if sin else 'fixed', [11, 15])  #start of sin component in kpc
	params[r'$\sigma_{R_{az}}$'] = (0.8, 'free' if sin else 'fixed', [0.1, 2.5]) # amplitude decay width in 1/kpc2 ([1, 6.5] for comp1)

else:
	params[r'$a_3$'] = (0.3, 'free' if sin else 'fixed', [0, 0.50]) #a_3 in kpc
	params[r'$R_{c}$'] = (7.54, 'free' if sin else 'fixed', [6, 9])  #start of sin component in kpc
	params[r'$P_0$'] = (1.9, 'free' if sin else 'fixed', [1.01, 5.0])  #period of sin component in kpc
	params[r'$P_1$'] = (0.2, 'free' if sin else 'fixed', [-0.2, 0.3])  #period increasement of sin component in kpc
	# Radial linear baseline 
	params[r'$a_{rad,0}$'] = (0.00, 'free' if sin else 'fixed', [-0.2, 0.2])  #baseline
	# Circumferential width
	params[r'$\phi_{circ}$'] = (33.5, 'fixed', [30, 38])  #phi center in deg
	params[r'$\sigma_{circ}$'] = (4, 'fixed', [1, 10]) #phi width in deg


#params[r'$\phi_{circ}$'] = (36.5, 'fixed', [28, 50])  #phi center in deg
#params[r'$\phi_{circ,width}$'] = (6, 'fixed', [1, 10]) #phi width in deg

#[1.39573585e-01, 7.74237334, 2.86256551, 1.57132833e-01,3.21995180e+01, 9.11163275, 1.06050455e-02]
#params[r'$a_5$'] = (0.0102, 'free', [0.001, 6])		#period of sin component in kpc
## [-0.07616907  0.10783739  7.79992154  0.9348509  -9.02005785]   constant xf
##sin [1.36492974e-01 8.47536353e+00 1.22706616e+00 4.68903037e-01 2.94463588e+01 9.84374395e+00 7.01067072e-04]
## [-0.069754    0.13771357  8.12746402  0.78466442 -8.47120866]   ##xfactor
###MCMC setting  [ -0.08204349   0.10188      7.77984319 -10.99553919  -0.11923537 13.79306634 -23.78646842]

nwalkers = 32 	#how many thread
burnin = 500#1500 	#number of iteration to reach local minimum (usually 20~30% of niter)
niter = 1000#3000	#number of iteration

'''
one component: fix a0=0, free (a1, Rw1, bw1, phiw1)
two component: fix a0=0, bw1=1, free (a1, Rw1, phiw1, a2, Rw2, bw2, phiw2)
phiw at diff R: free (phiW1) only
sin component: fix phisin=33.5, phisinwidth=4, free (a3, Rsin, P0, P1, a4)
'''


params[r'$a_1$'][1] = 'fixed'
params[r'$R_{w1}$'][1] = 'fixed'
params[r'$b_{w1}$'][1] = 'fixed'



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


def radcliffe_wave(d, A, delta, P, dmaxXgamma, phi, a4):
	'''
	Exact Radcliffe wave form from the paper:
	Δz(d) = A × exp[-δ×(d/kpc)²] × sin[(2πd/P) × (1 + d/(d_max×γ)) + φ]
	
	Parameters:
	d: distance from start point (kpc)
	A: amplitude
	delta: amplitude decay rate  
	P: period (kpc)
	gamma: period decay rate
	phi: phase
	d_max: maximum distance (kpc)
	'''
	wav = (A * np.sin((2 * np.pi * d / P) * (1 + d / dmaxXgamma *0) + phi) + a4)
	#amp = np.exp(-d**2/2/delta**2)
	amp = np.exp(-d/delta)
	return amp*wav


def function(x, free_params, warp=not sin, sin=sin):
	if sin and radwave:
		a0, a1, Rw1, bw1, PHIw1, a2, Rw2, bw2, PHIw2, Arad, Rrad, sigmaRrad, Period0, Period1, phase, Arad0, PHIrad, sigmaPHIrad = fp2p(free_params)
	elif sin and ringwave:
		a0, a1, Rw1, bw1, PHIw1, a2, Rw2, bw2, PHIw2, Acirc, PHIcirc, sigmaPHIcirc, Period0, Period1, Rcirc, sigmaRcirc = fp2p(free_params)
	else:
		a0, a1, Rw1, bw1, PHIw1, a2, Rw2, bw2, PHIw2, Arad, Rrad, Period0, Period1, Arad0, PHIcirc, sigmaPHIcirc = fp2p(free_params)

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
		if radwave:
			# radial
			Zrad = radcliffe_wave((R-Rrad), Arad, sigmaRrad, Period0, Period1, phase, Arad0)
			# circumferential
			Zaz = np.exp(-(PHI-PHIrad)**2/ 2 / sigmaPHIrad**2)
			y += Zrad * Zaz

		elif ringwave:
			Zrad = Acirc * np.sin(2 * np.pi * (PHI-PHIcirc) / Period0 * (1 + (PHI-PHIcirc)/10**Period1 *0))# * np.exp(-(PHI-PHIcirc)**2/2/sigmaPHIcirc**2)
			Zaz = np.exp(-(R-Rcirc)**2/ 2 / sigmaRcirc**2)
			y += Zrad * Zaz

		else:
			###sin component
			index = R>Rw1#Rrad
			###y[index] += ( a4*(R[index] - Rrad) + Arad*np.sin( (R[index] - Rrad) / (Period0 + Period1*(R[index]-Rrad)) * 2*np.pi ) ) * np.exp(-(PHI[index]-PHIcirc)**2/ 2 / PHIcircwid**2)
			# radial linear
			Zlin = Arad0 * (R[index] - Rrad)
			# radial sinusoidal
			Zsin = Arad * np.sin( (R[index] - Rrad) / (Period0 + Period1*(R[index]-Rrad)) * 2*np.pi )
			# circumferential width
			Zaz = np.exp(-(PHI[index]-PHIcirc)**2/ 2 / sigmaPHIcirc**2)
			y[index] += (Zlin+Zsin) * Zaz

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

	with multiprocessing.Pool() as pool:
		#mcmc sampler
		sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=data, \
			moves=emcee.moves.StretchMove(0.2), pool=pool)
		
		print("Initial = ", initial)

		print("Running burn-in...")
		p0, _, _ = sampler.run_mcmc(p0, burnin, progress=True)

		#print("After burn-in:\n")
		#for p in p0: print(' '.join(['%+.3f' % v for v in p]))

		sampler.reset()

		print("Running production...")
		pos, prob, state = sampler.run_mcmc(p0, niter, progress=True)
	
	return sampler, pos, prob, state


def _excludedXYZ(x, y, z):
	### remove excluded surface
	r, phi = XY2RPHI(x, y)
	l, b, d = gc2g(r, phi, z)
	idx = (l>165) & (l<199) & (x**2+y**2>16**2)
	x[idx] = np.nan
	y[idx] = np.nan
	z[idx] = np.nan
	return x, y, z



def convolveZ(x, y, z, w, sampleX, sampleY, kernel=1, useMask=True, generateMask=False):
	### sample and convolve on xy points
	### generate a high resolution mask for edge (generateMask=True)
	### apply mask (useMask=True)
	sigma = kernel/np.sqrt(8*np.log(2))
	clipRadius = 4*sigma
	tree = cKDTree(np.column_stack((x, y)))
	samplePoints = np.column_stack((sampleX.ravel(), sampleY.ravel()))  # flatten grid
	neighbors = tree.query_ball_point(samplePoints, r=clipRadius)

	### sample x,y
	sampleZ = np.full(sampleX.size, np.nan)
	for j, idxs in enumerate(neighbors):
		if len(idxs) == 0:
			continue
		dx = samplePoints[j,0] - x[idxs]
		dy = samplePoints[j,1] - y[idxs]
		d2 = dx**2 + dy**2
		kw = np.exp(-d2 / (2 * sigma**2))
		w_eff = w[idxs] * kw
		denom = np.sum(w_eff)
		if denom > 0:
			sampleZ[j] = np.sum(z[idxs] * w_eff) / denom


	### generate edge mask on high-resolution grid
	if useMask:
		fineStep = 0.1 # in kpc

		### finer grid
		axisX = np.arange(-10, 18, fineStep)
		axisY = np.arange(-15, 24, fineStep)

		if generateMask:
			# mesh grid
			gridX, gridY = np.meshgrid(axisX, axisY)
			# Flatten grid for KDTree query
			gridPoints = np.column_stack((gridX.ravel(), gridY.ravel()))

			# Build KDTree of data points
			tree = cKDTree(np.column_stack((x, y)))
			kernelRadius = kernel/2 # in kpc

			# Count neighbors within kernel
			neighbors = tree.query_ball_point(gridPoints, r=kernelRadius)
			numDen = np.array([len(idx) for idx in neighbors])

			# Reshape to grid
			numDen = numDen.reshape(gridX.shape)

			# >5 clouds in kernel
			mask = gaussian_filter(numDen, sigma=1/fineStep) >=5

			# mask 8 kpc ring
			smoothZero = gaussian_filter((numDen>0).astype(float), sigma=1/fineStep)
			mask[smoothZero<0.5] = False

			# mask l and Rgal range
			tol=3 #deg
			gridR, gridPHI = XY2RPHI(gridX, gridY)
			gridL, gridB, gridD = gc2g(gridR, gridPHI, np.zeros(gridR.shape))
			mask[(gridL<15-tol) | (gridL>229.75+tol)] = False # cloud range
			mask[(gridL>165+tol) & (gridL<195-tol)] = False # Kdist range
			mask[(gridL>165+tol) & (gridL<199-tol) & (gridR>16+0.2)] = False
			mask[gridR > 17] = False # Rgal range

			np.save('highres_mask.npy', mask)
		else:
			mask = np.load('highres_mask.npy')

		interp = RegularGridInterpolator((axisY, axisX), mask, method='linear', bounds_error=False, fill_value=0)
		sampleMask = interp(samplePoints[:,::-1])>0.5 #xy order need reverse here
		sampleZ[~sampleMask] = np.nan

	sampleZ = sampleZ.reshape(sampleX.shape)
	return sampleZ



def convolveZonGrid(x, y, z, w, step=None, polar=True, **kws):
	'''
	convolveZ on regular grid
	step: grid step for (dx, dy) or (dr, dphi)
	polar: grid in xy or polar
	'''
	if polar:
		if step is None: step=(.25, 2)	#
		gridR, gridPHI = np.meshgrid(np.arange(1, 22.1, step[0]), np.arange(0, 360.1, step[1]))
		gridX, gridY = RPHI2XY(gridR, gridPHI)
	else:
		if step is None: step=(.2, .2)
		gridX, gridY = np.meshgrid(np.arange(-10, 18, step[0]), np.arange(-15, 24, step[1]))
	
	gridZ = convolveZ(x, y, z, w, gridX, gridY, **kws)
	'''
	if not residual:
		gridZ = convolveZ(X, Y, Z, mass, gridX, gridY, **kws)
	else:
		gridZ = convolveZ(X, Y, Z-function((R, PHI), bestmed, warp=warp, sin=sin), mass, gridX, gridY, **kws)
	'''
	return gridX, gridY, gridZ



def XY2RPHI(x, y):
	x = np.asarray(x)
	y = np.asarray(y)
	r = np.sqrt(x**2 + y**2)
	phi = np.arctan2(x, y)/np.pi*180
	return r, phi

def RPHI2XY(r, phi):
	r = np.asarray(r)
	phi = np.asarray(phi)
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

def g2gc(l, b, d):
	#galactic to galactocentric
	d2r = np.pi/180
	X = d * np.cos(b*d2r) * np.cos(l*d2r) - 8.15
	Y = d * np.cos(b*d2r) * np.sin(l*d2r)
	Z = d * np.sin(b*d2r)

	Rgal = np.sqrt(X**2 + Y**2)
	Az = np.arctan2(Y, -X) / d2r #az=0 toward sun and counterclockwise
	return Rgal, Az, Z



if __name__ == '__main__':
	if 1:
		###load observation

		# exclude l within 100~136 on perseus arm.
		l,b,v,lwidth,distance,err_dist,size_pc,mas,\
		complete,surface_density,viral_para,gal_ridus,\
		xx,yy,zz,index,l_rms,b_rms,v_rms,beta,T_int,area,tpeak,tex,nh2,pint,fill,lp,bp,vp = \
		np.loadtxt('per_para.txt').T
		#idx = (l<100) | (l>136)
		per = np.loadtxt('per_para.txt')#[idx]

		l,b,v,lwidth,distance,err_dist,size_pc,mas,\
		complete,surface_density,viral_para,gal_ridus,\
		xx,yy,zz,index,l_rms,b_rms,v_rms,beta,T_int,area,tpeak,tex,nh2,pint,fill,lp,bp,vp = \
		np.hstack((np.loadtxt('osc2_para.txt').T, np.loadtxt('out_para.txt').T, per.T))
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
			ind = (mas > 0) & ~((l>195) & (l<199) & (xx**2+yy**2>16**2))# & (xx**2+yy**2<16**2)# & (mas>1e4))
		else:
			ind = (mas > 0)
		print('Total:', ind.size, 'After excluding:', ind.sum())
		X = xx[ind]
		Y = yy[ind]
		Z = zz[ind]
		mass = mas[ind] #/complete[ind]) #*gal_ridus[ind]
		norm_mass = mass/np.sum(mass)

		# to get statistical uncertainty, i need to recalculate R,Phi,z after adding noise
		l = l[ind]
		b = b[ind]
		distance = distance[ind]
		err_dist = err_dist[ind]

		R, PHI = XY2RPHI(X, Y)

		# estimate uncertainty
		#dR_ds = (distance * np.cos(np.deg2rad(b))**2 - 8.15 * np.cos(np.deg2rad(l)) * np.cos(np.deg2rad(b))) / R
		#dR = (np.sin(np.deg2rad(b)) - dz_dR * dR_ds)**2 * err_dist**2

		#data = data.T[mass>1e3].T
		'''
		data = np.loadtxt(cat).T
		R, PHI, Z, mass = data

		X, Y = RPHI2XY(R, PHI)
		'''

		### filter R range
		if 0:
			### warp model of annulus
			i = 1
			R_sep = [9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
			R1, R2 = R_sep[i:i+2]
			print(R1, R2)
			idx = (R>=R1) & (R<R2)
			R = R[idx]
			PHI = PHI[idx]
			Z = Z[idx]
			mass = mass[idx]
			norm_mass = mass/np.sum(mass)
			l = l[idx]
			b = b[idx]
			distance = distance[idx]
			err_dist = err_dist[idx]
			suffix = '_%i_%i' % (R1, R2)

		elif sin:
			### sin component
			if radwave:
				suffix = '_radwave%i' % radnum
			elif ringwave:
				suffix = '_ringwave'
			else:
				suffix = '_sin'
		else:
			### warp model
			suffix = ''
			#suffix = '_errD_mock0957'
			#suffix = '_exPer100136'
			#suffix = '_le16kpc'

		data = np.vstack([R, PHI, Z, mass])

		### load residual data
		if sin:
			data = np.load('residual_%icomp.npy' % (component))
			#data[3] = data[3] * data[0]**2

			R, PHI, Z, mass = data
			X, Y = RPHI2XY(R, PHI)
			print(Z)
			
		print('Data in shape:', data.shape)
	else:
		### fit median
		if sin:
			suffix = '_fitmedian'
			R, PHI, Z, mass = np.load('residual_%icomp.npy' % (component))
			X, Y = RPHI2XY(R, PHI)
			gridX, gridY, gridZ = convolveZonGrid(X, Y, Z, mass, step=(0.2, 0.2), polar=False, useMask=True, generateMask=True)
			gridR, gridPHI = XY2RPHI(gridX, gridY)
			gridMass = np.ones_like(gridZ)
			idx = np.isfinite(gridZ)
			#plt.scatter(gridX[idx], gridY[idx], c=gridZ[idx], cmap='coolwarm', vmin=-0.3, vmax=0.3, s=10)
			#plt.show()
			data = np.vstack([gridR[idx], gridPHI[idx], gridZ[idx], gridMass[idx]])
		'''
		suffix = '_star'
		from astropy.io import ascii
		t = ascii.read("star/apjs.webarchive")
		l = np.array(t['glon'])
		b = np.array(t['glat'])
		d = np.array(t['d-all'])/1e3
		w = mass = 1/(np.array(t['e_d-all'])/1e3)**2
		from shared import g2gc
		X, Y, Z = g2gc(l, b, d, xyz=True)
		R, PHI, Z = g2gc(l, b, d, xyz=False)
		norm_mass = mass/np.sum(mass)

		data = np.vstack([R, PHI, Z, mass])

		if sin:
			data = np.load('residual_%icomp%s.npy' % (component, suffix))
			R, PHI, Z, mass = data
			X, Y = RPHI2XY(R, PHI)
		'''

		'''
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
		'''

	###run MCMC
	if 0:
		if sn:
			for i in sn:
				print('iter %04i' % i)
				R, PHI, Z = g2gc(l,b,distance+np.random.normal(0,1,err_dist.size)*err_dist/5*7)#err_dist=dD/5km/s, convert to 7 km/s
				data = np.vstack([R, PHI, Z, mass])

				if 1:
					# R-phi
					R_sep = [9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
					for s in range(9):
						R1, R2 = R_sep[s:s+2]
						print('dealing range:', R1, R2)
						idx = (R>=R1) & (R<R2)
						data = np.vstack([R[idx], PHI[idx], Z[idx], mass[idx]])
						sampler, pos, prob, state = main(data)
						steps = sampler.flatchain
						probs = sampler.flatlnprobability
						suffix = '_%i_%i' % (R1, R2)
						np.save(os.path.join(path, 'steps%s_errD_mock%04i.npy' % (suffix, i)), steps)
						np.save(os.path.join(path, 'probs%s_errD_mock%04i.npy' % (suffix, i)), probs)
				else:
					# all samples
					sampler, pos, prob, state = main(data)
					steps = sampler.flatchain
					probs = sampler.flatlnprobability
					np.save(os.path.join(path, 'steps%s_errD_mock%04i.npy' % (suffix, i)), steps)
					np.save(os.path.join(path, 'probs%s_errD_mock%04i.npy' % (suffix, i)), probs)
		else:
			sampler, pos, prob, state = main(data)
			steps = sampler.flatchain
			probs = sampler.flatlnprobability
			np.save(os.path.join(path, 'steps%s.npy' % suffix), steps)
			np.save(os.path.join(path, 'probs%s.npy' % suffix), probs)
	else:
		steps = np.load(os.path.join(path, 'steps%s.npy' % suffix))
		probs = np.load(os.path.join(path, 'probs%s.npy' % suffix))

	### get median / best / best med
	if 1:
		def betterOutput(p):
			i = 0
			out = ['','']
			for k in params:
				if params[k][1] == 'free':
					out[0] += '%15s' % k
					out[1] += '%15.4f' % p[i]
					i+=1
			out[0]+='\n'
			return out
		print(probs.shape, steps.shape)
		med  = np.median(steps, axis=0)
		best = steps[np.argmax(probs)]
		bestmed = np.mean(steps[probs > np.percentile(probs, 98)], axis=0)
		print('The best fitting result is:\n', *betterOutput(best))
		print('The median fitting result is:\n', *betterOutput(med))
		print('The median of the best is:\n', *betterOutput(bestmed))

		if sn:
			mockmed = []
			for i in sn:
				if not os.path.exists(os.path.join(path, 'steps%s_errD_mock%04i.npy' % (suffix, i))): continue
				mocksteps = np.load(os.path.join(path, 'steps%s_errD_mock%04i.npy' % (suffix, i)))
				mockprobs = np.load(os.path.join(path, 'probs%s_errD_mock%04i.npy' % (suffix, i)))
				if np.isnan(mocksteps).any(): continue
				bestmedi = np.mean(mocksteps[mockprobs > np.percentile(mockprobs, 98)], axis=0)
				mockmed.append(bestmedi)
			mockmed = np.array(mockmed)

			# mocks give positive a2 and phi_w2+90 results, regulate them
			if component==2:
				circ = mockmed[:,6] > 25
				mockmed[circ, 3] *= -1	# reverse a2
				mockmed[circ, 6] -= 90	# move phi_w2 180/2
				mockmed = mockmed[(11.5<mockmed[:,4]) & (mockmed[:,4]<16.45)]
				#mockmed = mockmed[mockmed[:,5]>0.5]
			print(mockmed.shape, np.mean(mockmed))
			np.save('mockmed.npy', mockmed)
			mockmean = np.mean(mockmed, axis=0)
			bias = mockmean-bestmed
			p16, p84 = np.percentile(mockmed, [15.86, 84.14], axis=0)
			upper = p84 - mockmean
			lower = mockmean - p16
			err_upper = np.sqrt(upper**2 + np.max([np.zeros(bias.size), +bias], axis=0)**2)
			err_lower = np.sqrt(lower**2 + np.max([np.zeros(bias.size), -bias], axis=0)**2)
			print('error + is:\n', *betterOutput(err_upper))
			print('error - is:\n', *betterOutput(-err_lower))
			
	'''
	#bestmed=initial
	fig, ax = plt.subplots()
	idx = np.abs(PHI-bestmed[-2]) < 4
	ax.scatter(R[idx], Z[idx], s=(mass[idx]*8e-3), c=PHI[idx], cmap='coolwarm', alpha=0.5)

	Raxis = np.linspace(7, 22, 1000)
	PHIaxis = np.full(1000, bestmed[-2])
	mz = function((Raxis, PHIaxis), bestmed)
	ax.plot(Raxis, mz, 'r--')

	fig, ax = plt.subplots()
	Xg, Yg = np.meshgrid(np.linspace(-20,20,300), np.linspace(-20,20,300))
	Rg, PHIg = XY2RPHI(Xg, Yg)
	mz = function((Rg, PHIg), bestmed)
	ax.imshow(mz, extent=(-20,20,-20,20), origin='lower', cmap='coolwarm', vmin=-0.3, vmax=0.3)
	from shared import rad_sep
	for r in rad_sep: ax.plot([0,np.sin(np.deg2rad(r))*20],[0,np.cos(np.deg2rad(r))*20], 'k--', alpha=0.3)
	for r in np.arange(2, 20, 2): ax.plot(*RPHI2XY(np.full(360,r), np.linspace(0,360,360)), 'k--', alpha=0.3)
	#plt.show()
	'''

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
		#aicbic(data, bestmed)


	### corner plot
	if 0:
		from shared import textwidth, subfigureIndexFont
		figscale = 0.6 if not sin else 0.45
		figwidth = textwidth*figscale
		plt.rcParams['savefig.dpi'] = 280
		plt.rcParams['axes.linewidth'] = 0.8
		plt.rcParams['axes.labelsize'] = 20
		plt.rcParams['axes.labelweight'] = 'bold'
		plt.rcParams['xtick.labelsize'] = 18
		plt.rcParams['ytick.labelsize'] = 18
		plt.rcParams['xtick.direction'] = 'in'
		plt.rcParams['ytick.direction'] = 'in'
		plt.rcParams['xtick.top'] = True
		plt.rcParams['ytick.right'] = True
		plt.rcParams['xtick.minor.visible'] = True
		plt.rcParams['ytick.minor.visible'] = True
		plt.rcParams['legend.fontsize'] = 15

		def corners(steps, probs, bins=10, showtop=0.9, order=None,
			figure_kws = dict(figsize=(figwidth, figwidth)),
			hist_kws = dict(density=True, color = 'grey', histtype='step'),
			values = None,
			hist_value_kws = dict(pos='top', marker='v', color='tab:red'),
			hist2d_value_kws = dict(marker='+', color='tab:red', markersize=8, markeredgewidth=2),
			limits = [0.1586, 0.8414],
			hist_limit_kws = dict(linestyle='--', color='k', linewidth=1, alpha=0.5),
			labels = None,
			hist_label_kws = dict(fontsize=20),
			hist2d = True,
			hist2d_kws = dict(cmap='Greys'), 
			scatter = False,
			scatter_kws = dict(s=0.2, cmap='RdYlBu', alpha=0.3),
			contour = False,
			contour_smooth_sigma = 0,
			contour_kws = dict(colors='gray', levels=[0.2, 0.4, 0.6, 0.8], linewidths=1.5, alpha=0.8),
			output = 'test.out',
			):

			ndim = steps.shape[1]
			### plot corners
			def autoLimit(s):
				### extend according to s range
				ma, mi = np.nanmax(s), np.nanmin(s)
				return mi-(ma-mi)*0.1, ma+(ma-mi)*0.1
			def scientificNotationTitle(v, l, u, digit=4):
				### show value, lower, and upper uncertainty in scientific notation format
				e = np.floor(np.log10(np.abs(v)))
				### dont convert xx.xx to x.xxx*10^1
				if e == 1:
					e = 0
					mind = 2
				else: mind = 3
				vs = v/10**e
				us = u/10**e
				ls = l/10**e
				### show more digit if any error is zero
				for d in range(mind, 11):
					ustr = f'%+.{d}f' % us
					lstr = f'%+.{d}f' % ls
					if float(ustr) != 0.0 and float(lstr) != 0.0: break
					else: d = 10
				vstr = f'%+.{d}f' % vs
				ustr = f'%+.{d}f' % us
				lstr = f'%+.{d}f' % ls
				### print latex form
				vlatex = f'%+.{d-int(e)}f' % (v)
				ulatex = f'%+.{d-int(e)}f' % (u)
				llatex = f'%+.{d-int(e)}f' % (l)
				print('& %s & $^{%s}_{%s}$' % (vlatex, ulatex, llatex))#export latex table
				if e==0:
					return '$%s_{%s}^{%s}$' % (vstr, lstr, ustr)
				else:
					return '$%s_{%s}^{%s}$ x$10^{%i}$' % (vstr, lstr, ustr, e)

			### sort steps (or high prob steps might be put behind)
			steps = steps[np.argsort(probs)]
			probs = np.sort(probs)		
			### only plot the top% steps, clip low prob steps.
			topmask = probs >= np.quantile(probs, showtop)

			### get a square figure.
			fig, ax = plt.subplots(nrows=ndim, ncols=ndim, **figure_kws)
			if ndim == 1: ax = np.array(ax).reshape(1,1)
			plt.subplots_adjust(bottom=0.09, top=0.91, left=0.09, right=0.91, wspace=0.05, hspace=0.05)

			### decide where to put marker in hist
			valuePos = hist_value_kws.pop('pos', 'bartop')
			contourLevels = np.array(contour_kws.pop('levels', [0.2, 0.4, 0.6, 0.8]))

			### adjust order
			if order is not None:
				steps = steps[:, order]
				if values is not None: values = [values[o] for o in order]
				if labels is not None: labels = [labels[o] for o in order]

			### output values and limits
			f = open(output, 'w')
			for i in range(ndim):
				### hist at [i,i]
				h,xe,_ = ax[i,i].hist(steps[:,i][topmask], bins=bins, **hist_kws)

				### values
				if values is not None:
					if valuePos == 'bartop': 
						upperY = h[np.searchsorted(xe, values[i], side='right')-1] * 1.1
					elif valuePos == 'top':
						upperY = np.max(h)
				ax[i,i].plot(values[i], upperY, **hist_value_kws)

				### limits
				lim = np.quantile(steps[:,i], limits)
				for v in lim:
					ax[i,i].axvline(v, **hist_limit_kws)

				### labels
				if labels is not None:
					text = '%s=%s' % (labels[i], scientificNotationTitle(values[i], min(lim)-values[i], max(lim)-values[i]))
					ax[i,i].text(0, 1.02, text, transform=ax[i,i].transAxes, ha='left', va='bottom', **hist_label_kws)

				### output
				if (labels is not None) and (values is not None):
					f.write('%s %f %f %f\n' % (labels[i], values[i], min(lim)-values[i], max(lim)-values[i]))

				### hide x ticklabels except the bottom row
				if i<ndim-1: ax[i,i].xaxis.set_ticklabels([])
				### hide all y tickslables
				ax[i,i].set_yticklabels([])

				### ticks
				ax[i,i].set_autoscale_on(False)
				ax[i,i].minorticks_on()
				ax[i,i].tick_params(top=True, left=False, direction='in', labelrotation=45, length=5)
				ax[i,i].tick_params(which='minor', top=True, left=False, direction='in', length=2)
				ax[i,i].set_xlim(autoLimit(steps[:,i][topmask]))

				### plot hist2d
				for j in range(i):
					h, xe, ye = np.histogram2d(steps[:,j][topmask], steps[:,i][topmask], bins=bins)
					### hist2d
					if hist2d:
						ax[i,j].imshow(h.T, origin='lower', extent=[xe[0], xe[-1], ye[0], ye[-1]], aspect='auto', **hist2d_kws)
					### scatter
					if scatter:
						ax[i,j].scatter(steps[:,j][topmask], steps[:,i][topmask], c=probs[topmask], **scatter_kws)
					### contour
					if contour:
						if contour_smooth_sigma is not None:
							h = gaussian_filter(h, sigma=contour_smooth_sigma)
						xc = (xe[1:]+xe[:-1])/2
						yc = (ye[1:]+ye[:-1])/2
						ax[i,j].contour(xc, yc, h.T, levels=contourLevels*np.max(h), **contour_kws)

					### values
					ax[i,j].plot(values[j], values[i], **hist2d_value_kws)
					
					### hide x ticklabels except the bottom row
					if i<ndim-1: ax[i,j].xaxis.set_ticklabels([])
					### hide y ticklabels except the left column
					if j>0: ax[i,j].yaxis.set_ticklabels([])

					### ticks
					ax[i,i].set_autoscale_on(False)
					ax[i,j].minorticks_on()
					ax[i,j].set_yticks(ax[i,i].get_xticks()) ### make sure ytick is the same as the [i,i] xtick
					ax[i,j].tick_params(top=True, right=True, direction='in', labelrotation=45, length=5)
					ax[i,j].tick_params(which='minor', top=True, right=True, direction='in', length=2)
					ax[i,j].set_xlim(autoLimit(steps[:,j][topmask]))
					ax[i,j].set_ylim(autoLimit(steps[:,i][topmask]))

					### hide upper right
					ax[j,i].axes.set_axis_off()
			f.close()
			'''
			if sin:
				if component==1:
					if radwave:
						ax[0,0].text(-0.25, 1, 'a', ha='left', va='top', color='black', font=subfigureIndexFont, transform=ax[0,0].transAxes)
					else:
						ax[0,0].text(-0.25, 1, 'b', ha='left', va='top', color='black', font=subfigureIndexFont, transform=ax[0,0].transAxes)
				else:
					if radwave:
						ax[0,0].text(-0.32, 1, 'c', ha='left', va='top', color='black', font=subfigureIndexFont, transform=ax[0,0].transAxes)
					else:
						ax[0,0].text(-0.32, 1, 'd', ha='left', va='top', color='black', font=subfigureIndexFont, transform=ax[0,0].transAxes)
			'''
			return fig, ax

		if sin:
			if radwave:
				if radnum==0: order = (0,3,4,1,2,5,6)
				else: order = (0,2,1,3)#(0,3,1,2,4,5)
			elif ringwave: order = (0,2,1,3,4)
			else: order = None
		else: order = None

		if sn:
			fig, ax = corners(mockmed, np.ones(mockmed.shape[0]), bins=31, labels=free_params_name, values=bestmed, showtop=1)
		else:
			fig, ax = corners(steps, probs, bins=31, labels=free_params_name, values=bestmed, showtop=0.5, order=order)

		'''
		#avoid ticklabel overlap
		if (component==1) and sin and radwave and (radnum==0):
			for i in range(7): ax[i, 0].set_xlim(0.1973882544240669, 0.20249999)
		elif (component==2) and sin and radwave and (radnum==0):
			for i in range(7): ax[i, 1].set_xlim(4.4664515570296555, 4.537)
			for i in range(7): ax[i, 2].set_xlim(0.0776, 0.0819999999)
			for i in range(2): ax[2, i].set_ylim(0.07794477178721429, 0.08199999)
		plt.savefig('fig/corner_%icomp%s.%s' % (component, suffix, plt.rcParams['savefig.format']), bbox_inches='tight')
		plt.savefig('fig/corner_%icomp%s.png' % (component, suffix), bbox_inches='tight')
		#corner.corner(steps, show_titles=True, plot_contours=False, plot_datapoints=True, quantiles=[0.16, 0.5, 0.84], range=None)#[(v-0.5,v+0.5) for v in best])
		'''

	### export residual, only normal warp
	if not sin:
		data[2] -= function(data[:2], bestmed, warp=True, sin=False)
		#print(data[:2], function(data[:2], bestmed, warp=True, sin=False))
		np.save('residual_%icomp%s.npy' % (component, suffix), data)
		print('Export to residual_%icomp%s.npy' % (component, suffix))


	###prepare functions for 3D visualization
	if False:
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



	### plot face-on gplane
	if 0:
		def gplane(fig, ax, im=None, step=0.25, kernel=1.2, residual=False, arm=True, grid=True, draw_phi=True, rotation_arrow=False, \
			colorbar_kws={'rect':[0.28, 0.2, 0.02, 0.22], 'orientation':'vertical'}, **kws):

			if im is None:
				# if sin==False, use convolveZonGrid(X, Y, Z-Z_w, mass, ....
				from shared import function_warp, p_1comp, p_2comp, subfigureIndexFont, function_PoggioWarp
				#Z_w = function_warp((R,PHI), p=p_1comp if component==1 else p_2comp)
				#Z_w = function_PoggioWarp((R,PHI), cepheids=True, straight=True)

				Z_sin = function((R, PHI), bestmed)

				gridX, gridY, gridZ = convolveZonGrid(X, Y, Z, mass, step=(step, step), polar=False, useMask=True, generateMask=False)
				extent=(gridX[0,0]-step/2, gridX[0,-1]+step/2, gridY[0,0]-step/2, gridY[-1,0]+step/2)
				im = ax.imshow(gridZ, origin='lower', extent=extent, **kws)
				#ax.scatter(x, y, c='k', s=0.3, alpha=0.1)
				#im = ax.scatter(x, y, c=z, s=w, **kws)
			else:
				im = ax.imshow(im, **kws)
			ax.set_aspect('equal')
			#ax.plot(0,0,'ko')


			#spiral arms
			def spiral(ax, brange, bk, Rk, phi_lt, phi_gt, width=0, **kws):
				#brange=[180-230,180-12]
				def_kws = dict(color='k', linewidth=2.5, linestyle='-', alpha=0.3)
				def_kws.update(kws)
				b=np.deg2rad(np.arange(*brange,0.3))
				bk = np.deg2rad(bk)
				phi = np.deg2rad([phi_lt if v<=bk else phi_gt for v in b])
				for rk in [Rk,]:# Rk-width*1.1775, Rk+width*1.1775]:
					R = rk*np.exp(-(b-bk)*np.tan(phi))
					x=R*np.sin(b)
					y=R*np.cos(b)
					ax.plot(x, y, **def_kws)#'-' if rk==Rk else '--', 
					#ax.plot(Rk*np.sin(bk), Rk*np.cos(bk),'+', **kws)
				#from curvetext import CurvedText
				#CurvedText(x,y,label,va='bottom', axes=ax, fontsize=font['SMALL'], color=kws['color'])
			if arm:
				###Parameters in Reid2019
				#spiral(ax, [-150,190],15, 3.52, -4.2, -4.2, width=0.18, color='yellow')	#3kpc
				#spiral(ax, [-60,54],  18, 4.46, -1.0, 19.5, width=0.14, color='red')	#Norma
				#spiral(ax, [-330,81], 23, 4.91, 14.1, 12.1, width=0.23, color='green')	#Sct-Cen
				#spiral(ax, [-150,230],24, 6.04, 17.1,  1.0, width=0.27, color='magenta')#Sgt-Car
				#spiral(ax, [-30,34],   9, 8.26, 11.4, 11.4, width=0.31, color='cyan')	#local
				#spiral(ax, [-30,255], 40, 8.87, 10.3,  8.7, width=0.35, color='green')	#Perseus
				#spiral(ax, [-30,270], 18,12.24,  3.0,  9.4, width=0.65, color='green')	#Outer
				###Parameters adjusted
				#spiral(ax, [-140,18],15, 3.52, -4.2, -4.2, width=0.18, color='yellow')	#3kpc
				#spiral(ax, [45,200],90, 3.3, 4.2, 4.2, width=0.18, color='yellow')	#3kpc
				#spiral(ax, [-55,55],  18, 4.46, -0.5, 19.5, width=0.14, color='red')	#Norma
				#spiral(ax, [-54,81], 23, 4.91, 14.1, 12.1, width=0.23, color='blue')	#Sct-Cen
				#spiral(ax, [-340,-53],-54, 6.88, 10.0, 12.1, width=0.23, color='blue')
				#spiral(ax, [-33,230], 24, 6.04, 16.2,  7.1, width=0.27, color='magenta')#Sgt-Car
				#spiral(ax, [-170,-32],-33, 8.20, 9.9,  7.1, width=0.27, color='magenta')
				#spiral(ax, [-30,34],   9, 8.26, 11.4, 11.4, width=0.31, color='cyan')	#local
				#spiral(ax, [-30,255], 40, 8.87, 10.3, 13.0, width=0.35, color='black')	#Perseus
				#spiral(ax, [-30,305], 18,12.24,  4.9, 11.6, width=0.65, color='red')	#Outer
				###Parameters for Reid2019 figure1
				#spiral(ax, [-140,18],15, 3.52, -4.2, -4.2, width=0.18, label='3kpc', color='yellow')	#3kpc
				#spiral(ax, [45,200], 90,  3.3,  4.2,  4.2, width=0.18, color='yellow')	#3kpc
				#spiral(ax, [-20,55], 18, 4.46, -0.5, 19.5, width=0.14, label='Norma', color='darkred')	#Norma
				#spiral(ax, [-47,81], 23, 4.91, 14.1, 12.1, width=0.23, label='Scutum–Centaurus', color='darkblue')
				#spiral(ax, [-342,-192],-54, 6.88, 10.0, 12.1, width=0.23, label='Outer–Scutum–Centaurus', color='blue', linestyle='--')
				#spiral(ax, [-24,168],24, 6.04, 16.2,  7.1, width=0.27, label='Sagittarius', color='Purple')#Sgt-Car
				#spiral(ax, [-19,34],  9, 8.26, 11.4, 11.4, width=0.31, label='Local', color='darkcyan')	#local
				#spiral(ax, [-26,168],40, 8.87, 10.3, 13.0, width=0.35, label='Perseus', color='blue')	#Perseus
				#spiral(ax, [-27,168],18,12.24,  4.9, 11.6, width=0.65, label='Outer', color='blue')	#Outer
				### Parameters from Sun 2024
				spiral(ax, [-27,168],30.0, 10.1,  9.8,  9.8, color='k', linewidth=2.0, linestyle='-', label='CO')	#Outer
				spiral(ax, [-27,168],20.3, 13.3,  3.5, 11.1, color='k', linewidth=2.0, linestyle='-')	#Outer
				spiral(ax, [-27,168],47.0, 16.2, 12.3, 12.3, color='k', linewidth=2.0, linestyle='-')	#OSC
				### Reid 2019, brange clipped
				#spiral(ax, [-330,81], 23, 4.91, 14.1, 12.1, linewidth=0.23, color='green')	#Sct-Cen-OSC
				#spiral(ax, [-8,34], 9, 8.26, 11.4, 11.4, linewidth=0.23, color='green', label='Reid et al. 2019')	#Local
				#spiral(ax, [-23,115], 40, 8.87, 10.3,  8.7, linewidth=0.35, color='green')	#Perseus
				#spiral(ax, [-16, 71], 18,12.24,  3.0,  9.4, linewidth=0.65, color='green')	#Outer
				spiral(ax, [-342,-192],-54, 6.88, 10.0, 12.1, color='k', linewidth=2.0, linestyle=(0, (6,2)), label='HMSFRs') #OSC
				spiral(ax, [-19,34],     9, 8.26, 11.4, 11.4, color='k', linewidth=2.0, linestyle=(0, (6,2)))	#local
				spiral(ax, [-26,168],   40, 8.87, 10.3, 13.0, color='k', linewidth=2.0, linestyle=(0, (6,2)))	#Perseus
				spiral(ax, [-27,168],   18,12.24,  4.9, 11.6, color='k', linewidth=2.0, linestyle=(0, (6,2)))	#Outer
				### parameters from Drimmel 2025, or Poggio 2025
				spiral(ax, [-90,60], 0.0, np.exp(2.48), 20.2, 20.2, color='k', linewidth=3.0, linestyle=(0, (1, 1)), label='Classical Cepheids')	#Per (Drimmel et al. 2025)
				spiral(ax, [-90,60], 0.0, np.exp(2.21), 18.6, 18.6, color='k', linewidth=3.0, linestyle=(0, (1, 1)))	#Local
				ax.legend(loc=[0.02, 0.06], handletextpad=0.15, frameon=False, borderpad=0.1, labelspacing=0.1, fontsize=18)

			ax.set_autoscale_on(False)
			ax.set_xlim(-9, 16.5)
			ax.set_ylim(-15, 20)
			
			#Sun
			ax.plot(0, 8.15, color='k', marker='o', markerfacecolor='none', markersize=12, markeredgewidth=1.5, zorder=200)
			ax.plot(0, 8.15, color='k', marker='o', markersize=4, zorder=200)
			ax.text(0, 7.8, 'Sun', ha='center', va='top', fontsize=20, fontweight='normal', zorder=200)
			#GC
			ax.plot(0, 0, color='k', marker='*', markerfacecolor='k', markersize=18, zorder=200)
			ax.text(0, -0.5, 'GC', ha='center', va='top', fontsize=20, fontweight='normal', zorder=200)
			

			from matplotlib.patches import Ellipse
			ax.add_patch(Ellipse((0, 0), 1.5*2, 4.5*2, angle=-30.5, edgecolor='none', facecolor='gray', alpha=0.4))

			#grid
			if draw_phi:
				# draw dissigned phi
				from shared import rad_sep
				phi_sep = rad_sep
			else:
				# draw 10 deg separated phi
				phi_sep = range(0, 360, 10)
			if grid:
				for ray in phi_sep:
					ax.plot([0,np.sin(ray/180*np.pi)*30],[0,np.cos(ray/180*np.pi)*30], '--', color='gray', linewidth=1, alpha=0.3)
				for rad in range(2, 30, 2):
					ax.plot(rad*np.sin(np.linspace(0, 2*np.pi, 360)), rad*np.cos(np.linspace(0, 2*np.pi, 360)), '--', color='gray', linewidth=1, alpha=0.3)

			# rotation direction arrow
			if rotation_arrow:
				ax.plot(*RPHI2XY(np.full(100, 19.5), np.linspace(25, 50, 100)), 'k-', linewidth=1.5)
				ax.plot(*RPHI2XY([19.5, 19.9], [50, 47]), 'k-', linewidth=1.5)

			#colorbar
			#fig.colorbar(im, ax=ax, orientation='horizontal')
			ax.set_xlabel('X (kpc)')
			ax.set_ylabel('Y (kpc)')

			'''
			axins = inset_axes(ax,
				   width="4%",  # width = 10% of parent_bbox width
				   height="35%",  # height : 40%
				   loc='lower left',
				   bbox_to_anchor=(0.1, 0.1, 1, 1),
				   bbox_transform=ax.transAxes,
				   borderpad=0,
				   )
			'''
			if im:
				vis = colorbar_kws.pop('visible', True)
				if vis:
					ori = colorbar_kws.pop('orientation', 'vertical')
					cax = fig.add_axes(**colorbar_kws)
					cbar = fig.colorbar(im, cax=cax, orientation=ori)
					if ori=='vertical':
						cbar.ax.set_ylabel('$\Delta$Z (kpc)', fontsize=15, labelpad=0)
						cax.tick_params(axis='y', which='major', labelsize=15)
					else:
						cbar.ax.set_xlabel('$\Delta$Z (kpc)', fontsize=15)
						cax.tick_params(axis='x', which='major', labelsize=15)
			#cbar.ax.set_ylabel('Z') # 可选：为colorbar添加标签
		

		plt.rcParams['savefig.dpi'] = 280
		plt.rcParams['axes.linewidth'] = 0.8
		plt.rcParams['axes.labelsize'] = 20
		plt.rcParams['axes.labelweight'] = 'bold'
		plt.rcParams['xtick.labelsize'] = 18
		plt.rcParams['ytick.labelsize'] = 18
		plt.rcParams['xtick.direction'] = 'in'
		plt.rcParams['ytick.direction'] = 'in'
		plt.rcParams['xtick.top'] = True
		plt.rcParams['ytick.right'] = True
		plt.rcParams['xtick.minor.visible'] = True
		plt.rcParams['ytick.minor.visible'] = True
		plt.rcParams['legend.fontsize'] = 15
		if 0:
			figscale = 0.53
			figwidth = textwidth*figscale
			fontscale = 1.3
			# single plot without arm
			from shared import function_PoggioWarp

			fig, ax = plt.subplots(ncols=1, nrows=1, sharex=True, sharey=True, figsize=(figwidth, figwidth*1.33))
			plt.subplots_adjust(left=0.11, right=0.97, top=0.95, wspace=0.1, hspace=0.1)

			# image
			from shared import function_warp, p_1comp, p_2comp, subfigureIndexFont, p_circ1comp
			#Z_w = function_warp((R,PHI), p=p_1comp if component==1 else p_2comp)
			#Z -= function_warp((R,PHI), p=[0,0,0,0,0], circ=p_circ1comp)
			kws = dict(cmap='coolwarm', vmin=-0.3, vmax=0.3, colorbar_kws={'rect':[0.2, 0.14, 0.025, 0.22], 'orientation':'vertical'})
			gplane(fig, ax, im=None, arm=False, grid=True, rotation_arrow=True, **kws)

			#ring slice
			kws = dict(color='k', linewidth=1, linestyle='-')
			r = 12.7
			phi = np.linspace(-43, -21.5, 100)
			ax.plot(*RPHI2XY(r, phi), **kws)
			#ax.text(*RPHI2XY(11.8, -33), '12.5 kpc', rotation=33, ha='center', va='center', fontsize=14, fontweight='bold',color='red')
			ax.text(*RPHI2XY(r-0.7, -38), '%.1f' % r, rotation=38, ha='center', va='center', fontsize=20, fontweight='bold')
			ax.text(*RPHI2XY(r-0.7, -27), 'kpc', rotation=27, ha='center', va='center', fontsize=20, fontweight='bold')

			# phi names
			from shared import rad_sep
			theta_sep = rad_sep
			sep = -1
			labelR = [15.1]*8+[16.1, 17.2]+[18.3]*5
			for i in range(len(theta_sep)-2):
				### skip 11 ~ -11
				if i == 13: sep += 2
				else: sep += 1
				theta1, theta2 = theta_sep[sep:sep+2]
				theta = (theta1+theta2)/2
				ax.text(*RPHI2XY(labelR[i], theta), r'$\mathbf{\phi_{%i}}$' % (i+1), ha='center', va='center', rotation=0, fontsize=20, fontweight='bold')

			# phi ticks
			i, o = 3.2, 4.5
			ax.plot([-i, -o], [0, 0], 'k-')
			ax.text(0, -o-0.2, '$180^\circ$', va='top', ha='center', fontsize=16, fontweight='bold')
			ax.plot([+i, +o], [0, 0], 'k-')
			ax.text(0, +o+0.1, '$\phi=0^\circ$', va='bottom', ha='center', fontsize=16, fontweight='bold')
			ax.plot([0, 0], [-i, -o], 'k-')
			ax.text(-o-0.1, 0, '$-90^\circ$', va='center', ha='right', fontsize=16, fontweight='bold')
			ax.plot([0, 0], [+i, +o], 'k-')
			ax.text(+o+0.2, 0, '$90^\circ$', va='center', ha='left', fontsize=16, fontweight='bold')

			# subfigure index
			#ax.text(-0.1, 1, 'a', font=subfigureIndexFont, transform = ax.transAxes, ha='right', va='top')
			fig.savefig('fig/dZ_plane_%1icomp.%s' % (component, 'png'), bbox_inches='tight')
			fig.savefig('fig/dZ_plane_%1icomp.%s' % (component, 'pdf'), bbox_inches='tight')
		elif 1:
			figscale = 0.39 if component==1 else 0.63
			figwidth = textwidth*figscale
			fontscale = 1.3
			# single plot
			from shared import function_PoggioWarp

			fig, ax = plt.subplots(ncols=1, nrows=1, sharex=True, sharey=True, figsize=(figwidth, figwidth*1.33))
			plt.subplots_adjust(left=0.11, right=0.97, top=0.95, wspace=0.1, hspace=0.1)
			from shared import function_warp, p_1comp, p_2comp, subfigureIndexFont, p_circ1comp
			#Z_w = function_warp((R,PHI), p=p_1comp if component==1 else p_2comp)
			#Z -= function_warp((R,PHI), p=[0,0,0,0,0], circ=p_circ1comp)
			#Z += Z_w
			#Z_w = function_PoggioWarp((R,PHI))
			#Z -= Z_w
			kws = dict(cmap='coolwarm', vmin=-0.3, vmax=0.3, colorbar_kws={'rect':[0.18, 0.27, 0.025, 0.22], 'orientation':'vertical'})
			gplane(fig, ax, im=None, arm=True, grid=False, **kws)

			from matplotlib.patches import FancyArrowPatch
			ax.add_patch(FancyArrowPatch(RPHI2XY(15, 38), RPHI2XY(13.8, 38), arrowstyle='-|>', mutation_scale=12, color='k', linewidth=3))
			ax.add_patch(FancyArrowPatch(RPHI2XY(16, 4), RPHI2XY(14.8, 4), arrowstyle='-|>', mutation_scale=12, color='k', linewidth=3))
			ax.add_patch(FancyArrowPatch(RPHI2XY(15, -20), RPHI2XY(13.8, -20), arrowstyle='-|>', mutation_scale=12, color='k', linewidth=3))

			#quadrant ticks
			l = 1.8
			xlim = ax.get_xlim()
			ax.plot([xlim[0], xlim[0]+l], [8.15, 8.15], '-', color='#888888')
			ax.text(xlim[0]+0.1, 8.15-0.1, '$270^\circ$', va='top', ha='left', fontsize=12)
			ax.text(xlim[0]+l, 8.15+1, 'Q3', va='bottom', ha='center', fontsize=15)
			ax.text(xlim[0]+l, 8.15-1, 'Q4', va='top', ha='center', fontsize=15)

			ax.plot([xlim[1]-l, xlim[1]], [8.15, 8.15], '-', color='#888888')
			ax.text(xlim[1]-0.1, 8.15-0.1, '$90^\circ$', va='top', ha='right', fontsize=12)
			ax.text(xlim[1]-l, 8.15+1, 'Q2', va='bottom', ha='center', fontsize=15)
			ax.text(xlim[1]-l, 8.15-1, 'Q1', va='top', ha='center', fontsize=15)

			ylim = ax.get_ylim()
			ax.plot([0, 0], [ylim[0], ylim[0]+l], '-', color='#888888')
			ax.text(-0.1, ylim[0]+0.1, 'l$=0^\circ$', rotation=90, va='bottom', ha='right', fontsize=12, fontfamily='Georgia', fontstyle='italic')

			ax.plot([0, 0], [ylim[1]-l, ylim[1]], '-', color='#888888')
			ax.text(-0.1, ylim[1]-0.1, '$180^\circ$', rotation=90, va='top', ha='right', fontsize=12)

			#ax.text(-0.1, 1, 'a', font=subfigureIndexFont, transform = ax.transAxes, ha='right', va='top')
			fig.savefig('fig/dZ_plane_%1icomp_arm.%s' % (component, 'png'), bbox_inches='tight')
			fig.savefig('fig/dZ_plane_%1icomp_arm.%s' % (component, 'pdf'), bbox_inches='tight')
		elif 0:
			figscale = 0.53
			figwidth = textwidth*figscale
			fontscale = 1.3
			# single plot without arm
			from shared import function_PoggioWarp

			fig, ax = plt.subplots(ncols=1, nrows=1, sharex=True, sharey=True, figsize=(figwidth, figwidth*1.33))
			plt.subplots_adjust(left=0.11, right=0.97, top=0.95, wspace=0.1, hspace=0.1)
			from shared import function_warp, p_1comp, p_2comp, subfigureIndexFont, p_circ1comp
			
			kws = dict(cmap='coolwarm', vmin=-0.3, vmax=0.3, colorbar_kws={'rect':[0.2, 0.14, 0.025, 0.22], 'orientation':'vertical'})
			gplane(fig, ax, im=None, arm=False, grid=True, rotation_arrow=True, **kws)

			step = 0.05
			#Z_w = function_warp((R, PHI), p=[0,0,0,0,0], circ=p_circ1comp)
			#gridX, gridY, gridZ = convolveZonGrid(X, Y, Z_w, mass, step=(step, step), polar=False, useMask=True, generateMask=False)
			#extent=(gridX[0,0]-step/2, gridX[0,-1]+step/2, gridY[0,0]-step/2, gridY[-1,0]+step/2)
			#gplane(fig, ax, im=gridZ, arm=False, grid=True, extent=extent, origin='lower', rotation_arrow=True, **kws)
			#ax.contour(gridX, gridY, gridZ, levels=[-0.1, 0.1], colors=['blue', 'red'], linestyles=':')

			# phi names
			from shared import rad_sep
			theta_sep = rad_sep
			sep = -1
			labelR = [15.1]*8+[16.1, 17.2]+[18.3]*5
			for i in range(len(theta_sep)-2):
				### skip 11 ~ -11
				if i == 13: sep += 2
				else: sep += 1
				theta1, theta2 = theta_sep[sep:sep+2]
				theta = (theta1+theta2)/2
				ax.text(*RPHI2XY(labelR[i], theta), r'$\mathbf{\phi_{%i}}$' % (i+1), ha='center', va='center', rotation=0, fontsize=20, fontweight='bold')

			kws = dict(color='k', linewidth=4, linestyle=':', alpha=1)
			r = 12.3
			phi = np.linspace(-44, 52, 100)
			ax.plot(*RPHI2XY(r, phi), **kws)
			#ax.text(*RPHI2XY(11.8, -33), '12.5 kpc', rotation=33, ha='center', va='center', fontsize=14, fontweight='bold',color='red')
			ax.text(*RPHI2XY(11.6, -38), '12.3', rotation=38, ha='center', va='center', fontsize=20, fontweight='bold')
			ax.text(*RPHI2XY(11.6, -27), 'kpc', rotation=27, ha='center', va='center', fontsize=20, fontweight='bold')

			from matplotlib.patches import FancyArrowPatch
			ax.add_patch(FancyArrowPatch(RPHI2XY(15., 38), RPHI2XY(13.8, 38), arrowstyle='-|>', mutation_scale=12, color='k', linewidth=3))
			ax.add_patch(FancyArrowPatch(RPHI2XY(16., 4), RPHI2XY(14.8, 4), arrowstyle='-|>', mutation_scale=12, color='k', linewidth=3))
			ax.add_patch(FancyArrowPatch(RPHI2XY(15., -20), RPHI2XY(13.8, -20), arrowstyle='-|>', mutation_scale=12, color='k', linewidth=3))

			# subfigure index
			#ax.text(-0.1, 1, 'a', font=subfigureIndexFont, transform = ax.transAxes, ha='right', va='top')
			fig.savefig('fig/dZ_plane_%1icomp_greatwave.%s' % (component, 'pdf'), bbox_inches='tight')
			#fig.savefig('fig/dZ_plane_%1icomp.%s' % (component, 'pdf'), bbox_inches='tight')
		elif 0:
			# compare
			figscale = 0.53
			figwidth = textwidth*figscale
			fontscale = 1.3
			fig, ax = plt.subplots(ncols=1, nrows=3, sharex=True, sharey=True, figsize=(figwidth, figwidth*1.1))# figwidth*0.91))
			plt.subplots_adjust(left=0.09, right=0.88, bottom=0.08, top=0.95, wspace=0.04, hspace=0.02)

			plt.rcParams['savefig.format'] = 'pdf'
			plt.rcParams['savefig.dpi'] = 280
			plt.rcParams['axes.linewidth'] = 0.8
			plt.rcParams['axes.labelsize'] = 14
			plt.rcParams['axes.labelweight'] = 'bold'
			plt.rcParams['xtick.labelsize'] = 12
			plt.rcParams['ytick.labelsize'] = 12
			plt.rcParams['xtick.direction'] = 'in'
			plt.rcParams['ytick.direction'] = 'in'
			plt.rcParams['xtick.top'] = True
			plt.rcParams['ytick.right'] = True
			plt.rcParams['xtick.minor.visible'] = True
			plt.rcParams['ytick.minor.visible'] = True
			plt.rcParams['legend.fontsize'] = 12

			kws = dict(cmap='coolwarm', vmin=-0.3, vmax=0.3)

			from shared import function_warp, p_1comp, p_2comp, subfigureIndexFont, function_PoggioWarp

			#ax[0].set_title('MWISP Clouds - CO warp model')
			step = 0.25
			gridX, gridY, gridZ = convolveZonGrid(X, Y, Z, mass, step=(step, step), polar=False, useMask=True, generateMask=False)
			extent=(gridX[0,0]-step/2, gridX[0,-1]+step/2, gridY[0,0]-step/2, gridY[-1,0]+step/2)
			gplane(fig, ax[0], im = gridZ, extent=extent, origin='lower', arm=False, draw_phi=False, colorbar_kws={'rect':[0.20, 0.68, 0.015, 0.2]}, **kws)
			ax[0].set_xlabel('')
			ax[0].set_ylabel('')
			#ax[0].plot(0, 8.15, color='k', marker='o', markerfacecolor='none', markersize=12, zorder=200)
			#ax[0].plot(0, 8.15, color='k', marker='o', markersize=4, zorder=200)
			#ax[0].text(0, 7.2,, 'Sun', ha='center', va='top', fontsize=20, fontweight='normal', zorder=200)
			ax[0].set_yticks(np.arange(0,20,5))

			#ax[1].set_title('Young giants (Poggio et al. 2025)')
			imgPoggio = np.load('PoggioFigure/Young_giants_fitstraightLON_DZ_residuals_fit.npy')
			extent = [-15., 15., -15.+8.15, 15.+8.15]
			gplane(fig, ax[1], im = imgPoggio, extent=extent, arm=False, draw_phi=False, colorbar_kws={'visible':False}, **kws)
			ax[1].set_xlabel('')
			ax[1].set_ylabel('')
			#ax[1].tick_params(labelleft=False)
			#ax[1].plot(8.15, 0, color='k', marker='o', markerfacecolor='none', markersize=12, zorder=200)
			#ax[1].plot(8.15, 0, color='k', marker='o', markersize=4, zorder=200)
			#ax[1].text(7.2, -0.6, 'Sun', ha='center', va='top', fontsize=20, fontweight='normal', zorder=200)
			ax[1].set_yticks(np.arange(0,20,5))

			#ax[2].set_title('Cepheids (Poggio et al. 2025)')
			imgPoggio = np.load('PoggioFigure/Cepheids_fitstraightLON_DZ_residuals_fit.npy')
			extent = [-15., 15., -15.+8.15, 15.+8.15]
			gplane(fig, ax[2], im = imgPoggio, extent=extent, arm=False, draw_phi=False, colorbar_kws={'visible':False}, **kws)
			#ax[2].plot(8.15, 0, color='k', marker='o', markerfacecolor='none', markersize=12, zorder=200)
			#ax[2].plot(8.15, 0, color='k', marker='o', markersize=4, zorder=200)
			#ax[2].text(7.2, -0.6, 'Sun', ha='center', va='top', fontsize=20, fontweight='normal', zorder=200)
			ax[2].set_yticks(np.arange(0,20,5))
			ax[2].set_xlabel('X (kpc)')
			ax[2].set_ylabel('Y (kpc)')

			
			from matplotlib.patches import FancyArrowPatch
			for a in ax.ravel():
				a.add_patch(FancyArrowPatch(RPHI2XY(15., 38), RPHI2XY(13.8, 38), arrowstyle='-|>', mutation_scale=12, color='k', linewidth=3))
				a.add_patch(FancyArrowPatch(RPHI2XY(15.5, 4), RPHI2XY(14.3, 4), arrowstyle='-|>', mutation_scale=12, color='k', linewidth=3))
				a.add_patch(FancyArrowPatch(RPHI2XY(15., -20), RPHI2XY(13.8, -20), arrowstyle='-|>', mutation_scale=12, color='k', linewidth=3))			
			ax[0].set_xlim(-12, 12)
			ax[0].set_ylim(4.5, 15.5)
			#ax[2].set_ylabel('')
			#ax[1].set_xlim(-15, 15)
			#ax[1].set_ylim(-15+8.15, 19)
			fig.savefig('fig/compare_dZ_plane_%1icomp2v.%s' % (component, 'pdf'), bbox_inches='tight')
			#ax.text(-0.1, 1, 'a', font=subfigureIndexFont, transform = ax.transAxes, ha='right', va='top')			
		else:
			# compare
			figscale = 0.53
			figwidth = textwidth*figscale
			fontscale = 1.3
			fig, ax = plt.subplots(ncols=2, nrows=2, sharex=True, sharey=True, figsize=(figwidth, figwidth*1.0))# figwidth*0.91))
			plt.subplots_adjust(left=0.09, right=0.88, bottom=0.08, top=0.95, wspace=0.04, hspace=0.02)

			plt.rcParams['savefig.format'] = 'pdf'
			plt.rcParams['savefig.dpi'] = 280
			plt.rcParams['axes.linewidth'] = 0.8
			plt.rcParams['axes.labelsize'] = 14
			plt.rcParams['axes.labelweight'] = 'bold'
			plt.rcParams['xtick.labelsize'] = 12
			plt.rcParams['ytick.labelsize'] = 12
			plt.rcParams['xtick.direction'] = 'in'
			plt.rcParams['ytick.direction'] = 'in'
			plt.rcParams['xtick.top'] = True
			plt.rcParams['ytick.right'] = True
			plt.rcParams['xtick.minor.visible'] = True
			plt.rcParams['ytick.minor.visible'] = True
			plt.rcParams['legend.fontsize'] = 12

			kws = dict(cmap='coolwarm', vmin=-0.3, vmax=0.3)

			from shared import function_warp, p_1comp, p_2comp, subfigureIndexFont, function_PoggioWarp

			ax[0,0].set_title('MWISP Clouds - CO warp model')
			step = 0.25
			Z_w = function_warp((R,PHI), p=p_1comp if component==1 else p_2comp)
			gridX, gridY, gridZ = convolveZonGrid(X, Y, Z, mass, step=(step, step), polar=False, useMask=True, generateMask=False)
			extent=(gridX[0,0]-step/2, gridX[0,-1]+step/2, gridY[0,0]-step/2, gridY[-1,0]+step/2)
			gplane(fig, ax[0,0], im = gridZ, extent=extent, origin='lower', arm=True, draw_phi=False, colorbar_kws={'visible':False}, **kws)
			ax[0,0].set_xlabel('')
			ax[0,0].set_ylabel('')

			ax[0,1].set_title('MWISP Clouds - Cepheids warp model')
			Z_w = function_PoggioWarp((R,PHI), cepheids=True, straight=True)
			gridX, gridY, gridZ = convolveZonGrid(X, Y, Z, mass, step=(step, step), polar=False, useMask=True, generateMask=False)
			extent=(gridX[0,0]-step/2, gridX[0,-1]+step/2, gridY[0,0]-step/2, gridY[-1,0]+step/2)
			gplane(fig, ax[0,1], im = gridZ, extent=extent, origin='lower', arm=False, draw_phi=False, colorbar_kws={'visible':False}, **kws)
			ax[0,1].set_xlabel('')
			ax[0,1].set_ylabel('')
			ax[0,1].tick_params(labelleft=False)

			ax[1,0].set_title('Cepheids (Poggio et al. 2025)')
			imgPoggio = np.load('PoggioFigure/Cepheids_fitstraightLON_DZ_residuals_fit.npy')
			extent = [-15., 15., -15.+8.15, 15.+8.15]
			gplane(fig, ax[1,0], im = imgPoggio, extent=extent, arm=False, draw_phi=False, colorbar_kws={'visible':False}, **kws)

			ax[1,1].set_title('Young giants (Poggio et al. 2025)')
			imgPoggio = np.load('PoggioFigure/Young_giants_fitstraightLON_DZ_residuals_fit.npy')
			extent = [-15., 15., -15.+8.15, 15.+8.15]
			gplane(fig, ax[1,1], im = imgPoggio, extent=extent, arm=False, draw_phi=False, colorbar_kws={'rect':[0.89, 0.138, 0.015, 0.304]}, **kws)
			ax[1,1].set_xlabel('')
			ax[1,1].set_ylabel('')
			ax[1,1].tick_params(labelleft=False)


			from matplotlib.patches import FancyArrowPatch
			for a in ax.ravel():
				a.add_patch(FancyArrowPatch(RPHI2XY(15.5, 38), RPHI2XY(13.8, 38), arrowstyle='-|>', mutation_scale=12, color='k', linewidth=3))
				a.add_patch(FancyArrowPatch(RPHI2XY(16.5, 4), RPHI2XY(14.8, 4), arrowstyle='-|>', mutation_scale=12, color='k', linewidth=3))
				a.add_patch(FancyArrowPatch(RPHI2XY(15.5, -20), RPHI2XY(13.8, -20), arrowstyle='-|>', mutation_scale=12, color='k', linewidth=3))
			ax[1,1].set_xlim(-16, 15)
			ax[0,0].set_ylim(-13, 18)
			ax[0,1].set_ylim(5, 18)
			#ax[1,0].set_ylim(-5, 18)
			#ax[1,1].set_ylim(-5, 18)
			#ax[2].set_ylabel('')
			#ax[1].set_xlim(-15, 15)
			#ax[1].set_ylim(-15+8.15, 19)
			fig.savefig('fig/compare_dZ_plane_%1icomp.%s' % (component, 'pdf'), bbox_inches='tight')
			#ax.text(-0.1, 1, 'a', font=subfigureIndexFont, transform = ax.transAxes, ha='right', va='top')



	### 3D visualization with plotly (finally adopted)
	if 0 and not sin:
		import plotly.graph_objs as go
		figscale = 0.5
		fontscale = 1.3

		### Data
		#points = go.Scatter3d(x=X, y=Y, z=Z - function((R,PHI), bestmed, sin=False), mode='markers', \
		#	marker=dict(size=3, color=Z, colorscale='RdYlBu_r', cmin=-1, cmax=1, opacity=0.2, colorbar=dict(thickness=20)))
		#points = go.Scatter3d(x=X, y=Y, z=Z, mode='markers', \
		#	marker=dict(size=10, color=np.log10(mass), colorscale='gray', opacity=0.3))#, colorbar=dict(thickness=20)

		
		### median
		gridX, gridY, gridZ = convolveZonGrid(X, Y, Z, mass, step=(0.75, 5), polar=True, useMask=False, generateMask=False)
		gridR, gridPHI = XY2RPHI(gridX, gridY)

		#gridZ = gridZ - function((gridR, gridPHI), bestmed, warp1=True, sin=False)
		#gridR, gridPHI, gridZ = convolveZonGrid(X, Y, Z, mass, step=(0.5, 0.5), polar=False)
		points = go.Scatter3d(x=gridX.ravel(), y=gridY.ravel(), z=gridZ.ravel(), mode='markers', \
			marker=dict(size=8, color=gridZ.ravel(), colorscale='RdYlBu_r', cmin=-0.5, cmax=0.5, opacity=1.0, \
			colorbar=dict(thickness=30, title=dict(text='Z (kpc)', font=dict(size=int(30*fontscale), weight='bold')), \
				x=0.2, y=0.08, len=0.3, orientation='h', tickvals=np.arange(-0.5,1,0.25), tickfont=dict(size=int(25*fontscale)))))


		### wireframe
		wires = []
		line_marker = dict(color='#000000', width=5)#, colorscale='RdYlBu_r', cmin=-0.5, cmax=0.5)
		### R axis
		from shared import rad_sep, to_subscript
		for phi in list(range(323, 155, -12))+rad_sep:
			#for phi in np.arange(0, 360, 15):
			lineR = np.arange(8.5, 20.6, 0.5)
			linePHI = np.repeat(phi, lineR.size)
			lineX, lineY = RPHI2XY(lineR, linePHI)
			lineZ = function((lineR, linePHI), bestmed, warp=True, sin=False)
			#line_marker['color'] = lineZ
			wires.append(go.Scatter3d(x=lineX, y=lineY, z=lineZ, mode='lines', line=line_marker))

		### phi axis 
		for r in [8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 17.5, 20.5]:
			linePHI = np.linspace(0, 360, int(360*r)//40)
			lineR = np.repeat(r, linePHI.size)
			lineX, lineY = RPHI2XY(lineR, linePHI)
			lineZ = function((lineR, linePHI), bestmed, warp=True, sin=False)
			#line_marker['color'] = lineZ
			wires.append(go.Scatter3d(x=lineX, y=lineY, z=lineZ, mode='lines', line=line_marker))

		### phi sector text
		sectorPhi = np.array([(p1+p2)/2 for p1,p2 in zip(rad_sep[:-1], rad_sep[1:]) if p1*p2>0])#avoid az=0
		sectorR = np.repeat(22.5, len(sectorPhi))
		sectorZ = function((sectorR, sectorPhi), bestmed)
		sectorZ[sectorZ>2.21] = 2.21
		sectorX, sectorY = RPHI2XY(sectorR, sectorPhi)
		#sectorT = ['ϕ'+to_subscript(i+1) for i in range(len(sectorPhi))] # different format
		sectorT = ['$\huge \phi_{%i}$' % (i+1) for i in range(len(sectorPhi))]
		sectorText = [dict(x=sectorX[i], y=sectorY[i], z=sectorZ[i], text=sectorT[i], showarrow=False, \
				xanchor = 'center', xshift=5, yanchor='middle', font=dict(color='#444444', size=int(30*fontscale))) for i in range(len(sectorT))]

		### radius text
		annulusR = np.array([(p1+p2)/2 for p1,p2 in zip([8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 17.5], [9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 17.5, 20.5])])
		annulusPhi = np.repeat(-40, len(annulusR))
		annulusZ = function((annulusR, annulusPhi), bestmed)
		annulusX, annulusY = RPHI2XY(annulusR, annulusPhi)
		#annulusT = np.array(['R'+to_subscript(i+1) for i in range(len(annulusR))]) # different format
		annulusT = ['$\huge R_{%i}$' % (i+1) for i in range(len(annulusR))]
		annulusText = [dict(x=annulusX[i], y=annulusY[i], z=annulusZ[i], text='$\huge R_{%i}$' % (i+1), showarrow=False, \
			xanchor = 'center', xshift=5, yanchor='middle', font=dict(color='#444444', size=int(20*fontscale))) for i in [0,8]]
		
		sun = go.Scatter3d(x=[0], y=[8.15], z=[0], mode='markers', \
				marker=dict(color='#ff6600', size=6))
		sunText = dict(x=0, y=8.15, z=0, text='Sun', showarrow=False,
			xanchor='left', xshift=5, yanchor='bottom', font=dict(color='#ff6600', size=int(28*fontscale)))

		gc = go.Scatter3d(x=[0], y=[0], z=[0], mode='markers', \
			marker=dict(color='#555555', size=6))
		gcText = dict(x=0, y=0, z=0, text='GC', showarrow=False,
			xanchor='left', xshift=5, yanchor='bottom', font=dict(color='#666666', size=int(28*fontscale)))

		eyeD = 1.6
		eyeAz = -35
		eyeEl = 26
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
				bestmedi = np.mean(steps[probs > np.percentile(probs, 95)], axis=0)
				_r.append(i+0.5)
				_phi.append(bestmedi[-1])
				_z.append(0)
			print(_r, _phi)
			_x, _y = RPHI2XY(np.array(_r), np.array(_phi))
			_rphiw = go.Scatter3d(x=_x, y=_y, z=_z, mode='lines+markers',\
				marker=dict(symbol='square', size=6, color='darkgreen', opacity=1), line=dict(color='darkgreen', width=12))
			_r = [8.5, 20.5]
			_phi = [bestmed[-1]]*2
			_x, _y = RPHI2XY(np.array(_r), np.array(_phi))
			_phiw_1comp = go.Scatter3d(x=_x, y=_y, z=[0, 0], mode='lines',\
				marker=dict(symbol='square', size=0.1, color='darkgreen', opacity=1), line=dict(dash='longdash', color='darkgreen', width=10))

		### layout format
		title_font = dict(size=int(28*fontscale), weight='bold')
		tick_kws = dict(showticklabels=True, \
			ticks='outside', tickcolor='#ffffff', tickwidth=8, ticklen=20, tickfont=dict(size=int(18*fontscale)), \
			gridcolor='#dddddd', gridwidth=5,
			showspikes = False, showbackground=False, backgroundcolor='white')

		layout = go.Layout(
			scene=dict(
				xaxis=dict(title = '',#dict(text='X (kpc)', font=title_font),
					range = [-24.01, 24.01],
					tickvals = np.arange(-24, 25, 8),
					#ticktext = ['%+i' % x for x in range(-24, 25, 8)],
					**tick_kws),
				yaxis=dict(title = '',#dict(text='Y (kpc)', font=title_font),
					range = [-24.01, 24.01],
					tickvals = np.arange(-24, 25, 8), 
					#ticktext = ["%+i" % y for y in range(-24, 25, 8)],
					**tick_kws),
				zaxis=dict(title = '',#dict(text='Z (kpc)', font=title_font),
					range = [-2, 2 if component==1 else 2.21],
					tickvals = np.arange(-1, 3, 1),
					#ticktext = ["%+i" % z for z in range(-1, 3, 1)],
					zerolinecolor='#000000', zerolinewidth=5, **tick_kws),
				annotations=[sunText, gcText, *sectorText, *annulusText],
				aspectmode='manual',
				aspectratio=dict(x=1, y=1, z=.45),
				camera = dict(
					up=dict(x=0, y=0, z=1),
					center=dict(x=0.11, y=0, z=-0.3),
					eye=eyePos
					),
				),
			showlegend=False,
			width=int(18*400*figscale/2),	#textwidth * dpi * texfigurescale
			height=int(18*400*figscale*0.8/2),
			margin=dict(r=0, l=0, b=0, t=0, pad=20),
			#scene_aspectmode='manual',
			#scene_aspectratio=dict(x=1, y=1, z=.45),
		)


		# Create a Figure object and add the 3d objects to it
		if component==1:
			fig = go.Figure(data=[points, *wires, sun, gc, _phiw_1comp], layout=layout)#_rphiw
		else:
			fig = go.Figure(data=[points, *wires, sun, gc], layout=layout)

		### axis title
		fig.add_annotation(x=0.24, y=0.22, textangle=34, text='X (kpc)', xref='paper', yref='paper', \
			showarrow=False, font=dict(color='rgb(37,63,98)', size=int(40*fontscale), weight='bold'))
		fig.add_annotation(x=0.93, y=0.27, textangle=-55, text='Y (kpc)', xref='paper', yref='paper', \
			showarrow=False, font=dict(color='rgb(37,63,98)', size=int(40*fontscale), weight='bold'))
		fig.add_annotation(x=0.02, y=0.61, textangle=74, text='Z (kpc)', xref='paper', yref='paper', \
			showarrow=False, font=dict(color='rgb(37,63,98)', size=int(40*fontscale), weight='bold'))

		### panel ID
		if component==1:
			from shared import subfigureIndexFont
			#fig.add_annotation(x=0.05, xref='paper', y=0.95, yref='paper', text='a', showarrow=False, font=dict(color='black', family="Arial Black", size=int(35/figscale)))

		### Show the plot
		if 0: fig.show()
		else:
			fig.write_image("fig/median_model_%icomp.pdf" % component)#, scale=3)
			#fig.write_html("fig/median_model_%icomp.html" % component)#, scale=3)



	### 3D visualization of residual with plotly (finally adopted)
	if 0 and sin:
		figscale = 0.62
		fontscale = 1.3

		from shared import *
		import plotly.graph_objs as go
		
		
		### residual
		resX, resY, resZ = convolveZonGrid(X, Y, Z, mass, step=(0.25, 1), polar=True)

		##resX, resY, resZ = _excludedXYZ(resX, resY, resZ)
		residual = go.Surface(x=resX, y=resY, z=resZ, surfacecolor=resZ, colorscale='RdYlBu_r', cmin=-0.2, cmax=0.2, opacity=1,
			colorbar = dict(thickness=30, title=dict(text='ΔZ (kpc)',font=dict(size=int(60*figscale), weight='bold')), \
				x=0.72, y=0.12, len=0.3, orientation='h', tickvals=np.arange(-0.2,0.5,0.1), tickfont=dict(size=int(54*figscale))))
		
		#residual = go.Scatter3d(x=resX.ravel(), y=resY.ravel(), z=resZ.ravel(), mode='markers', \
		#	marker=dict(size=5, color=resZ.ravel(), colorscale='RdYlBu_r', cmin=-0.2, cmax=0.2, opacity=1.0, \
		#	colorbar=dict(thickness=20, title='Z (kpc)', x=0.24, y=0.7, len=0.05, len=0.3, orientation='h')))
		
		
		### residual wireframe
		wires = []
		line_marker = dict(color='#223322', width=5)
		#for phi in np.arange(0, 360, 15):
		from shared import rad_sep
		for phi in list(range(323, 155, -12))+rad_sep:
			lineR = np.arange(0, 22.1, 0.25)
			linePHI = np.repeat(phi, lineR.size)
			lineX, lineY = RPHI2XY(lineR, linePHI)
			lineZ = convolveZ(X, Y, Z, mass, lineX, lineY)
			if phi in [29,38]: line_marker['width']=15
			else: line_marker['width']=5
			wires.append(go.Scatter3d(x=lineX, y=lineY, z=lineZ+0.01, mode='lines', line=line_marker))

		for r in np.arange(8, 22.1, 1):
			#linePHI = np.linspace(0, 360, int(360*r)//40)
			linePHI = np.arange(0, 360.1, 1)
			lineR = np.repeat(r, linePHI.size)
			lineX, lineY = RPHI2XY(lineR, linePHI)
			lineZ = convolveZ(X, Y, Z, mass, lineX, lineY)
			if r in [12,]: line_marker['width']=15
			else: line_marker['width']=5
			wires.append(go.Scatter3d(x=lineX, y=lineY, z=lineZ+0.01, mode='lines', line=line_marker))
		


		### phi sector text
		sectorPhi = np.array([(p1+p2)/2 for p1,p2 in zip(rad_sep[:-1], rad_sep[1:]) if p1*p2>0])
		sectorR = np.linspace(17, 19, sectorPhi.size)
		sectorZ = np.zeros(sectorPhi.size)
		sectorX, sectorY = RPHI2XY(sectorR, sectorPhi)
		#sectorT = ['ϕ'+to_subscript(i+1) for i in range(len(sectorPhi))]
		sectorT = ['$\huge \phi_{%i}$' % (i+1) for i in range(len(sectorPhi))]
		#for i in range(len(sectorPhi)):
		sectorText = [dict(x=sectorX[i], y=sectorY[i], z=sectorZ[i], text=sectorT[i], showarrow=False, \
				xanchor = 'center', xshift=5, yanchor='middle', font=dict(color='#000000', size=int(60*fontscale))) for i in range(len(sectorPhi))] 
		

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
			xanchor='left', xshift=5, yanchor='bottom', font=dict(color='#666666', size=int(30*fontscale)))

		sun = go.Scatter3d(x=[0], y=[8.15], z=[0], mode='markers', \
				marker=dict(color='#ff6600', size=8))
		sunText = dict(x=0, y=8.15, z=0, text='Sun', showarrow=False,
			xanchor='left', xshift=5, yanchor='bottom', font=dict(color='#ff6600', size=int(30*fontscale)))

		gc = go.Scatter3d(x=[0], y=[0], z=[0], mode='markers', \
			marker=dict(color='#555555', size=8))
		gcText = dict(x=0, y=0, z=0, text='GC', showarrow=False,
			xanchor='left', xshift=5, yanchor='bottom', font=dict(color='#666666', size=int(30*fontscale)))

		### arrow
		arrowR = 21
		a = np.linspace(24.5, 53, 200) /180*np.pi
		arrowX = np.concatenate( (arrowR * np.sin(a), [(arrowR+0.5) * np.sin(a[-20])]) )
		arrowY = np.concatenate( (arrowR * np.cos(a), [(arrowR+0.5) * np.cos(a[-20])]) )
		arrowZ = np.zeros(arrowX.shape)
		arrow = go.Scatter3d(x=arrowX, y=arrowY, z=arrowZ, mode='lines', line=dict(color='black', width=4))

		eyeD = 2.0
		eyeAz = -67
		eyeEl = 42
		eyePos = dict(x = eyeD*np.cos(eyeEl*np.pi/180)*np.sin(eyeAz*np.pi/180),
			y = eyeD*np.cos(eyeEl*np.pi/180)*np.cos(eyeAz*np.pi/180),
			z = eyeD*np.sin(eyeEl*np.pi/180))

		### layout format
		title_font = dict(size=int(35*fontscale), weight='bold')
		tick_kws = dict(showticklabels=True, \
			ticks='outside', tickcolor='#ffffff', tickwidth=8, ticklen=30, tickfont=dict(size=int(18*fontscale)), \
			gridcolor='#dddddd', gridwidth=5,
			showspikes = False, showbackground=False, backgroundcolor='white')

		#title_font = dict(size=19, weight='bold')
		#tick_kws = dict(ticks='outside', tickwidth=5, gridcolor='#dddddd', #zerolinecolor='#000000',zerolinewidth=3,
		#			showspikes = False, showbackground=False, backgroundcolor='white', tickfont=dict(size=15))

		layout = go.Layout(
			scene=dict(
				xaxis=dict(title = dict(text='', font=title_font), #X (kpc)
					range=[-8.01, 20.01], tickvals=np.arange(-8,20.1,4), **tick_kws),
				yaxis=dict(title=dict(text='', font=title_font), #Y (kpc)
					range=[-16.01, 20.01], tickvals=np.arange(-16,20.1,4), **tick_kws),
				zaxis=dict(title=dict(text='', font=title_font), #ΔZ (kpc)
					range=[-1.01, 1.01], tickvals=np.arange(-0.5, 2, 0.5), zerolinecolor='#000000',zerolinewidth=5, **tick_kws),
				annotations=[
					sunText,
					gcText,
					barText,
					*sectorText
					],
				camera = dict(
					up=dict(x=0, y=0, z=1),
					center=dict(x=0.10, y=-0.06, z=-0.32),
					eye=eyePos
					#eye=dict(x=-0.4, y=-1.6, z=1.2)
					),
				aspectmode='manual',
				aspectratio=dict(x=1, y=45/35., z=.4),
			),
			showlegend=False,
			width=int(18*400*figscale/2),	#textwidth * dpi * texfigurescale
			height=int(18*400*figscale*0.9/2),
			margin=dict(r=0, l=0, b=0, t=0),
		)

		# Create a Figure object and add the 3d objects to it
		fig = go.Figure(data=[residual, *wires, bar, sun, gc, arrow], layout=layout)#

		### panel ID
		#if component == 1:
		#	fig.add_annotation(x=0.05, y=0.88, xref='paper', yref='paper', text='a', \
		#		showarrow=False, font=dict(color='black', size=int(70), family="Arial Black"))

		### axis title
		fig.add_annotation(x=0.08, y=0.35, textangle=73, text='X (kpc)', xref='paper', yref='paper', \
			showarrow=False, font=dict(color='rgb(37,63,98)', size=int(40*fontscale), weight='bold'))
		fig.add_annotation(x=0.69, y=0.19, textangle=-24, text='Y (kpc)', xref='paper', yref='paper', \
			showarrow=False, font=dict(color='rgb(37,63,98)', size=int(40*fontscale), weight='bold'))
		fig.add_annotation(x=0.003, y=0.70, textangle=70, text='ΔZ (kpc)', xref='paper', yref='paper', \
			showarrow=False, font=dict(color='rgb(37,63,98)', size=int(40*fontscale), weight='bold'))

		### Show the plot
		if 0: fig.show()
		else: fig.write_image("fig/residual_%icomp.pdf" % (component))#, scale=3)


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



	### 3d visualization with pyvista
	if 0:
		import pyvista as pv

		### plot axes
		p = pv.Plotter(window_size=(1600,1600))#, lighting=None)

		if 1:
			### data + model

			### DATA points
			#p.add_points(pv.PointSet(np.column_stack(X, Y, Z)), style='points_gaussian', scalars=np.log10(mass), cmap='gray', clim=[0,5], opacity=0.2, point_size=5, show_scalar_bar=False)
			#p.add_points(pv.PointSet(np.column_stack(X, Y, Z)), scalars=Z, style='points_gaussian', render_points_as_spheres=False, emissive=False, cmap='RdYlBu_r', clim=[-1,1], opacity=0.3, point_size=5, show_scalar_bar=False)

			### DATA regrid surface
			### use a small step for surface to make it smooth, eg. (0.25, 2)
			### use a larger step for wireframe to make it clear to see, eg. (0.5, 4)
			#gridX, gridY, gridZ = convolveZonGrid(X, Y, Z, mass, step=(0.75, 5), polar=True, warp=False, sin=False)
			#gridX, gridY, gridZ = _dataConvolve(step=(0.5, 4), polar=True)
			#p.add_mesh(pv.StructuredGrid(gridX, gridY, gridZ), style='points_gaussian', color='k', clim=[-1,1], opacity=1.0, point_size=14, render_points_as_spheres=False)
			
			#p.add_points(pv.PointSet(np.column_stack(gridX.ravel(), gridY.ravel(), gridZ.ravel())), style='points_gaussian', scalars=gridZ.ravel(), cmap='RdYlBu_r', clim=[-1,1], opacity=1.0, point_size=3, render_points_as_spheres=False)
			#p.add_mesh(pv.StructuredGrid(gridX, gridY, gridZ), style='surface', scalars=gridZ.T, cmap='RdYlBu_r', clim=[-1,1], opacity=0.9, point_size=5, show_scalar_bar=False)

			### MODEL wireframe
			#gridX, gridY, gridZ = _modelGrid(sin=False)
			grid = pv.StructuredGrid(gridX, gridY, gridZ)
			grid.texture_map_to_plane(inplace=True)

			p.add_mesh(grid, style='surface', opacity=1.0, line_width=1, texture = pv.read_texture('ssc2008-10a_ext.jpg'))
			#p.add_mesh(pv.StructuredGrid(gridX, gridY, gridZ), style='surface', scalars=gridZ.T, cmap='RdYlBu_r', clim=[-1,1], point_size=10)

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
				lineResidual = convolveZ(X, Y, Z, mass, lineX, lineY, kernel=0.5) - function((lineR, linePHI), best, sin=False)
				xyr = np.vstack((lineX, lineY, lineResidual)).T
				p.add_lines(xyr, color='k', width=2)
				#p.add_mesh(pv.StructuredGrid(lineX, lineY, lineResidual), style='wireframe', scalars=lineResidual, color='RdYlBu_r', clim=[-.05,.05], opacity=1.0, line_width=5 if i%2==1 else 0)

			for r in np.arange(8, 22.1, 1):
				linePHI = np.arange(-30, 180, 1)
				lineR = np.repeat(r, linePHI.size)
				lineX, lineY = RPHI2XY(lineR, linePHI)
				lineResidual = convolveZ(X, Y, Z, mass, lineX, lineY, kernel=0.5) - function((lineR, linePHI), best, sin=False)
				xyr = np.vstack((lineX, lineY, lineResidual)).T
				p.add_lines(xyr, color='k', width=2)
				#p.add_mesh(pv.StructuredGrid(lineX, lineY, lineResidual), style='wireframe', scalars=lineResidual, color='RdYlBu_r', clim=[-.05,.05], opacity=1.0, line_width=5 if r%1==0 else 0)
			
			gridX, gridY, gridZ = _residualGrid(step=(0.2, 1), polar=True)
			#p.add_mesh(pv.StructuredGrid(gridX, gridY, gridZ), style='wireframe', scalars=gridZ.T, cmap='RdYlBu_r', clim=[-.1,.1], opacity=1.0, line_width=3)
			p.add_mesh(pv.StructuredGrid(gridX, gridY, gridZ), style='surface', scalars=gridZ.T, cmap='RdYlBu_r', clim=[-.1,.1], opacity=0.7)
		
			#gridX, gridY, gridZ = _modelGrid(warp=False)
			#p.add_mesh(pv.StructuredGrid(gridX, gridY, gridZ), style='wireframe', color='k', opacity=1.0, line_width=2)
			#p.add_mesh(pv.StructuredGrid(gridX, gridY, gridZ), style='surface', scalars=gridZ.T, cmap='RdYlBu_r', clim=[-1,1], point_size=10)

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



	### 3d visualization with mayavi
	if 0:
		from mayavi import mlab

		zscale = 5
		gridX, gridY, gridZ = convolveZonGrid(X, Y, Z, mass, step=(0.75, 5), polar=True)


		# Start a new figure
		mlab.figure(bgcolor=(1, 1, 1), size=(800, 600))

		# Plot scatter points with color
		pts = mlab.points3d(gridX, gridY, gridZ*zscale, gridZ, mode='sphere', scale_mode='none', scale_factor=0.7, \
			colormap='jet', vmin=-0.5, vmax=0.5, opacity=1)

		cmap = colormaps['RdYlBu_r']
		clist = (cmap(np.arange(256))[:,:]*255).astype(int)
		pts.module_manager.scalar_lut_manager.lut.table=clist
		'''
		scene = mlab.gcf()
		scene.scene.camera.view_up = (0, 0, 1)  # optional, set view up
		scene.scene.render()

		# Access the root transform (important!)
		transform = pts.scene.children[0].children[0].children[0].children[0].actor.actor.scale
		transform[:] = [1, 1, 2]  # x, y, z scale
		'''


		### R grid
		from shared import rad_sep
		for phi in list(range(323, 155, -12))+rad_sep:
			#for phi in np.arange(0, 360, 15):
			lineR = np.arange(8.5, 20.6, 0.1)
			linePHI = np.repeat(phi, lineR.size)
			lineX, lineY = RPHI2XY(lineR, linePHI)
			lineZ = function((lineR, linePHI), bestmed, warp=True)
			mlab.plot3d(lineX, lineY, lineZ*zscale, color=(0,0,0), tube_radius=None, line_width=2)

		### phi grid
		for r in [8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 17.5, 20.5]:
			linePHI = np.linspace(0, 360, int(360*r)//40)
			lineR = np.repeat(r, linePHI.size)
			lineX, lineY = RPHI2XY(lineR, linePHI)
			lineZ = function((lineR, linePHI), bestmed, warp=True)
			mlab.plot3d(lineX, lineY, lineZ*zscale, color=(0,0,0), tube_radius=None, line_width=2)

		
		### axis
		x0,x1,dx = -24,+24,8
		y0,y1,dy = -24,+24,8
		z0,z1,dz = -2,2,1
		for x in np.arange(x0, x1+0.1, dx):
			mlab.plot3d([x, x], [y0, y1], [z0*zscale, z0*zscale], color=(0.7,0.7,0.7), line_width=1, tube_radius=None)
			mlab.plot3d([x, x], [y0, y0], [z0*zscale, z1*zscale], color=(0.7,0.7,0.7), line_width=1, tube_radius=None)
			mlab.text3d(x, y1+4, z0*zscale, '%i' % x, color=(0,0,0), line_width=2)
		for y in np.arange(y0, y1+0.1, dy):
			mlab.plot3d([x0, x1], [y, y], [z0*zscale, z0*zscale], color=(0.7,0.7,0.7), line_width=1, tube_radius=None)
			mlab.plot3d([x1, x1], [y, y], [z0*zscale, z1*zscale], color=(0.7,0.7,0.7), line_width=1, tube_radius=None)
			mlab.text3d(x0-2, y, z0*zscale, '%i' % y, color=(0,0,0), line_width=2)
		for z in np.arange(z0, z1+0.1, dz):
			mlab.plot3d([x0, x1], [y0, y0], [z*zscale, z*zscale], color=(0.7,0.7,0.7) if z!=0 else (0,0,0), line_width=1, tube_radius=None)
			mlab.plot3d([x1, x1], [y0, y1], [z*zscale, z*zscale], color=(0.7,0.7,0.7) if z!=0 else (0,0,0), line_width=1, tube_radius=None)
			mlab.text3d(x1, y1+2, z*zscale, '%i' % (z), color=(0,0,0), line_width=2)
		

		### r phi relation
		if component==1:
			_r, _phi, _z = [], [], []
			for i in range(9,18):
				suffix = '_%i_%i' % (i, i+1)
				steps = np.load(os.path.join(path, 'steps%s.npy' % suffix))
				probs = np.load(os.path.join(path, 'probs%s.npy' % suffix))
				bestmedi = np.mean(steps[probs > np.percentile(probs, 95)], axis=0)
				_r.append(i+0.5)
				_phi.append(bestmedi[-1])
				_z.append(0)
			_x, _y = RPHI2XY(np.array(_r), np.array(_phi))
			mlab.plot3d(_x, _y, _z, color=(0.3,0.6,0.3), line_width=2, tube_radius=None)

			### plot zero line
			_r = [8.5, 20.5]
			_phi = [bestmed[-1]]*2
			_x, _y = RPHI2XY(np.array(_r), np.array(_phi))
			mlab.plot3d(_x, _y, [0,0], color=(0.3,0.6,0.3), line_width=2, tube_radius=None)

		
		### phi sector text
		from shared import rad_sep
		phi_sep = rad_sep[:-3]+[None]+rad_sep[-3:]
		label = 0
		for i in range(len(phi_sep)-1):
			if phi_sep[i] is None or phi_sep[i+1] is None: continue
			label += 1
			_r = 22.5
			_phi = (phi_sep[i] + phi_sep[i+1])/2
			_x, _y = RPHI2XY(_r, _phi)
			_z = function((np.array([_r]), np.array([_phi])), bestmed, warp=True)[0]
			if _z > 2.21: _z = 2.21
			mlab.text3d(_x, _y, _z*zscale, '₁', color=(0,0,0), line_width=2)


		### radius text
		sectorR = np.array([(p1+p2)/2 for p1,p2 in zip([8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 17.5], [9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 17.5, 20.5])])
		#sectorR = np.array([(p1+p2)/2 for p1,p2 in zip([8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 18.5], [9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 18.5, 20.5])])
		sectorPhi = np.repeat(-40, len(sectorR))
		sectorZ = function((sectorR, sectorPhi), bestmed)
		sectorX, sectorY = RPHI2XY(sectorR, sectorPhi)
		for i in [0, 8]:#range(len(sectorR)):
			mlab.text3d(sectorX[i], sectorY[i], sectorZ[i], text='$\huge R_{%i}$' % (i+1))

		#draw_grid()
		mlab.view(azimuth=125, elevation=90-26)

		# Show the plot
		mlab.show()



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
	


	### interactive online figure
	if 0 and not sin:
		from shared import *
		import plotly.graph_objs as go
		figscale = 1
		fontscale = 1.3

		### Data
		#points = go.Scatter3d(x=X, y=Y, z=Z - function((R,PHI), bestmed, sin=False), mode='markers', \
		#	marker=dict(size=3, color=Z, colorscale='RdYlBu_r', cmin=-1, cmax=1, opacity=0.2, colorbar=dict(thickness=20)))
		markersize = np.log10(mass)*3
		markersize -= np.nanmin(markersize)-1
		clouds = go.Scatter3d(x=X, y=Y, z=Z, mode='markers', visible=False, \
			marker=dict(size=markersize, color=Z, colorscale='RdYlBu_r', cmin=-0.5, cmax=0.5, opacity=0.5, line=dict(width=0), \
			colorbar=dict(thickness=25, title=dict(text='Z (kpc)', font=dict(size=int(20*fontscale), weight='bold')), \
				x=0.24, y=0.12, len=0.3, orientation='h', tickvals=np.arange(-0.5,1,0.25), tickfont=dict(size=int(18*fontscale)))))

		### median
		gridX, gridY, gridZ = convolveZonGrid(X, Y, Z, mass, step=(0.75, 5), polar=True, useMask=False)
		gridR, gridPHI = XY2RPHI(gridX, gridY)

		#gridZ = gridZ - function((gridR, gridPHI), bestmed, warp1=True, sin=False)
		#gridR, gridPHI, gridZ = convolveZonGrid(X, Y, Z, mass, step=(0.5, 0.5), polar=False)
		meds = go.Scatter3d(x=gridX.ravel(), y=gridY.ravel(), z=gridZ.ravel(), mode='markers', \
			marker=dict(size=6, color=gridZ.ravel(), colorscale='RdYlBu_r', cmin=-0.5, cmax=0.5, opacity=1.0, \
			colorbar=dict(thickness=25, title=dict(text='Z (kpc)', font=dict(size=int(20*fontscale), weight='bold')), \
				x=0.24, y=0.12, len=0.3, orientation='h', tickvals=np.arange(-0.5,1,0.25), tickfont=dict(size=int(18*fontscale)))))


		### wireframe
		wires = []
		line_marker = dict(color='#000000', width=5)#, colorscale='RdYlBu_r', cmin=-0.5, cmax=0.5)
		### R grid
		from shared import rad_sep
		for phi in list(range(323, 155, -12))+rad_sep:
			#for phi in np.arange(0, 360, 15):
			lineR = np.arange(8.5, 20.6, 0.5)
			linePHI = np.repeat(phi, lineR.size)
			lineX, lineY = RPHI2XY(lineR, linePHI)
			lineZ = function((lineR, linePHI), bestmed, warp=True, sin=False)
			#line_marker['color'] = lineZ
			wires.append(go.Scatter3d(x=lineX, y=lineY, z=lineZ, mode='lines', line=line_marker))

		### phi grid
		for r in [8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 17.5, 20.5]:
			linePHI = np.linspace(0, 360, int(360*r)//40)
			lineR = np.repeat(r, linePHI.size)
			lineX, lineY = RPHI2XY(lineR, linePHI)
			lineZ = function((lineR, linePHI), bestmed, warp=True, sin=False)
			#line_marker['color'] = lineZ
			wires.append(go.Scatter3d(x=lineX, y=lineY, z=lineZ, mode='lines', line=line_marker))


		### phi sector text
		sectorPhi = np.array([(p1+p2)/2 for p1,p2 in zip(rad_sep[:-1], rad_sep[1:]) if p1*p2>0])
		sectorR = np.repeat(22.5, len(sectorPhi))
		sectorZ = function((sectorR, sectorPhi), bestmed)
		sectorX, sectorY = RPHI2XY(sectorR, sectorPhi)
		sectorZ[sectorZ>2.21] = 2.21
		sectorT = ['ϕ'+to_subscript(i+1) for i in range(len(sectorPhi))]
		sectorText = go.Scatter3d(x=sectorX, y=sectorY, z=sectorZ, mode='text', \
			text=sectorT, textposition='middle center', textfont=dict(color='#444444', size=int(18*fontscale)))

		### radius text
		annulusR = np.array([(p1+p2)/2 for p1,p2 in zip([8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 17.5], [9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 17.5, 20.5])])
		annulusPhi = np.repeat(-40, len(annulusR))
		annulusZ = function((annulusR, annulusPhi), bestmed)
		annulusX, annulusY = RPHI2XY(annulusR, annulusPhi)
		annulusT = np.array(['R'+to_subscript(i+1) for i in range(len(annulusR))])
		useIdx = [0,-1]
		annulusText = go.Scatter3d(x=annulusX[useIdx], y=annulusY[useIdx], z=annulusZ[useIdx], mode='text', \
			text=annulusT[useIdx], textposition='middle center', textfont=dict(color='#444444', size=int(15*fontscale)))

		
		sun = go.Scatter3d(x=[0], y=[8.15], z=[0], mode='markers', marker=dict(color='#ff6600', size=6))
		sunText = dict(x=0, y=8.15, z=0, text='Sun', showarrow=False,
			xanchor='left', xshift=5, yanchor='bottom', font=dict(color='#ff6600', size=int(18*fontscale)))

		gc = go.Scatter3d(x=[0], y=[0], z=[0], mode='markers', marker=dict(color='#555555', size=6))
		gcText = dict(x=0, y=0, z=0, text='GC', showarrow=False,
			xanchor='left', xshift=5, yanchor='bottom', font=dict(color='#666666', size=int(18*fontscale)))

		
		### R vs phi_w
		if component==1:
			_r, _phi, _z = [], [], []
			for i in range(9,18):
				suffix = '_%i_%i' % (i, i+1)
				steps = np.load(os.path.join(path, 'steps%s.npy' % suffix))
				probs = np.load(os.path.join(path, 'probs%s.npy' % suffix))
				bestmedi = np.mean(steps[probs > np.percentile(probs, 95)], axis=0)
				_r.append(i+0.5)
				_phi.append(bestmedi[-1])
				_z.append(0)
			_x, _y = RPHI2XY(np.array(_r), np.array(_phi))
			_rphiw = go.Scatter3d(x=_x, y=_y, z=_z, mode='lines+markers',\
				marker=dict(symbol='square', size=6, color='darkgreen', opacity=1), line=dict(color='darkgreen', width=12))
			_r = [8.5, 20.5]
			_phi = [bestmed[-1]]*2
			_x, _y = RPHI2XY(np.array(_r), np.array(_phi))
			_phiw_1comp = go.Scatter3d(x=_x, y=_y, z=[0, 0], mode='lines',\
				marker=dict(symbol='square', size=0.1, color='darkgreen', opacity=1), line=dict(dash='longdash', color='darkgreen', width=10))


		### eye pos
		eyeD = 2
		eyeAz = -35
		eyeEl = 26
		eyePos = dict(x = eyeD*np.cos(eyeEl*np.pi/180)*np.sin(eyeAz*np.pi/180),
			y = eyeD*np.cos(eyeEl*np.pi/180)*np.cos(eyeAz*np.pi/180),
			z = eyeD*np.sin(eyeEl*np.pi/180))

		### layout format
		title_font = dict(size=int(20*fontscale), weight='bold')
		tick_kws = dict(showticklabels=True, \
			ticks='outside', tickcolor='#ffffff', tickwidth=8, ticklen=15, tickfont=dict(size=int(12*fontscale)), \
			gridcolor='#dddddd', gridwidth=5,
			showspikes = False, showbackground=False, backgroundcolor='white')

		layout = go.Layout(
			scene=dict(
				xaxis=dict(title = dict(text='X (kpc)', font=title_font),
					range = [-24.01, 24.01],
					tickvals = np.arange(-24, 25, 8),
					#ticktext = ['%+i' % x for x in range(-24, 25, 8)],
					**tick_kws),
				yaxis=dict(title = dict(text='Y (kpc)', font=title_font),
					range = [-24.01, 24.01],
					tickvals = np.arange(-24, 25, 8), 
					#ticktext = ["%+i" % y for y in range(-24, 25, 8)],
					**tick_kws),
				zaxis=dict(title = dict(text='Z (kpc)', font=title_font),
					range = [-2.1, 2 if component==1 else 2.21],
					tickvals = np.arange(-1, 3, 1),
					#ticktext = ["%+i" % z for z in range(-1, 3, 1)],
					zerolinecolor='#000000', zerolinewidth=5, **tick_kws),
				annotations=[sunText, gcText],
				aspectmode='manual',
				aspectratio=dict(x=1, y=1, z=.45),
				camera = dict(
					up=dict(x=0, y=0, z=1),
					center=dict(x=0.11, y=0, z=-0.3),
					eye=eyePos
					),
				),
			showlegend=False,
			width=int(6*400*figscale/2),	#textwidth * dpi * texfigurescale
			height=int(6*400*figscale*0.8/2),
			margin=dict(r=0, l=0, b=0, t=0, pad=20),
			#scene_aspectmode='manual',
			#scene_aspectratio=dict(x=1, y=1, z=.45),
		)

		# Create a Figure object and add the 3d objects to it
		if component==1:
			trace = [clouds, meds, *wires, sectorText, annulusText, sun, gc, _phiw_1comp]#_rphiw
		else:
			trace = [clouds, meds, *wires, sectorText, annulusText, sun, gc]
		fig = go.Figure(data=trace, layout=layout)#


		fig.update_layout(
			updatemenus=[dict(
				type="buttons",
				direction="right",
				buttons=[
					dict(label='Mean', method='restyle', args=[dict(visible=[False, True] + [True]*(len(trace)-2))]),
					dict(label='Data', method='restyle', args=[dict(visible=[True, False] + [True]*(len(trace)-2))]),
					],
				showactive=True,
				x=0.0,
				y=1.15,
				xanchor="right",
				yanchor="top"
			)],
			margin=dict(l=0, r=0, t=50, b=0),
		)


		### Show the plot
		if 0: fig.show()
		else:
			fig.write_image("fig/median_model_%icomp.png" % component)#, scale=3)
			fig.write_html("fig/median_model_%icomp.html" % component)#, scale=3)


	if 0 and sin:
		figscale = 1
		fontscale = 1.3

		from shared import *
		import plotly.graph_objs as go
		
		
		### residual
		resX, resY, resZ = convolveZonGrid(X, Y, Z, mass, step=(0.25, 1), polar=True)
		#resX, resY, resZ = _excludedXYZ(resX, resY, resZ)
		residual = go.Surface(x=resX, y=resY, z=resZ, surfacecolor=resZ, colorscale='RdYlBu_r', cmin=-0.2, cmax=0.2, opacity=1,
			colorbar=dict(thickness=25, title=dict(text='Z (kpc)', font=dict(size=int(20*fontscale), weight='bold')), \
				x=0.72, y=0.12, len=0.3, orientation='h', tickvals=np.arange(-0.2,0.5,0.1), tickfont=dict(size=int(18*fontscale))))
		#residual = go.Scatter3d(x=resX.ravel(), y=resY.ravel(), z=resZ.ravel(), mode='markers', \
		#	marker=dict(size=5, color=resZ.ravel(), colorscale='RdYlBu_r', cmin=-0.2, cmax=0.2, opacity=1.0, \
		#	colorbar=dict(thickness=20, title='Z (kpc)', x=0.24, y=0.7, len=0.05, len=0.3, orientation='h')))
		
		
		### residual wireframe
		wires = []
		line_marker = dict(color='#223322', width=5)
		#for phi in np.arange(0, 360, 15):
		from shared import rad_sep
		for phi in list(range(323, 155, -12))+rad_sep:
			lineR = np.arange(8, 22.1, 0.25)
			linePHI = np.repeat(phi, lineR.size)
			lineX, lineY = RPHI2XY(lineR, linePHI)
			lineZ = convolveZ(X, Y, Z, mass, lineX, lineY)
			if phi in [29,38]: line_marker['width']=15
			else: line_marker['width']=5
			wires.append(go.Scatter3d(x=lineX, y=lineY, z=lineZ+0.01, mode='lines', line=line_marker))

		for r in np.arange(8, 22.1, 1):
			#linePHI = np.linspace(0, 360, int(360*r)//40)
			linePHI = np.arange(0, 360.1, 1)
			lineR = np.repeat(r, linePHI.size)
			lineX, lineY = RPHI2XY(lineR, linePHI)
			lineZ = convolveZ(X, Y, Z, mass, lineX, lineY)
			if r in [12,]: line_marker['width']=15
			else: line_marker['width']=5
			wires.append(go.Scatter3d(x=lineX, y=lineY, z=lineZ+0.01, mode='lines', line=line_marker))
		
		
		#sectorPhi = np.arange(165, -35, -15)-7.5
		sectorPhi = np.array([(p1+p2)/2 for p1,p2 in zip(rad_sep[:-1], rad_sep[1:]) if p1*p2>0])
		sectorR = np.linspace(17, 19, sectorPhi.size)
		sectorZ = np.zeros(sectorPhi.size)
		sectorX, sectorY = RPHI2XY(sectorR, sectorPhi)
		sectorT = ['ϕ'+to_subscript(i+1) for i in range(len(sectorPhi))]
		sectorText = go.Scatter3d(x=sectorX, y=sectorY, z=sectorZ, mode='text', \
			text=sectorT, textposition='middle center', textfont=dict(color='#444444', size=int(18*fontscale)))
		sectorText = []
		for i in range(len(sectorPhi)):
			sectorText.append(dict(x=sectorX[i], y=sectorY[i], z=sectorZ[i], text='ϕ'+to_subscript(i+1), showarrow=False, \
				xanchor = 'center', xshift=5, yanchor='middle', font=dict(color='#000000', size=int(18*fontscale))))
		

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
			xanchor='left', xshift=5, yanchor='bottom', font=dict(color='#666666', size=int(18*fontscale)))

		sun = go.Scatter3d(x=[0], y=[8.15], z=[0], mode='markers', \
				marker=dict(color='#ff6600', size=8))
		sunText = dict(x=0, y=8.15, z=0, text='Sun', showarrow=False,
			xanchor='left', xshift=5, yanchor='bottom', font=dict(color='#ff6600', size=int(18*fontscale)))

		gc = go.Scatter3d(x=[0], y=[0], z=[0], mode='markers', \
			marker=dict(color='#555555', size=8))
		gcText = dict(x=0, y=0, z=0, text='GC', showarrow=False,
			xanchor='left', xshift=5, yanchor='bottom', font=dict(color='#666666', size=int(18*fontscale)))

		### arrow
		arrowR = 21
		a = np.linspace(24.5, 53, 200) /180*np.pi
		arrowX = np.concatenate( (arrowR * np.sin(a), [(arrowR+0.5) * np.sin(a[-20])]) )
		arrowY = np.concatenate( (arrowR * np.cos(a), [(arrowR+0.5) * np.cos(a[-20])]) )
		arrowZ = np.zeros(arrowX.shape)
		arrow = go.Scatter3d(x=arrowX, y=arrowY, z=arrowZ, mode='lines', line=dict(color='black', width=4))

		eyeD = 2.2
		eyeAz = -67
		eyeEl = 42
		eyePos = dict(x = eyeD*np.cos(eyeEl*np.pi/180)*np.sin(eyeAz*np.pi/180),
			y = eyeD*np.cos(eyeEl*np.pi/180)*np.cos(eyeAz*np.pi/180),
			z = eyeD*np.sin(eyeEl*np.pi/180))

		### layout format
		title_font = dict(size=int(20*fontscale), weight='bold')
		tick_kws = dict(showticklabels=True, \
			ticks='outside', tickcolor='#ffffff', tickwidth=8, ticklen=20, tickfont=dict(size=int(12*fontscale)), \
			gridcolor='#dddddd', gridwidth=5,
			showspikes = False, showbackground=False, backgroundcolor='white')

		#title_font = dict(size=19, weight='bold')
		#tick_kws = dict(ticks='outside', tickwidth=5, gridcolor='#dddddd', #zerolinecolor='#000000',zerolinewidth=3,
		#			showspikes = False, showbackground=False, backgroundcolor='white', tickfont=dict(size=15))

		layout = go.Layout(
			scene=dict(
				xaxis=dict(title = dict(text='X (kpc)', font=title_font), #
					range=[-8.01, 20.01], tickvals=np.arange(-8,20.1,4), **tick_kws),
				yaxis=dict(title=dict(text='Y (kpc)', font=title_font), #
					range=[-16.01, 20.01], tickvals=np.arange(-16,20.1,4), **tick_kws),
				zaxis=dict(title=dict(text='ΔZ (kpc)', font=title_font), #
					range=[-1.01, 1.01], tickvals=np.arange(-0.5, 2, 0.5), zerolinecolor='#000000',zerolinewidth=5, **tick_kws),
				annotations=[
					sunText,
					gcText,
					barText,
					*sectorText
					],
				camera = dict(
					up=dict(x=0, y=0, z=1),
					center=dict(x=0.10, y=-0.06, z=-0.32),
					eye=eyePos
					#eye=dict(x=-0.4, y=-1.6, z=1.2)
					),
				aspectmode='manual',
				aspectratio=dict(x=1, y=45/35., z=.4),
			),
			showlegend=False,
			width=int(6*400*figscale/2),	#textwidth * dpi * texfigurescale
			height=int(6*400*figscale*0.9/2),
			margin=dict(r=0, l=0, b=0, t=0),
		)

		# Create a Figure object and add the 3d objects to it
		fig = go.Figure(data=[residual, *wires, bar, sun, gc, arrow], layout=layout)#


		### Show the plot
		if 0: fig.show()
		else:
			fig.write_image("fig/dZ_%icomp.png" % (component))#, scale=3)
			fig.write_html("fig/dZ_%icomp.html" % (component))#, scale=3)


	plt.show()

