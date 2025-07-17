## this is a test script to plot faceon Galactic plane with arms.
## the test results will apply to warp plot of each arm
import numpy as np
import pylab as pl
import scipy.interpolate as sp_interp
from scipy.stats import norm
import scipy.optimize as sp_opt
from  Distance import KinUnit, Distance
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from shared import *


def RPHI2XY(r, phi):
	phi = np.array(phi)
	x = r*np.sin(phi/180*np.pi)
	y = r*np.cos(phi/180*np.pi)
	return x, y

def plotArm(ax, phiRange, p, **kws):
	phi = np.linspace(*phiRange, 1400)
	r = function_arm(phi, p)
	print(r)
	x = r*np.sin(phi/180*np.pi)
	y = r*np.cos(phi/180*np.pi)
	ax.plot(x, y, **kws)

fig, ax = plt.subplots()

arm_kws = dict(lw=0.5)
plotArm(ax, (-170, 220), [0, 6.5, 12], color='k', **arm_kws) #sag
plotArm(ax, ( -50, 45), [0, 8.15, 12], color='pink', **arm_kws) #loc
plotArm(ax, ( -50, 270), [0, 10, 12], color='darkgreen', **arm_kws) #per
plotArm(ax, ( -48, 420), [0, 12, 12], color='darkred', **arm_kws) #out
plotArm(ax, ( -45, 450), [0, 18, 12], color='darkblue', **arm_kws) #osc

#plotArm(ax, (-22.7, 80.4), [90, 8.5, 12], lw=3, color='r')
#plotArm(ax, (-29.8, 153.4), [180, 8.5 ,12], lw=3, color='g')
#plotArm(ax, (-30.3, 158.8), [270, 8.5, 12], lw=3, color='b')

ax.plot(0, 8.15, 'ro')
ax.text(0, 8.15, 'Sun', ha='left', va='bottom', color='r')

ax.plot(0, 0, 'ko')
ax.text(0, 0, 'GC', ha='left', va='bottom')

### solar circle
ax.plot(*RPHI2XY(np.repeat(8.15, 200), np.linspace(0, 360, 200)), '--')


ax.set_aspect('equal')

plt.show()