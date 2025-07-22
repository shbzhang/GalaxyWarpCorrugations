import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from scipy.optimize import root_scalar
from scipy.interpolate import interp1d
import pyvista as pv
from scipy.spatial import cKDTree


def cmapStar():
	# a cmap of stars in HR diagram
	c = np.array([(176, 30, 27), (248, 50, 24), (251, 155, 39), (253, 206, 51), (255, 246, 136), (255, 249, 193), (209, 246, 254), (190, 255, 255)], dtype=float)
	c /= 255
	return mcolors.LinearSegmentedColormap.from_list('Star', c, N=256)


def cmapCloud():
	# a cmap of clouds
	c = np.array([(218, 115, 99), (143, 121, 117), (175, 208, 230), (105, 120, 159), (0,0,0)], dtype=float)
	c /= 255
	return mcolors.LinearSegmentedColormap.from_list('Cloud', c, N=256)


import numpy as np

def pdf(x, n):
	return ((n + 1) / n) * (1 - x**n)

def cdf(x, n):
	return ((n + 1)*x - x**(n + 1)) / n


def generateStars(n, k=2, length=100, width=0.1, height=0.1, ageK=3):
	# Distribution along arm
	# f(x) ~ 1-x^k
	# pdf(x) = (k+1)/k * (1-x^k)
	# cdf(x) = ( (k+1) * x - x^(k+1) ) / k
	_x = np.linspace(0,1,100)
	_y = ( (k+1) * _x - _x**(k+1) ) / k
	inverseCDF = interp1d(_y, _x, kind='linear', bounds_error=False, fill_value=(0, 1))

	u = np.random.uniform(0, 1, n)
	normX = u#inverseCDF(u)

	X = normX * length
	Y = np.random.normal(0, width, n)# * (normX*2+1) #Y get wider
	Z = np.random.normal(0, height, n)# * (normX*2+1) #Z get wider
	stars = np.vstack((X, Y, Z)).T
	ages = np.random.uniform(0, 1, n)**ageK #ages=0 is young

	return stars, ages


def generateRegions(nStars=2000, nRegions=7, **kws):
	# move stars toward SF regions
	# return moved stars
	stars, ages = generateStars(nStars, **kws)
	regions, _ = generateStars(nRegions, **kws)

	#find nearest neighbour
	tree = cKDTree(regions)
	distances, indices = tree.query(stars)
	
	#vector to nearest neighbour
	distancesVector = regions[indices] - stars
	clusteredStars = stars# + distancesVector * (1-ages[:, np.newaxis])
	return clusteredStars, ages


def generateArm1(R0, phi0, pitch, **kws):
	# transform stars in arm coordiante to galactocentric coordinate system
	# R = R0 * exp(-(phi-phi0) * tan(pitch))
	#   = exp(phi0*tan(pitch)+lnR0) * exp(-phi*tan(pitch))

	#phi = phi0 + x/R0
	#R = 
	stars, ages = generateRegions(**kws)

	dPhi = (phi - phi0)/180*np.pi
	lnRtoR0 = -dPhi * np.tan(pitch/180*np.pi)

	lnR = lnRtoR0 + np.log(R0)
	R = np.exp(lnR)

	X = R * np.sin(phi)
	Y = R * np.cos(phi)
	Z = stars[:,2]
	armStars = np.vstack((X,Y,Z))

	return armStars, ages


def generateArm(R0, phi0, pitch_deg, **kws):
	#
	stars, ages = generateRegions(**kws)
	pitch = np.deg2rad(pitch_deg)
	# Position along the spiral
	phi = phi0 + stars[:,0] * np.tan(pitch) / R0
	R = R0 * np.exp(stars[:,0] * np.tan(pitch) / R0)
	
	# Add perpendicular offset y
	X = R * np.cos(phi) - stars[:,1] * np.sin(phi)
	Y = R * np.sin(phi) + stars[:,1] * np.cos(phi)
	armStars = np.vstack((X,Y,stars[:,2])).T

	return armStars, ages

def spiral_to_galactic(x, y, z, R_start, phi_start, pitch_deg):
	pitch = np.deg2rad(pitch_deg)
	# Position along the spiral
	phi = phi_start + x * np.tan(pitch) / R_start
	R = R_start * np.exp(x * np.tan(pitch) / R_start)
	
	# Add perpendicular offset y
	X = R * np.cos(phi) - y * np.sin(phi)
	Y = R * np.sin(phi) + y * np.cos(phi)
	
	return X, Y, z



from scipy.integrate import quad

def spiral_length(phi_start, phi_end, R0, pitch_deg):
	pitch = np.deg2rad(pitch_deg)
	tan_psi = np.tan(pitch)
	cos_psi = np.cos(pitch)
	
	L, err = quad(integrand, phi_start, phi_end)
	return L
def arm(R0, phi0, pitch):
	d2r = np.pi/180
	def arm(phi):
		return R0 * np.exp((phi - phi0) * np.tan(pitch*d2r))
	def armLength(phi):
		return R0 * np.exp((phi - phi0) * np.tan(pitch*d2r)) / np.cos(pitch*d2r)
	phiAxis = np.linspace(phi0, phi0+720, 180)
	lenAxis = [quad(armLength, phi0, v)[0] for v in phiAxis]
	len2phi = interp1d(lenAxis, phiAxis, kind='linear', bounds_error=False, fill_value=(0, np.inf))

	stars, ages = generateStars(1000)
	phiArm = len2phi(stars[:,0])
	RArm = arm(phiArm)

	R = RArm + stars[:,1] * np.cos(pitch*d2r)
	phi = phiArm + 1

	x = R * np.cos(phi)

	print(stars[:,0][np.isnan(phi)])


	plt.plot(phiAxis, lenAxis, '.-')
	plt.show()

arm(5,0,5)

p = pv.Plotter(window_size=(1600,1600))#
p.set_background('black')

stars, ages = generateRegions()
X,Y,Z = spiral_to_galactic(stars[:,0], stars[:,1], stars[:,2], 5, 0, 8)
stars = np.vstack((X,Y,Z)).T
#stars, ages = generateArm(0, 6.5, 12)
print(stars)
plt.scatter(stars[:,0], stars[:,1], c=ages, s=2, cmap=cmapStar())
plt.show()
#polyArm = pv.PolyData(stars)
polyArm = pv.PointSet(stars)
polyArm['color'] = ages
#polyArm['']
print(ages)

#p.add_mesh(polyArm, scalars='color', style='points', opacity=1, cmap=cmapStar())
p.add_points(polyArm, scalars='color', style='points', opacity=1, cmap=cmapStar())
p.show()

