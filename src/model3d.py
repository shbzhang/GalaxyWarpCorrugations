import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from scipy.optimize import root_scalar
from scipy.interpolate import interp1d
import pyvista as pv
from scipy.spatial import cKDTree
from scipy.integrate import quad, cumtrapz


d2r = np.pi/180
np.random.seed(45)

def cmapStar():
	# a cmap of stars in HR diagram
	c = np.array([(176, 30, 27), (248, 50, 24), (251, 155, 39), (253, 206, 51), (255, 246, 136), (255, 249, 193), (209, 246, 254), (190, 255, 255)], dtype=float)
	c /= 255
	return mcolors.LinearSegmentedColormap.from_list('Star', c, N=256)


def cmapCloud():
	# a cmap of clouds
	c = np.array([(218, 115, 99), (143, 121, 117), (175, 208, 230), (105, 120, 159), (105, 120, 159), (0,0,0)], dtype=float)
	c /= 255
	return mcolors.LinearSegmentedColormap.from_list('Cloud', c, N=256)



def pdf(x, n):
	return ((n + 1) / n) * (1 - x**n)

def cdf(x, n):
	return ((n + 1)*x - x**(n + 1)) / n


def generateStars(n=1000, k=1000, length=100, width=lambda x: 1, height=lambda x: 1, ageK=1):
	'''
	generate stars along a spiral arm
	input:
		n: number of stars
		k: index controls number density along arm
		   k>0 gives a fade tail, more uniform for k>>0
		   generate based on:
		      f(x) ~ 1-x^k
		      pdf(x) = (k+1)/k * (1-x^k)
		      cdf(x) = ( (k+1) * x - x^(k+1) ) / k
		length: length of arm
		width, height: a number or a function of width/height profile
		ageK: index controls age
			k>1 get more young stars
	output:
		coordinate and age of stars
		x length to the start point of arm
		y offset in plane
		z vertical offset
	'''
	_x = np.linspace(0,1,100)
	_y = ( (k+1) * _x - _x**(k+1) ) / k
	inverseCDF = interp1d(_y, _x, kind='linear', bounds_error=False, fill_value=(0, 1))

	u = np.random.uniform(0, 1, n)
	normX = inverseCDF(u) #range 0~1

	X = normX * length
	#Y = np.random.normal(0, 1, n)
	Y = np.random.uniform(0, 1, n)
	Y[Y<0.5] = np.log(2*Y[Y<0.5])
	Y[Y>=0.5] = -np.log(2*(1-Y[Y>=0.5]))
	if callable(width): Y *= width(normX) #width function
	else: Y *= width

	### density = Gaussian(z)
	#Z = np.random.normal(0, 1, n)
	### density = exp(-abs(x)), unit thick
	Z = np.random.uniform(0, 1, n)
	Z[Z<0.5] = np.log(2*Z[Z<0.5])
	Z[Z>=0.5] = -np.log(2*(1-Z[Z>=0.5]))
	#height is scaled when arm is generated, so don't set height for arms
	if callable(height): Z *= height(normX) #height function
	else: Z *= height
	stars = np.vstack((X, Y, Z)).T
	ages = np.random.uniform(0, 1, n)**ageK #ages=0 is young

	return stars, ages


def generateRegions(nRegions=10, moveK=0.2, **kws):
	'''
	generate SF regions in a spiral arm, move young star towards those regions
	input:
		nRegions: number of SF regions
		moveK:
			index controls distance to move
			higher k, old stars move more
	output:
		coordinate and age of moved stars
	'''
	stars, ages = generateStars(**kws)

	### SF regions
	if nRegions<=0: return stars, ages
	kws.pop('n')
	# shrink width and height
	width = kws.get('width')
	if width is None: kws['width'] = 0.6
	elif callable(width): kws['width'] = lambda x: width(x) * 0.6
	else: kws['width'] = width * 0.6
	height = kws.get('height')
	if height is None: kws['height'] = 0.6
	elif callable(height): kws['height'] = lambda x: height(x) * 0.6
	else: kws['height'] = height * 0.6
	regions, _ = generateStars(n=nRegions, **kws)

	#find nearest neighbour
	tree = cKDTree(regions)
	distances, indices = tree.query(stars)
	
	#vector to nearest neighbour
	distancesVector = regions[indices] - stars
	clusteredStars = stars + distancesVector * (1-ages[:, np.newaxis])**moveK
	return clusteredStars, ages


def generateArm(R0, phi_deg, pitch_deg, **kws):
	'''
	bend star to generate a spiral arm, move young star towards those regions
	view the Milky way from above
	positive Y axis from GC to Sun
	positive X axis from GC to right
	phi is counted from X axis, counterclockwise is positive
	Az is counted from Y axis, clockwise is positive
	thus: phi = 90 - Az
	input:
		R0: radius to the start of arm
		phi_deg: a list of phi in deg, indicate where the arm pitch angle changes
		pitch_deg: a list of pitch in deg, indicate the corresponding pitch angles
	output:
		coordinate and age of stars in arm
	'''
	phi = np.array(phi_deg) * d2r
	pitch = np.array(pitch_deg) * d2r

	def armFunc(phi, R0, phi0, pitch):
		return R0 * np.exp((phi - phi0) * np.tan(pitch))

	### get numerical length, interpolate inversely
	phiAxis = []
	RAxis = []
	pitchAxis = []
	lenAxis = []
	for i in range(len(phi)):
		if i < len(phi)-1:
			phiAxisI = np.linspace(phi[i], phi[i+1], 200)
		else:
			phiAxisI = np.linspace(phi[i], phi[i]+np.pi*2, 200)
		phiAxis += list(phiAxisI)

		if i == 0:
			RAxisI = armFunc(phiAxisI, R0, phi[i], pitch[i])
		else:
			RAxisI = armFunc(phiAxisI, RAxisI[-1], phi[i], pitch[i])
		RAxis += list(RAxisI)
		pitchAxis += [pitch[i]]*200

		dArc = RAxisI * np.sqrt(1 + np.tan(pitch[i])**2)

		if i == 0:
			lenAxisI = cumtrapz(dArc, phiAxisI, initial=0)
		else:
			lenAxisI = cumtrapz(dArc, phiAxisI, initial=0)+lenAxis[-1]
		lenAxis += list(lenAxisI)
	#plt.plot(lenAxis, np.array(phiAxisI)/d2r)
	#plt.show()

	len2phi = interp1d(lenAxis, phiAxis, kind='linear', bounds_error=False, fill_value=(phiAxis[0], phiAxis[-1]))
	len2R = interp1d(lenAxis, RAxis, kind='linear', bounds_error=False, fill_value=(RAxis[0], RAxis[-1]))
	len2pitch = interp1d(lenAxis, pitchAxis, kind='linear', bounds_error=False, fill_value=(pitchAxis[0], pitchAxis[-1]))

	### bend along arm
	stars, ages = generateRegions(**kws)
	#stars, ages = generateRegions(**kws)
	X, Y, Z = stars.T
	phi = len2phi(X)
	R = len2R(X)

	Xc = R * np.cos(phi)
	Yc = R * np.sin(phi)

	theta_perp = phi - len2pitch(X)
	armX = Xc + Y * np.cos(theta_perp)
	armY = Yc + Y * np.sin(theta_perp)
	### flare disk
	def flare(r):
		h = np.empty(r.shape)
		h[r<=8.15] = 0.25
		h[r>8.15] = 0.25*np.exp((r[r>8.15]-8.15)/5)
		return h
	armZ = Z * flare(R)
	armStars = np.vstack((armX, armY, armZ)).T

	return armStars, ages


def generateGalaxy():
	'''
	generate stars in a galaxy
	'''
	galaxy = []
	#galaxyStars, galaxyAges = [], []
	
	### ARMS
	#sag
	galaxy.append(generateArm(3.6, [-140, -90, -45, 0, 70, 200], [13, 7, 4, 7, 16, 10], length=65, n=6500, nRegions=65//2, \
		width=lambda x: np.exp(x/1.2)-0.95, k=1, ageK=0.8, moveK=2.5))
	#per
	galaxy.append(generateArm(3.6, [-170, -90, -40, 0, 150], [14, 13, 10, 11, 15], length=54, n=5400, nRegions=54//2, \
		width=lambda x: np.exp(x/1.2)-0.95, k=1, ageK=0.8, moveK=2.5))

	#out
	galaxy.append(generateArm(3.4, [-320, -240, -170, -90, 80], [9, 6, 15, 10, 8], length=65, n=6500, nRegions=65//2, \
		width=lambda x: np.exp(x/1.2)-0.95, k=1, ageK=0.8, moveK=2.5))

	#osc
	galaxy.append(generateArm(3.6, [-350, -300, -180, -120, -80], [12, 14, 12, 11, 8], length=80, n=8000, nRegions=80//2, \
		width=lambda x: np.exp(x/1.2)-0.95, k=1, ageK=0.8, moveK=2.5))
	
	
	### SPURS
	spurKws = dict(width=0.25, k=5, ageK=1, moveK=3)
	# Local spur
	galaxy.append(generateArm(7.4, [30, 70], [3, 7], length=16, n=300, nRegions=4, **spurKws))
	# Scutum–Centaurus Spur + Sagittarius–Carina Spur
	galaxy.append(generateArm(3.7, [30, 90, 150], [28, 0, 34], length=23, n=800, nRegions=12, **spurKws))
	galaxy.append(generateArm(4.5, [-210,], [32,], length=5, n=100, nRegions=10, **spurKws))
	galaxy.append(generateArm(4.6, [-180,], [30,], length=7, n=100, nRegions=10, **spurKws))
	galaxy.append(generateArm(5.8, [-150,], [30,], length=8, n=100, nRegions=10, **spurKws))
	galaxy.append(generateArm(7.0, [-100,], [25,], length=9, n=100, nRegions=10, **spurKws))
	galaxy.append(generateArm(8.0, [-50,], [25,], length=10, n=100, nRegions=10, **spurKws))
	galaxy.append(generateArm(10.0, [10,], [23,], length=12, n=100, nRegions=10, **spurKws))
	galaxy.append(generateArm(6.2, [-50,], [23,], length=5, n=100, nRegions=10, **spurKws))
	galaxy.append(generateArm(7.5, [-170,], [35,], length=10, n=100, nRegions=10, **spurKws))
	
	### RING
	ringStars, ringAges = generateArm(3.1, [0], [0], n=3000, nRegions=3, length=19.5, width=0.4, height=2,  k=1000, ageK=0.1, moveK=3)
	galaxy.append((ringStars, ringAges*0.8))

	### BAR
	barStars, barAges = generateStars(length=6, n=20000, \
		width=lambda x: np.sqrt(1-(x-0.5)**2*4)*0.6, height=lambda x: np.sqrt(1-(x-0.5)**2*4)*0.6, k=1000, ageK=0.1)
	x = (barStars[:,0]-3) * np.cos(60*d2r) + barStars[:,1] * np.sin(60*d2r)
	y = (barStars[:,0]-3) * np.sin(60*d2r) - barStars[:,1] * np.cos(60*d2r)
	barStars[:,0] = x
	barStars[:,1] = y
	galaxy.append((barStars, barAges*0.8))
	
	return np.vstack([a[0] for a in galaxy]), np.concatenate([a[1] for a in galaxy])


def generateGas(stars, ages, p0=(-5,-5,-2), p1=[+5,+5,+2], spacing=(0.5, 0.5, 0.5)):
	x, y, z = np.meshgrid(np.arange(p0[0], p1[0]+0.01, spacing[0]), np.arange(p0[1], p1[1]+0.01, spacing[1]), np.arange(p0[2], p1[2]+0.01, spacing[2]), indexing='ij')
	xyz = np.stack((x,y,z), axis=-1)
	density = np.zeros_like(x)#np.exp(-(x**2 + y**2 + z**2)/2)
	### add volume
	for s, a in zip(stars, ages):
		density += np.exp(-((xyz-s)**2).sum(axis=-1) /2/0.4**2) * (2-a)
	return density


def bgTexture(plotter):
	def textureOUV(imgSize=(10,10), gc=(0,0), sun=(0,1), dSun=8.15):
		### calculate physical position of lower-left (O) / lower-right (U) / upper-left (V) corner of the image according to GC and Sun pixel.
		L = np.sqrt( (sun[0]-gc[0])**2 + (sun[1]-gc[1])**2 )
		scale = dSun / L
		unitY = ( (sun[0]-gc[0])/L, (sun[1]-gc[1])/L )
		unitX = ( -(sun[1]-gc[1])/L, (sun[0]-gc[0])/L )

		_o = [1-gc[0], imgSize[1]-gc[1]]
		O = (np.dot(_o, unitX)*scale,  np.dot(_o, unitY)*scale, 0)
		_u = [imgSize[0]-gc[0], imgSize[1]-gc[1]]
		U = (np.dot(_u, unitX)*scale,  np.dot(_u, unitY)*scale, 0)
		_v = [1-gc[0], 1-gc[1]]
		V = (np.dot(_v, unitX)*scale,  np.dot(_v, unitY)*scale, 0)

		return O, U, V
	X, Y = np.meshgrid(np.linspace(-30, 30, 10), np.linspace(-30, 30, 10))
	Z = X*0-1
	plane = pv.StructuredGrid(X, Y, Z)
	if pic==0:
		texture = pv.read_texture('art_Hurt_ext.jpg')
		O,U,V = textureOUV(texture.dimensions, gc=(1500, 1660), sun=(1500, 2280))
	else:
		texture = pv.read_texture('art_Zheng.png')
		O,U,V = textureOUV(texture.dimensions, gc=(1270, 1210), sun=(1266, 756))
	plane.texture_map_to_plane(origin=O, point_u=U, point_v=V, inplace=True)
	return plane, texture






p = pv.Plotter(window_size=(1600,1600))#
p.set_background('black')


#stars, ages = generateRegions(nStars=4000, nRegions=30, length=100, width=0.5, height=0.5, k=1.01)
#stars, ages = generateArm(3,0,-8, n=4000, length=100, width=0.1, height=0.1, k=1.01)


# Stars
stars, ages = generateGalaxy()
print('Star Range: Max', stars.max(axis=0), 'Min', stars.min(axis=0))

starPoints = pv.PointSet(stars)
starPoints['color'] = ages
p.add_points(starPoints, scalars='color', style='points', opacity=1, cmap=cmapStar()) # points_gaussian


# Gas
p0 = (-16, -16, -5)
p1 = (+16, +22, +5)
spacing = (0.3, 0.3, 0.3)
if 0:
	density = generateGas(stars, ages, p0=p0, p1=p1, spacing=spacing)
	np.save('volume.npy', density)
else:
	density = np.load('volume.npy')
print(density.max())
#density = np.sqrt(density)

gasGrid = pv.UniformGrid()
gasGrid.dimensions = np.array(density.shape)+1
gasGrid.spacing = spacing
gasGrid.origin = [p-s/2 for p,s in zip(p0, spacing)]
gasGrid.cell_data["values"] = np.log(density.ravel(order="F")+1)

opacity = [0, 0, 0, 0.1, 0.3, 0.6, 1]
#p.add_volume(gasGrid, cmap=cmapCloud())#, opacity=opacity)


### pic of background
#plane, texture = bgTexture(p, pic=0)
#p.add_mesh(plane, texture=texture)

sun = pv.Sphere(radius=0.2, center=(0, 8.15, -1), theta_resolution=30, phi_resolution=30)
p.add_mesh(sun, color=(255, 255, 0))
cyg = pv.Sphere(radius=0.2, center=(2, 8.15, -1), theta_resolution=30, phi_resolution=30)
p.add_mesh(cyg, color=(255, 0, 0))

eyeD = 50
eyeAz = np.pi
eyeEl = np.pi/2
p.camera_position = [[eyeD*np.sin(eyeAz)*np.cos(eyeEl), eyeD*np.cos(eyeAz)*np.cos(eyeEl), eyeD*np.sin(eyeEl)+2], [0, 0, 2], [0, 0, 1]]
#p.camera.zoom(1)

p.show()
