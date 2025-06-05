import pyvista as pv
import numpy as np
import matplotlib.image as mpimg
from shared import *
import matplotlib.pyplot as plt
from mcmcFittingFunction0129 import function


artist = 0
pWarp = p_1comp
zScale = 8

### Plotter
plotter = pv.Plotter(window_size=(2400,1600))

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


def arm_wave(R, Phi, arm='out', width=0.7, rScale=1, spline=0.01):
	### return a smoothed spline function of outer arm wave

	### load clouds
	if arm=='per':
		data = np.loadtxt('per_para.txt',comments='#')
		data = data[data[:,19]<76.6]
		p = best_per
	elif arm=='out':
		data = np.loadtxt('out_para.txt',comments='#')
		p = best_out
	else:
		data = np.loadtxt('osc2_para.txt',comments='#')
		p = best_osc
	vv = data[:,2]
	rr = data[:,11]
	zz = data[:,14]
	az = data[:,19]
	mass = data[:,7]/data[:,8]

	### get residual
	rr_on_arm = function_arm(az, p)
	zz_res = zz - function_warp( (rr_on_arm, az), p=pWarp)
	### bin
	azcen, azrms, zcen_res, zrms_res = cal_zcen_zrms(az, zz_res, weights=mass, binsize=4, nbin=60)
	idx = np.isfinite(zcen_res)
	### add zero edge
	x = [np.max(azcen[idx])+8, np.max(azcen[idx])+4, *azcen[idx], np.min(azcen[idx])-4, np.min(azcen[idx])-8][::-1]
	y = [0, 0, *zcen_res[idx], 0, 0][::-1]

	### interpolate
	#from scipy.interpolate import CubicSpline
	#wave = CubicSpline(azcen[idx][::-1], zcen_res[idx][::-1])

	#from scipy.interpolate import splrep, splev
	#wave = splrep(azcen[idx][::-1], zcen_res[idx][::-1], s=0.03)

	from scipy.interpolate import UnivariateSpline
	wave = UnivariateSpline(x, y)
	wave.set_smoothing_factor(spline)

	if 0:
		### check interpolated result
		xx = np.linspace(min(x), max(x), 500)
		plt.plot(x, y, 's')
		plt.plot(xx, wave(xx))
		plt.show()

	### get r of arm on grid
	phi_arm = Phi[:,0]
	r_arm = function_arm(phi_arm, p)
	z_res_arm = wave(phi_arm)
	z_res_arm[phi_arm<min(x)] = 0
	z_res_arm[phi_arm>max(x)] = 0
	Z_res = np.zeros(R.shape)
	for i in range(len(phi_arm)):
		Z_res[i] += z_res_arm[i] * np.exp(-(R[i]-r_arm[i]*rScale)**2 / 2 / width**2)

	return Z_res


### Create mesh
def createPlane(r = [8,9], phi = [34.5, 40]):
	if r.ndim==1:
		R, Phi = np.meshgrid(r, phi)
	else:
		R, Phi = r, phi
	X = R*np.sin(Phi/180*np.pi)
	Y = R*np.cos(Phi/180*np.pi)

	### warp component
	Z = function_warp((R, Phi), p=pWarp) * zScale

	### GC
	#Z += np.exp(-R**2/2/1**2) *0.3

	### arm component
	Z += arm_wave(R, Phi, arm='per', width=0.5, rScale=0.95 if artist else 0.9, spline=0.01) * zScale
	Z += arm_wave(R, Phi, arm='out', width=0.7, rScale=1.0 if artist else 0.9, spline=0.01) * zScale
	Z += arm_wave(R, Phi, arm='osc', width=1.1, rScale=1.1 if artist else 1.04, spline=1.2) * zScale
	return X, Y, Z, R, Phi


def sigmoid(x, x0=0, width=1):
	return 1/(1+np.exp(-(x-x0)/width)) 

### Create a structured grid
X, Y, Z, R, Phi = createPlane(r = np.linspace(0, 31, 105*2), phi = np.linspace(-90, 270, 360*2))
plane = pv.StructuredGrid(X, Y, Z)

### texture
if artist:
	texture = pv.read_texture('art_Hurt_ext.jpg')
	O,U,V = textureOUV(texture.dimensions, gc=(1500, 1660), sun=(1500, 2280))
	plane.texture_map_to_plane(origin=O, point_u=U, point_v=V, inplace=True)
else:
	if 0:
		texture = pv.read_texture('art_Zheng.jpg')
		O,U,V = textureOUV(texture.dimensions, gc=(1266, 1206), sun=(1266, 756))
		plane.texture_map_to_plane(origin=O, point_u=U, point_v=V, inplace=True)
	else:
		'''
		### transparent texture
		from PIL import Image
		img = Image.open('art_Zheng_HiRes.jpg').convert('RGB')
		img_np = np.array(img)
		#alpha = np.array(img.convert("L"))
		#alpha = (sigmoid(alpha.astype(float), 40, 5)*220).astype(np.uint8)
		#alpha = np.where(np.all(img_np == [0, 0, 0], axis=-1), 0, 255).astype(np.uint8)
		alpha = np.zeros(img_np[...,0].shape, dtype=np.uint8)+255
		#y,x = np.meshgrid(np.linspace(-150,150,301), np.linspace(-150,150,301))
		alpha[2935-150:2935+151, 4042-150:4042+151] = 100#(255*np.exp(-(x**2+y**2)/2/30**2)).astype(np.uint8)	#carve a hole near sun
		#print(alpha.dtype)
		#plt.imshow(alpha)
		#plt.show()
		img_rgba = np.dstack((img_np, alpha))
		texture = pv.Texture(img_rgba)
		'''
		texture = pv.read_texture('art_Zheng_HiRes.jpg')
		O,U,V = textureOUV(texture.dimensions, gc=(4042, 3976), sun=(4042, 2935))
		plane.texture_map_to_plane(origin=O, point_u=U, point_v=V, inplace=True)

plotter.add_mesh(plane, texture=texture)


c=0
### ARM cutaway
phi = np.linspace(-50, 170, 220)
r = function_arm(phi, p=best_out) * (1.0 if artist else 0.9)
Phi = phi[:,np.newaxis] * np.ones(200)
R = r[:, np.newaxis] + np.linspace(-5, 5, 200)
delta = np.zeros(R.shape) + np.linspace(-5, 5, 200)
X, Y ,Z, R, Phi = createPlane(r=R, phi=Phi)
armCut = pv.StructuredGrid(X, Y, Z+0.1)
### half transparent
alpha = np.exp(-delta**2/0.2**2) * (sigmoid(Phi, -23, 2) - sigmoid(Phi, 150, 2))
color = np.zeros((X.size, 4), dtype=np.uint8)
color[:, 0:3] = [176,196,222] if c else [220,200,180]
color[:, 3] = alpha.T.ravel() * 150#120
armCut.point_data['RGBA'] = color
plotter.add_mesh(armCut, scalars='RGBA', rgba=True)

'''
### WARP cutaway
X, Y, Z, R, Phi = createPlane(r = np.linspace(4, 15, 100), phi = np.linspace(89.3, 89.3+30, 30))
warpCut = pv.StructuredGrid(X, Y, Z+0.1)
### half transparent
alpha = np.exp(-(R-9)**2/4.5**2 -(Phi-89.3)**2/5**2)
color = np.zeros((X.size, 4), dtype=np.uint8)
color[:, 0:3] = [160, 200, 250]
color[:, 3] = alpha.T.ravel() * 190
warpCut.point_data['RGBA'] = color
plotter.add_mesh(warpCut, scalars='RGBA', rgba=True)
'''

### CORRUGATION cutaway
X, Y, Z, R, Phi = createPlane(r = np.linspace(4, 18, 100), phi = np.linspace(34.5, 34.5+30, 100))
sinCut = pv.StructuredGrid(X, Y, Z+0.05)
### half transparent
alpha = np.exp(-(R*(Phi-34.5))**2/80**2) * (sigmoid(R, 8, 0.5) - sigmoid(R, 16, 0.5))
color = np.zeros((X.size, 4), dtype=np.uint8)
color[:, 0:3] = [120,240,200] if c else [176,196,222] # [160, 200, 250]
color[:, 3] = alpha.T.ravel() * 140#120
sinCut.point_data['RGBA'] = color
plotter.add_mesh(sinCut, scalars='RGBA', rgba=True, specular=1.0, specular_power=100, smooth_shading=0)



# line from gc to sun to verify texture coordiante
#gc2sun = pv.Line((0,0,0), (0,8.15,0))
#plotter.add_mesh(gc2sun, line_width=20, color='red')

### Sun
sun_core = pv.Sphere(radius=0.2, center=(0, 8.15, 0), theta_resolution=30, phi_resolution=30)
plotter.add_mesh(sun_core, color=(255, 255, 0))
#sun_glow = pv.Sphere(radius=0.3, center=(0, 8.15, 0), theta_resolution=30, phi_resolution=30)
#plotter.add_mesh(sun_glow, color=(255, 255, 255), opacity=0.3, smooth_shading=True)


### annotation
plotter.add_lines(np.array([[8,-4.2,0], [8,-4.2,8]]), color='w', width=3)
plotter.add_point_labels([[8,-4.2,8]], ['Outer Arm\nCorrugation'], font_size=40, text_color='w', shape=False)
plotter.add_lines(np.array([[14,0,0], [14,0,8]]), color='w', width=3)
plotter.add_point_labels([[14,0,8]], ['Warp'], font_size=40, text_color='w', shape=False)
plotter.add_lines(np.array([[8.5,12,0], [8.5,12,8]]), color='w', width=3)
plotter.add_point_labels([[8.5,12,8]], ['Bar-direction\nCorrugation'], font_size=40, text_color='w', shape=False)

# setting
plotter.set_scale(xscale=1, yscale=1, zscale=1)   #scale z by 10
eyeD = 200
eyeAz = 265 /180*np.pi
eyeEl = 22 /180*np.pi
'''
plotter.set_position([np.sin(eyeAz)*np.cos(eyeEl), np.cos(eyeAz)*np.cos(eyeEl), np.sin(eyeEl)])   #position of eye
plotter.set_focus([0, 0, 0])  #position to look at
print(help(plotter.set_focus))
plotter.set_viewup([0, 0, 1]) #z is up
'''
plotter.camera_position = [[eyeD*np.sin(eyeAz)*np.cos(eyeEl), eyeD*np.cos(eyeAz)*np.cos(eyeEl), eyeD*np.sin(eyeEl)+2], [0, 0, 2], [0, 0, 1]]
plotter.camera.zoom(4)

plotter.set_background('black')
'''
plotter.show_bounds(bounds=[-30, 30, -30, 30, -5, 5], grid=True, bold=True, \
	font_size=10, color='gray', location='outer', padding=0, \
	xlabel='X (kpc)', \
	ylabel='Y (kpc)', \
	zlabel='Z - model (kpc)')
'''
plotter.show()
#plotter.save_graphic('fig/schematic.eps')
plotter.screenshot('fig/schematic.png')


