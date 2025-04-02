import numpy as np
#from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5 
from astropy.coordinates import Angle, Latitude, Longitude

class KinUnit(dict):
	"""

	"""
	def __init__(self, Ro = 8.34, To = 240.0, dTdR = -0.2,
	                   Uo = 10.7, Vo = 15.6 , Wo   = 8.9 ,
	                   Us = 2.9 , Vs = -1.6 , Ws   = -3.0,
	                   err_v_lsr = 5.0):
		self['Ro']         = Ro                   #
		self['To']         = To                   #
		self['dTdR']       = dTdR                 #
		self['Uo']         = Uo                   #
		self['Vo']         = Vo                   #
		self['Wo']         = Wo                   #
		self['Us']         = Us                   #
		self['Vs']         = Vs                   #
		self['Ws']         = Ws                   #
		self['err_v_lsr']  = err_v_lsr            #
	def replear(self, unit):
		if not isinstance(unit,KinUnit) :raise TypeError(unit,type(unit))
		self['Ro']         = unit['Ro']
		self['To']         = unit['To']
		self['dTdR']       = unit['dTdR']
		self['Uo']         = unit['Uo']
		self['Vo']         = unit['Vo']
		self['Wo']         = unit['Wo']
		self['Us']         = unit['Us']
		self['Vs']         = unit['Vs']
		self['Ws']         = unit['Ws']
		self['err_v_lsr']  = unit['err_v_lsr']

A5 = KinUnit( Ro = 8.15, To = 236  , dTdR = -0.2, Uo = 10.6, Vo = 10.7, Wo   = 7.6 ,Us = 6.0 , Vs = -4.3 , Ws   = -3.0,err_v_lsr = 5.0)


class Distance(object):
	"""




	"""
	def __init__(self, skyc = None, v_lsr = None, farnear = 0, unit = KinUnit()):
		"""


		"""
		if not isinstance(skyc,SkyCoord)      : raise TypeError(skyc,type(skyc))
		skyc    = skyc.reshape(skyc.size)
		v_lsr   = np.array([v_lsr])
		v_lsr   = v_lsr.reshape([v_lsr.size])
		farnear = np.array([farnear])
		farnear = farnear.reshape([farnear.size])
		if len(skyc)  != len(v_lsr)   : raise ValueError("len(skyc) != len(v_lsr)")
		if len(skyc)  != len(farnear) : raise ValueError("len(skyc) != len(farnear)")
		if len(v_lsr) != len(farnear) : raise ValueError("len(v_lsr) != len(farnear)")
		self.skyc       = skyc
		self.v_lsr      = v_lsr
		self.farnear    = farnear
		self.len        = len(skyc)
		self.Ro         = unit['Ro']
		self.To         = unit['To']
		self.dTdR       = unit['dTdR']
		self.Uo         = unit['Uo']
		self.Vo         = unit['Vo']
		self.Wo         = unit['Wo']
		self.Us         = unit['Us']
		self.Vs         = unit['Vs']
		self.Ws         = unit['Ws']
		self.err_v_lsr  = unit['err_v_lsr']

	def Kdist(self):
		"""
		Calculates kinematic distances using the revised prescription given in 
		Reid et al 2009 (Paper VI):  
		Version 2 (projects Vlsr to Galactic Plane)
		Version 3 (fixes bug: used Theta instead of Theta_0 in rootterm)
		return : vlsr,distance,error_high,error_low
		"""
#		Convert coordinates to radians and then to hours, deg
		ra_rad   = self.skyc.fk5.ra.rad                              #call hmsrad ( ra_hhmmss, ra_rad )
		dec_rad  = self.skyc.fk5.dec.rad                             #call dmsrad (dec_ddmmss,dec_rad ) 
		ra       = self.skyc.fk5.ra/15                               #ra = ra_rad * 12  / np.pi
		dec      = self.skyc.fk5.dec                                 #dec=dec_rad * 180 / np.pi
#		Get Galactic coordinates of source
		ell,bee  = self.skyc.galactic.l.deg,self.skyc.galactic.b.deg #call radec_to_galactic ( ra, dec,  ell, bee )
#		Calculate "adjusted" Vlsr, including effects of (Us,Vs).
#		Since applying (Us,Vs) needs Dk, so do both at same time.
		v_lsr_rev, Dk_best = self.calc_Dk(self.v_lsr)
#		-------------------------------------------------------
#		Estimation of (asymmetric) uncertainties based only on
#		the uncertainty in Vlsr supplied in "parameter_file.inp".
#		Assumes Galactic & Source parameters have no uncertainties.

#		Check if v_lsr exceeds tangent point speed.  If so,
#		use tangent point speed as center for velocity offsets
#		used for error estimate.

#		Calculate tangent point speed (T_tp)
		sin_l = np.sin( np.deg2rad(ell) )
		R_tp = self.Ro * sin_l                             # kpc
		T_tp = self.To + self.dTdR * (R_tp - self.Ro)      # km/s

#		Check of v_lsr_rev exceeds tangent point speed
		del_v = np.zeros([self.len])+0.
		v_tp  = np.zeros([self.len])-9999.9
#		if ( ell < 90.  ) :                     # quadrant 1
#			v_tp = +T_tp - self.To*sin_l          # km/s 
#			if ( v_lsr_rev > v_tp ) : del_v = v_lsr_rev - v_tp
		temp = np.where(ell<90)[0]
		v_tp[temp] = +T_tp[temp] - self.To*sin_l[temp]
		temp2 = np.where(v_lsr_rev[temp] > v_tp[temp])[0]
		del_v[temp[temp2]] = v_lsr_rev[temp[temp2]] - v_tp[temp[temp2]]

#		if ( ell > 270. ) :                     # quadrant 4
#			v_tp = -T_tp - self.To*sin_l          # km/s 
#			if ( v_lsr_rev < v_tp ) :  del_v = v_lsr_rev - v_tp
		temp = np.where(ell > 270.)[0]
		v_tp[temp] = -T_tp[temp] - self.To*sin_l[temp]
		temp2 = np.where(v_lsr_rev[temp] > v_tp[temp])[0]
		del_v[temp[temp2]] = v_lsr_rev[temp[temp2]] - v_tp[temp[temp2]]
		
#		Adjust v_lsr (if needed) so as not to exceed tangent point
		v_lsr_adj = self.v_lsr - del_v

#		Get lower velocity distance uncertainty
		v_lsr_low = v_lsr_adj - self.err_v_lsr
		v_lsr_rev_low, D1 = self.calc_Dk(v_lsr_low)
#		Flag if close to tangent point velocity
#		if ( abs(v_lsr_rev_low - v_tp) < self.err_v_lsr ) : D1=Dk_best
		temp = np.where( abs(v_lsr_rev_low - v_tp) < self.err_v_lsr )[0]
		D1[temp] = Dk_best[temp]

#		Get higher velocity distance uncertainty
		v_lsr_high = v_lsr_adj + self.err_v_lsr
		v_lsr_rev_high, D2 = self.calc_Dk(v_lsr_high)
#		Flag if close to tangent point velocity
#		if ( abs(v_lsr_rev_high - v_tp) < self.err_v_lsr ) : D2=Dk_best
		temp = np.where( abs(v_lsr_rev_high - v_tp) < self.err_v_lsr )[0]
		D2[temp] = Dk_best[temp]

#		Calculate + (high) and - (low) errors
#		if ( D2 > D1 ) : 
#		     d_err_low  = D1  - Dk_best
#		     d_err_high = D2  - Dk_best
#		else :
#		     d_err_low  = D2  - Dk_best
#		     d_err_high = D1  - Dk_best
		d_err_low  = np.array([D1  - Dk_best,D2  - Dk_best]).min(0)
		d_err_high = np.array([D1  - Dk_best,D2  - Dk_best]).max(0)


#		If flagged distance (error=0), use other error estimate
#		if ( abs(d_err_low) == 0 ) : d_err_low  = -d_err_high
		temp = np.where( abs(d_err_low) == 0 )[0]
		d_err_low[temp]  = -d_err_high[temp]

#		if ( d_err_high     == 0 ) : d_err_high = -d_err_low
		temp = np.where( abs(d_err_high) == 0 )[0]
		d_err_high[temp]  = -d_err_low[temp]

#		Check for pathalogical cases...
#		if ( Dk_best < 0.0 ) : 
#		   Dk_best    = 0.0
#		   d_err_low  = 0.0
#		   d_err_high = 0.0
		temp = np.where( Dk_best < 0.0)[0]
		Dk_best[temp]    = 0.0
		d_err_low[temp]  = 0.0
		d_err_high[temp] = 0.0

#		if ( Dk_best+d_err_low < 0.0 ) : d_err_low = -Dk_best
		temp = np.where( Dk_best+d_err_low < 0.0 )[0]
		d_err_low[temp] = -Dk_best[temp]

		self.v_lsr_rev  = v_lsr_rev
		self.Dk_best    = Dk_best
		self.d_err_high = d_err_high
		self.d_err_low  = d_err_low
		return v_lsr_rev,Dk_best,d_err_high,d_err_low

	def calc_Dk(self,vlsr):
		"""
		Calculate revised Vlsr by converting standard Vlsr back to
		heliocentric, apply modern Solar Motion values (Uo,Vo,Wo), and 
		remove effects of average source non-ciruclar motion (Us,Vs,Ws).  
		Then calculate kinematic distance using the linear rotation 
		curve specified by Ro, To, and dTdR
		"""
		n_iter_max = 1000

#     ============================================================
#     Get source galactic coordinates...

		gal_long, gal_lat  = self.skyc.galactic.l.deg,self.skyc.galactic.b.deg

		gal_long_rad = np.deg2rad( gal_long )      # radians
		cos_l = np.cos( gal_long_rad )
		sin_l = np.sin( gal_long_rad )

		gal_lat_rad  = np.deg2rad( gal_lat )        # radians
		cos_b = np.cos( gal_lat_rad )
		sin_b = np.sin( gal_lat_rad )

#     -------------------------------------------------------------
#     Convert to true Heliocentric frame (v_rad)
#     Add back Sun's peculiar motion to radial velocity
#     Use old Standard Solar Motion (Sun moves at 20 km/s toward
#     18h, +30deg in 1900 coordinates), since this has been used to 
#     define V_lsr:

		Uo_IAU = 10.27            # km/s precessed to J2000
		Vo_IAU = 15.32
		Wo_IAU =  7.74

		v_helio = vlsr - (Vo_IAU*sin_l + Uo_IAU*cos_l)*cos_b - Wo_IAU*sin_b

#     -------------------------------------------------------------
#     Make "new" V(LSR) using best Solar Motion
#      eg, Hipparcos (Dehnen & Binney 1997) gives 
#		Uo = 10.00d0                # km/s
#		Vo =  5.25d0
#		Wo =  7.17d0
  
		v_newlsr = v_helio + (self.Vo*sin_l + self.Uo*cos_l)*cos_b + self.Wo*sin_b

#     --------------------------------------------------------
#     Remove effects of common peculiar motions specified in
#     "parameter_file.inp"

#     If dTdR.ne.0, need to know distance to get source Galactocentric 
#     radius (Rs) to evalute rotation curve.   So must iterate...
		n_iter =  0
		del_d  = np.zeros([self.len])+99.
		Dk     = np.zeros([self.len])+ 3.

		while ( (del_d > 0.01).all() and n_iter < n_iter_max ) :

#        Save old value of kinematic distance
			Dk_old = Dk

#        Calculate "gamma" angle and projected Galactocentric radius 
			d_proj = Dk * cos_b                     # kpc in Gal Plane
			r_sq   = self.Ro**2 + d_proj**2 - 2.*self.Ro * d_proj * cos_l
			r_proj = ( r_sq )**0.5                  # kpc in Gal Plane

#        Calculate Galactocentric longitude (beta in paper)...
			sin_beta =   d_proj  * sin_l     / r_proj
			cos_beta = ( self.Ro - d_proj*cos_l ) / r_proj
			beta     = np.arctan2( sin_beta, cos_beta )  # radians
			beta_deg = np.rad2deg( beta )                # deg

#        Calculate Sun-Maser-GC angle...
			gamma = np.pi - gal_long_rad - beta        # radians
			cos_gamma = np.cos( gamma )
			sin_gamma = np.sin( gamma )

			v_fixed = v_newlsr - (self.Vs*sin_gamma - self.Us*cos_gamma)*cos_b - self.Ws*sin_b  # km/s

#        -----------------------------------------------------------------
#        Calculate a kinematic distance using best Ro, To and dTdR

			Rs            = r_proj
			V_proj        = v_fixed * cos_b
			D_near, D_far = self.kinematic_distance ( V_proj, gal_long, self.Ro, self.To, r_proj, self.dTdR)

			Dk = D_near
#			if ( self.farnear != 0 ) : Dk = D_far
			Dk = np.zeros([self.len])+D_near
			temp = np.where(self.farnear!=0)[0]
			Dk[temp] = D_far[temp]

#        Ignore "farnear" flag if one of the values is zero
#			if ( D_near < 0 and D_far  > 0 ) : Dk = D_far
#			if ( D_far  < 0 and D_near > 0 ) : Dk = D_near
			temp = np.where((D_near < 0) * (D_far  > 0))
			Dk[temp] = D_far[temp]
			temp = np.where((D_near > 0) * (D_far  <0))
			Dk[temp] = D_near[temp]

			del_d = abs( Dk - Dk_old )
			n_iter = n_iter + 1

		v_lsr_rev = v_fixed

		return v_lsr_rev,Dk

	def kinematic_distance(self, Vlsr, gal_long, Ro, To, Rs, dTdR):
		"""
		Caluclate kinematic distance given Vlsr (projected in Gal plane)
		and information required to construct the kinematic model.  
		Returns both near and far distances (when appropriate)
		"""

		glongrad = np.deg2rad( gal_long )

		cos_l = np.cos(glongrad)
		sin_l = np.sin(glongrad)

		Rosinl = Ro * sin_l
		Rocosl = Ro * cos_l

#		Non-flat rotation curve allowed...
		Theta = To + dTdR*(Rs - Ro)
		Tosinl= To    * sin_l
		Tsinl = Theta * sin_l


#		Tsinl / ( Tosinl/Ro + Vlsr/Ro ) = Rs
		rootterm = Rocosl**2 + ( Tsinl / ( Tosinl/Ro + Vlsr/Ro ) )**2 - Ro**2
#		if ( rootterm < 0 ) : rootterm = 0
		rootterm*= (rootterm > 0)                        #test
		D_near   = np.zeros([self.len])             #test
		D_far    = np.zeros([self.len])             #test
#		if ( gal_long >= 0 and gal_long < 90 ) :
#			D_near = Rocosl - ( rootterm )**0.5
#			D_far  = Rocosl + ( rootterm )**0.5
		temp     = np.where( (gal_long >= 0) * (gal_long < 90) > 0)[0]
		D_near[temp] = Rocosl[temp] - ( rootterm[temp] )**0.5
		D_far[temp]  = Rocosl[temp] + ( rootterm[temp] )**0.5

#		if ( gal_long >= 90 and gal_long <= 270 ) :
#			D_near = Rocosl + ( rootterm )**0.5
#			D_far  = D_near
		temp     = np.where( (gal_long >= 90) * (gal_long <= 270) > 0)[0]
		D_near[temp] = Rocosl[temp] + ( rootterm[temp] )**0.5
		D_far[temp]  = D_near[temp]
		
#		if ( gal_long > 270 and gal_long < 360 ) :
#			D_near = Rocosl - ( rootterm )**0.5
#			D_far  = Rocosl + ( rootterm )**0.5
		temp     = np.where( (gal_long > 270) * (gal_long < 360) > 0)[0]
		D_near[temp] = Rocosl[temp] - ( rootterm[temp] )**0.5
		D_far[temp]  = Rocosl[temp] + ( rootterm[temp] )**0.5

		return D_near,D_far

def Kdpdf( distance,Dk_best,d_err_high,d_err_low, weigh =1):
	"""

	"""
	gauss = lambda x, a, mu, sigma: a*np.e**-(( x - mu )**2/(2 * sigma**2))
	pdf   = np.zeros([len(distance)], dtype = float)
	pdf[:np.where(distance>=Dk_best)[0][1]] = gauss(distance[:np.where(distance>=Dk_best)[0][1]], weigh, Dk_best, d_err_low)
	pdf[np.where(distance>=Dk_best)[0][1]:] = gauss(distance[np.where(distance>=Dk_best)[0][1]:], weigh, Dk_best, d_err_high)
	return pdf

