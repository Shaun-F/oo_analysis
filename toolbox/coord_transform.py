import astropy.coordinates as coord
import astropy.units as u
import astropy.time
from datetime import datetime as dt
from dateutil.parser import parse
import numpy as np


class transformer():
	def __init__(self, **params):
		self.cenpa_latlon = {"lon":-122.302753, "lat":47.660530}
		
	def experiment_vel_GalacticFrame(self, timestamp):
		"""
		Method finds the position and velocity of the experimental apparatus in the ITRS coordinate system (geocentric) at the given timestamp
		"""
		cenpa_location_xyz = coord.EarthLocation.from_geodetic(lon=self.cenpa_latlon["lon"], lat = self.cenpa_latlon["lat"])
		time = astropy.time.Time(parse(timestamp))
		
		cenpa_ITRS = coord.ITRS(cenpa_location_xyz.x, cenpa_location_xyz.y, cenpa_location_xyz.z, representation_type = coord.CartesianRepresentation, v_x=0*(u.km/u.s), v_y=0*(u.km/u.s), v_z=0*(u.km/u.s), differential_type = coord.CartesianDifferential, obstime=time)
		
		cenpa_Galactic_centered = cenpa_ITRS.transform_to(coord.Galactocentric)
		
		#get tangential velocity
		x = cenpa_Galactic_centered.x.value
		y = cenpa_Galactic_centered.y.value
		vx = cenpa_Galactic_centered.v_x.value
		vy = cenpa_Galactic_centered.v_y.value
		
		cenpa_tangential_velocity = (y*vx - x*vy)/(np.sqrt(x**2 + y**2)) + np.arctan(y/x)*(x*vx + y*vy)/(np.sqrt(x**2 + y**2))
		return cenpa_tangential_velocity #Tangential velocity to galactic center
		
		
	def earth_posvel_ICRS(self, timestamp):
		"""
		Method finds the position and velocity of the earth in the ICRS coordinate system (barycenter of Solar System) at the given timestamp.
		"""
		time = astropy.time.Time(parse(timestamp))
		
		position_xyz, velocity_xyz = coord.solar_system._get_body_barycentric_posvel('earth', time, get_velocity=True)
		position_xyz = position_xyz*(u.AU.to(u.km))*(u.km/u.AU)
		velocity_xyz = velocity_xyz*(u.AU.to(u.km)/(u.d.to(u.s)))*(u.km/u.s)*(u.d/u.AU)
		return position_xyz, velocity_xyz
		
	def earth_vel_GalacticFrame(self, timestamp):
		"""
		Method returns the velocity of the earth in the frame of the galaxy at the given timestamp
		"""
		time = astropy.time.Time(parse(timestamp))
		pos, vel = self.earth_posvel_ICRS(timestamp)
		
		gvel = coord.ICRS(x = pos.x.value*u.km, y = pos.y.value*u.km, z = pos.z.value*u.km,v_x = vel.x.value*u.km/u.s, v_y = vel.y.value*u.km/u.s, v_z = vel.z.value*u.km/u.s,representation_type = coord.CartesianRepresentation,differential_type = coord.CartesianDifferential).transform_to(coord.Galactocentric)
		
		x = gvel.x.value
		y = gvel.y.value
		vx = gvel.v_x.value
		vy = gvel.v_y.value
		
		earth_tangential_velocity = (y*vx - x*vy)/(np.sqrt(x**2 + y**2)) + np.arctan(y/x)*(x*vx + y*vy)/(np.sqrt(x**2 + y**2))
		return earth_tangential_velocity
		
	def solar_vel_GalacticFrame(self):
		return coord.Galactocentric().galcen_v_sun.norm().value











