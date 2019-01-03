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
		
		cenpa_Galactic = cenpa_ITRS.transform_to(coord.Galactic)
		cenpa_velocity = cenpa_Galactic.velocity
		gvel_sun = coord.Galactocentric().galcen_v_sun
		
		cenpa_vel_delta_vectored = (cenpa_velocity + gvel_sun)
		cenpa_vel_delta = cenpa_vel_delta_vectored.norm() - gvel_sun.norm()
		
		return cenpa_vel_delta.value
		
		
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
		
		gvel = coord.ICRS(x = pos.x.value*u.km, y = pos.y.value*u.km, z = pos.z.value*u.km,v_x = vel.x.value*u.km/u.s, v_y = vel.y.value*u.km/u.s, v_z = vel.z.value*u.km/u.s,representation_type = coord.CartesianRepresentation,differential_type = coord.CartesianDifferential).transform_to(coord.Galactic).velocity
		gvel_sun = coord.Galactocentric().galcen_v_sun
		
		gvel_GalacticIntertialFrame = gvel_sun + vel
		gvel_GalacticIntertialFrame_delta = gvel_GalacticIntertialFrame.norm() - gvel_sun.norm()
		
		return gvel_GalacticIntertialFrame_delta.value
		
	def solar_vel_GalacticFrame(self):
		return coord.Galactocentric().galcen_v_sun.norm().value











