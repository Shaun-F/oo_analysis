import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
import astropy.units as u
import astropy.time
from datetime import datetime as dt
from dateutil.parser import parse
import numpy


"""
See http://docs.astropy.org/en/stable/api/astropy.coordinates.Galactocentric.html for more information on 
"""

class transformer():
	def __init__(self, **params):
		self.cenpa_latlon = {"lon":-122.302753, "lat":47.660530}
		
	def experiment_vel_GalacticFrame(self, timestamp):
		"""
		Method finds the position and velocity of the experimental apparatus in the ITRS coordinate system (geocentric) at the given timestamp
		"""
		cenpa_location_xyz = coord.EarthLocation.from_geodetic(lon=self.cenpa_latlon["lon"], lat = self.cenpa_latlon["lat"])
		times = self.parse_timestamps(timestamp)
		
		if isinstance(timestamp, numpy.ndarray) or isinstance(timestamp, list):
			x_values = [cenpa_location_xyz.x for iter in range(len(timestamp))]
			y_values = [cenpa_location_xyz.y for iter in range(len(timestamp))]
			z_values = [cenpa_location_xyz.z for iter in range(len(timestamp))]
			v_x_values = numpy.zeros(len(timestamp))*u.km/u.s
			v_y_values = numpy.zeros(len(timestamp))*u.km/u.s
			v_z_values = numpy.zeros(len(timestamp))*u.km/u.s
			#transform to ITRS coordinate based object then assign velocity differentials
			cenpa_ITRS = SkyCoord(x_values, y_values, z_values, frame='itrs', 
								v_x = v_x_values, v_y = v_y_values, v_z=v_z_values,
								obstime = times)
	
			cenpa_Galactic_centered = cenpa_ITRS.transform_to(coord.Galactocentric)  #Current bottleneck
		
			#returns dictionary of vectors
			return {timestamp[i]: [cenpa_Galactic_centered.v_x.value[i], cenpa_Galactic_centered.v_y.value[i], cenpa_Galactic_centered.v_z.value[i]] for i in range(len(timestamp))}
		
		elif isinstance(timestamp, str):
			x_values = cenpa_location_xyz.x
			y_values = cenpa_location_xyz.y
			z_values = cenpa_location_xyz.z
			v_x_values = 0*u.km/u.s
			v_y_values = 0*u.km/u.s
			v_z_values = 0*u.km/u.s
			#transform to ITRS coordinate based object then assign velocity differentials
			cenpa_ITRS = SkyCoord(x_values, y_values, z_values, frame='itrs', 
									v_x = v_x_values, v_y = v_y_values, v_z=v_z_values,
									obstime = times)
		
			cenpa_Galactic_centered = cenpa_ITRS.transform_to(coord.Galactocentric)  #Current bottleneck
			
			#returns vector
			return numpy.asarray([cenpa_Galactic_centered.v_x.value, cenpa_Galactic_centered.v_y.value, cenpa_Galactic_centered.v_z.value])
		
		elif isinstance(timestamp, dict):
			x_values = [cenpa_location_xyz.x for iter in range(len(timestamp))]
			y_values = [cenpa_location_xyz.y for iter in range(len(timestamp))]
			z_values = [cenpa_location_xyz.z for iter in range(len(timestamp))]
			v_x_values = numpy.zeros(len(timestamp))*u.km/u.s
			v_y_values = numpy.zeros(len(timestamp))*u.km/u.s
			v_z_values = numpy.zeros(len(timestamp))*u.km/u.s
			#transform to ITRS coordinate based object then assign velocity differentials
			cenpa_ITRS = SkyCoord(x_values, y_values, z_values, frame='itrs', 
								v_x = v_x_values, v_y = v_y_values, v_z=v_z_values,
								obstime = times)
	
			cenpa_Galactic_centered = cenpa_ITRS.transform_to(coord.Galactocentric)  #Current bottleneck
		
			#returns dictionary of vectors
			return {list(timestamp.values())[i]: [cenpa_Galactic_centered.v_x.value[i], cenpa_Galactic_centered.v_y.value[i], cenpa_Galactic_centered.v_z.value[i]] for i in range(len(timestamp))}
			
		
		
		#coord.ITRS(cenpa_location_xyz.x, cenpa_location_xyz.y, cenpa_location_xyz.z, representation_type = coord.CartesianRepresentation, v_x=0*(u.km/u.s), v_y=0*(u.km/u.s), v_z=0*(u.km/u.s), differential_type = coord.CartesianDifferential, obstime=time)
		
	def experiment_posvel_GCRS(self, timestamp):
		"""
		Method finds the position and velocity of the experimental apparatus in the GCRS coordinate system (earth centered) at the given timestamp.
		Typically of the order ~0.3 km/s (about 1/100 the orbital velocity of earth)
		"""
		cenpa_location_xyz = coord.EarthLocation.from_geodetic(lon=self.cenpa_latlon["lon"], lat = self.cenpa_latlon["lat"])
		times = self.parse_timestamps(timestamp)
		
		posvel = cenpa_location_xyz.get_gcrs_posvel(obstime=times)
		if isinstance(timestamp, dict):
			return {list(timestamp.values())[i]: [posvel[1].x[i].to(u.km/u.s).value, posvel[1].y[i].to(u.km/u.s).value, posvel[1].z[i].to(u.km/u.s).value] for i in range(len(timestamp))}
		elif isinstance(timestamp, list) or isinstance(timestamp. numpy.ndarray):
			return {timestamp[i]: [posvel[1].x[i].to(u.km/u.s).value, posvel[1].y[i].to(u.km/u.s).value, posvel[1].z[i].to(u.km/u.s).value] for i in range(len(timestamp))}
		else:
			return {timestamp: [posvel[1].x.to(u.km/u.s).value, posvel[1].y.to(u.km/u.s).value, posvel[1].z.to(u.km/u.s).value]}
	
	def earth_posvel_ICRS(self, timestamp):
		"""
		Method finds the position and velocity of the earth in the ICRS coordinate system (barycenter of Solar System) at the given timestamp.
		"""
		time = self.parse_timestamps(timestamp)
		
		position_xyz, velocity_xyz = coord.solar_system._get_body_barycentric_posvel('earth', time, get_velocity=True)
		position_xyz = position_xyz*(u.AU.to(u.km))*(u.km/u.AU)
		velocity_xyz = velocity_xyz*(u.AU.to(u.km)/(u.d.to(u.s)))*(u.km/u.s)*(u.d/u.AU)
		return position_xyz, velocity_xyz

	def earth_vel_GalacticFrame(self, timestamp):
		"""
		Method returns the velocity of the earth in the frame of the galaxy at the given timestamp
		Typically of the order ~30 km/s
		"""
		
		time = self.parse_timestamps(timestamp)
		pos, vel = self.earth_posvel_ICRS(timestamp)
		
		gvel = coord.ICRS(x = pos.x.value*u.km, y = pos.y.value*u.km, z = pos.z.value*u.km,v_x = vel.x.value*u.km/u.s, v_y = vel.y.value*u.km/u.s, v_z = vel.z.value*u.km/u.s,representation_type = coord.CartesianRepresentation,differential_type = coord.CartesianDifferential).transform_to(coord.Galactocentric)
		
		if isinstance(timestamp, list) or isinstance(timestamp, numpy.ndarray):
			return {timestamp[i]: [gvel.v_x[i].to(u.km/u.s).value,gvel.v_y[i].to(u.km/u.s).value,gvel.v_z[i].to(u.km/u.s).value] for i in range(len(timestamp))}
		elif isinstance(timestamp, dict):
			return {list(timestamp.values())[i]: [gvel.v_x[i].to(u.km/u.s).value,gvel.v_y[i].to(u.km/u.s).value,gvel.v_z[i].to(u.km/u.s).value] for i in range(len(timestamp))}
		else:
			return {timestamp: [gvel.v_x.to(u.km/u.s).value,gvel.v_y.to(u.km/u.s).value,gvel.v_z.to(u.km/u.s).value]}
		
	def solar_vel_GalacticFrame(self, timestamp=None):
		galcen_v_sun = coord.Galactocentric().galcen_v_sun
		
		if isinstance(timestamp, list) or isinstance(timestamp, numpy.ndarray):
			return {time: numpy.asarray([galcen_v_sun.d_x.value, galcen_v_sun.d_y.value, galcen_v_sun.d_z.value]) for time in timestamp}
		
		elif isinstance(timestamp, dict):
			return {time:numpy.asarray([galcen_v_sun.d_x.value, galcen_v_sun.d_y.value, galcen_v_sun.d_z.value]) for time in list(timestamp.values())}
		
		if timestamp==None:
			return numpy.asarray([galcen_v_sun.d_x.value, galcen_v_sun.d_y.value, galcen_v_sun.d_z.value])
		
	def parse_timestamps(self, timestamp):
		"""
		Converts input timestamps into astropy time objects. Input type can be string, list, numpy array, or dictionary
		"""
		if isinstance(timestamp, str):	
			times = astropy.time.Time(numpy.asarray([parse(timestamp, dayfirst=True)])) 
		
		elif isinstance(timestamp, list) or isinstance(timestamp, numpy.ndarray):
			times = astropy.time.Time(numpy.asarray([parse(x, dayfirst=True) for x in timestamp]))
		
		elif isinstance(timestamp, dict):
			times = astropy.time.Time(numpy.asarray([parse(x, dayfirst=True) for x in timestamp.values()]))
		
		return times










