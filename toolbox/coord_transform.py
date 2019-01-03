import astropy.coordinates as coord
import astropy.units as u
import astropy.time
from datetime import datetime as dt
from dateutil.parser import parse


class transformer():
	def __init__(self, **params):
		self.cenpa_latlon = {"lon":-122.302753, "lat":47.660530}
		self.l = coord.Galactocentric().galcen_coord.ra #RA of galactic center
		self.b = coord.Galactocentric().galcen_coord.dec #DEC of galactic center
		self.dist = coord.Galactocentric().galcen_distance #Distance to galactic center from solar system barycenter
		self.con_ICRS =  coord.ICRS(ra=self.l, dec = self.b, distance = self.dist) #put galactic center coordinates into container
		self.con_GAL = self.con_ICRS.transform_to(coord.Galactocentric)
		self.solar_vector = self.con_GAL.galcen_v_sun
		self.solar_vector_mag = self.solar_vector.norm()
		
	def experiment_vel_galacticframe(self, timestamp):
		"""
		Method finds the position and velocity of the experimental apparatus in the ITRS coordinate system (geocentric) at the given timestamp
		"""
		cenpa_location_xyz = coord.EarthLocation.from_geodetic(lon=self.cenpa_latlon["lon"], lat = self.cenpa_latlon["lat"])
		time = astropy.time.Time(parse(timestamp))
		cenpa_ITRS = coord.ITRS(cenpa_location_xyz.x, cenpa_location_xyz.y, cenpa_location_xyz.z, representation_type = coord.CartesianRepresentation, v_x=0*(u.km/u.s), v_y=0*(u.km/u.s), v_z=0*(u.km/u.s), differential_type = coord.CartesianDifferential, obstime=time)
		
		cenpa_Galactic = cenpa_ITRS.transform_to(coord.Galactic)
		cenpa_velocity = cenpa_Galactic.velocity
		gvel_sun = coord.Galactocentric().galcen_v_sun
		gvel_GalacticIntertialFrame = gvel_sun + cenpa_velocity
		
		return gvel_GalacticIntertialFrame
		
		
	def earth_posvel_ICRS(self, timestamp):
		"""
		Method finds the position and velocity of the earth in the ICRS coordinate system (barycenter of Solar System) at the given timestamp.
		"""
		time = astropy.time.Time(parse(timestamp))
		
		position_xyz, velocity_xyz = coord.solar_system._get_body_barycentric_posvel('earth', time, get_velocity=True)
		position_xyz = position_xyz*(u.AU.to(u.km))*(u.km/u.AU)
		velocity_xyz = velocity_xyz*(u.AU.to(u.km)/(u.d.to(u.s)))*(u.km/u.s)*(u.d/u.AU)
		return position_xyz, velocity_xyz
		
	def earth_vel_galacticframe(self, timestamp):
		"""
		Method returns the velocity of the earth in the frame of the galaxy at the given timestamp
		"""
		time = astropy.time.Time(parse(timestamp))
		pos, vel = self.earth_posvel_ICRS(timestamp)
		
		gvel = coord.ICRS(x = pos.x.value*u.km, y = pos.y.value*u.km, z = pos.z.value*u.km,v_x = vel.x.value*u.km/u.s, v_y = vel.y.value*u.km/u.s, v_z = vel.z.value*u.km/u.s,representation_type = coord.CartesianRepresentation,differential_type = coord.CartesianDifferential).transform_to(coord.Galactic).velocity
		gvel_sun = coord.Galactocentric().galcen_v_sun
		gvel_GalacticIntertialFrame = gvel_sun + vel
		
		return gvel_GalacticIntertialFrame

cl = transformer()
st = "2019-1-2 13:48:43.644409"
vel = cl.earth_posvel_ICRS(st)
print(vel)












