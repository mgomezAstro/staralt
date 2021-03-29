import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_moon, get_sun
import pytz
from datetime import datetime
import logging as _log


class ObjVisibility(object):
	'''
	A class to obtain the visibility of an astrophysical object in an specific date.

	...

	Attributes
	----------
	date-obs: str
		Date in which the observations will take place.
	location: str
		Observatory location. One of those listed in astropy.coordinates.EarthLocation.get_site_names().
	ra: float or list of floats
		Right Ascension (J2000, ICRS)
	dec: float or list of floats
		Declination (J2000, ICRS)
	names: str or list of str
		Name of the objects.
	saveFig: bool (optional; default: True)
		If True the program will save the figure.
	formatFig: str (optinal; default: pdf)
		Output format of the figure. One of the matplotlib supported formats.

	Methods
	-------
	staralt()
		Plots the staralt visibility plot of the list of objects.
	'''
	def __init__(self, date_obs, location, ra, dec, names, saveFig=True, formatFig="pdf"):
		'''
		Obtain all th enecessary data to construct the visibility plot of the desired objects at an specific date.

		Parameters
		----------
			date-obs: str
				Date in which the observations will take place.
			location: str
				Observatory location. One of those listed in astropy.coordinates.EarthLocation.get_site_names().
			ra: float or list of floats
				Right Ascension (J2000, ICRS)
			dec: float or list of floats
				Declination (J2000, ICRS)
			names: str or list of str
				Name of the objects.
			saveFig: bool (optional; default: True)
				If True the program will save the figure.
			formatFig: str (optinal; default: pdf)
				Output format of the figure. One of the matplotlib supported formats.

		Returns
		-------
			None.
		'''
		if not isinstance(names, list):
			self.names = [names]
		else:
			self.names = names
		if not isinstance(ra, list):
			ra = [ra]
			dec = [dec]
		self.coords = SkyCoord(ra, dec, unit="deg", frame="icrs")
		self.name_location = location
		self.location = EarthLocation.of_site(location)
		self.date_obs = Time(date_obs + " 00:00:00")
		self.time_zone = self.location.info.meta["timezone"]
		self.saveFig = saveFig
		self.formatFig = formatFig

	def __airmass_secz(self, degrees):
		'''
		Calculates the airmess assuming aplane-parallel atmosphere.

		Formula: sec(z)

		Parameters
		----------
			degrees: float
				Degrees from the zenithal position.

		Returns
		-------
			airmass: float
				The airmass in a determined angle.

		'''
		return 1. / np.cos(np.deg2rad(degrees))

	def staralt(self):
		'''
		Function to plot the staralt visibility map with moon position.

		Parameters
		----------
			None

		Returns
		-------
			Interactive plot if saveFigure=False or a file date_location_.formatFig.
		'''
		tz = pytz.timezone(self.time_zone)
		utcoffset = tz.utcoffset(datetime.now()).total_seconds() / 60. / 60. * u.hour
		print(f"UTC: {utcoffset}")

		# Plotting interval
		midnight = self.date_obs - utcoffset
		delta_midnight = np.linspace(-8.0, 8.0, 100) * u.hour
		day_range = midnight + delta_midnight
		frame_date = AltAz(obstime=day_range, location=self.location)

		# Get Moon
		moon = get_moon(day_range)
		moonaltz = moon.transform_to(frame_date)

		# Get Sun
		sun = get_sun(day_range)
		sunaltz = sun.transform_to(frame_date)

		# # Plotting
		# plt.plot(delta_midnight, objairmass)
		# plt.show()
		fig, ax = plt.subplots(tight_layout=True)
		ax.plot(delta_midnight, moonaltz.alt, color="k", ls="--", label="Moon")
		ax.fill_between(delta_midnight, 0*u.deg, 90*u.deg, sunaltz.alt < -0*u.deg, color="0.8", zorder=0)
		ax.fill_between(delta_midnight, 0.0*u.deg, 90.0*u.deg, sunaltz.alt < -18*u.deg, color="0.5", zorder=0)
		ax.set_ylim(0.0, 90.0)
		ax.set_xlim(-8.0, 8.0)
		ax.set_xlabel("Local Time [h]")
		ax.set_ylabel("Altitude [deg]")
		ax.set_title(f"Location: {self.name_location}. UTC offset: {utcoffset}. Date: {self.date_obs.fits.split('T')[0]}")
		ax.grid()

		## Check if multiple objects
		if len(self.names) == 1:
			objaltz = self.coords[0].transform_to(frame_date)
			moon_dist = moon.separation(objaltz)
			objplot = ax.scatter(delta_midnight, objaltz.alt, c=moon_dist, label=self.names[0], s=8, cmap="viridis")
			fig.colorbar(objplot, ax=ax, label="Moon distance [deg]", pad=0.07) #Only if there is one object
		elif len(self.names) > 1:
			for i in range(len(self.names)):
				objaltz = self.coords[i].transform_to(frame_date)
				objplot = ax.plot(delta_midnight, objaltz.alt, linestyle="-.", label=self.names[i], ms=8)
		
		# Legend
		ax.legend()

		# Second axes
		ax2 = ax.twinx()
		ax2.set_ylabel(r"Airmass [sec$(z)$]", fontsize=8.0, rotation=270, labelpad=10.0)
		yyl = np.asarray(ax.get_yticks())
		new_ylabels = np.round(self.__airmass_secz(yyl[::-1]), 2)
		new_ylabels = new_ylabels.tolist()
		new_ylabels[0] = " "
		ax2.set_yticks(yyl)
		ax2.set_yticklabels(new_ylabels, rotation="vertical", fontsize=8.0)

		ax3 = ax.twiny()
		ax3.set_xlabel("UT [h]", fontsize=8.0)
		xx1 = np.asarray(ax.get_xticks())
		ax3.set_xticks(xx1)
		new_upxlabels = xx1 - utcoffset.value + 24
		for i in range(len(new_upxlabels)):
			if new_upxlabels[i] > 24:
				new_upxlabels[i] -= 24
		ax3.set_xticklabels(new_upxlabels, fontsize=8.0)

		# Changing xticks marks to local time
		ll = ax.get_xticks()
		new_xlabels = np.asarray(ll) + 24
		for i in range(len(new_xlabels)):
			if new_xlabels[i] > 24:
				new_xlabels[i] -= 24
		ax.set_xticklabels(new_xlabels)
		
		if self.saveFig:
			fig.savefig(self.date_obs.fits.split('T')[0] + "_" + self.name_location + "_visibility." + self.formatFig)
		elif not self.saveFig:
			plt.show()


if __name__ == "__main__":
	c = SkyCoord.from_name("NGC 3242")
	c2 = SkyCoord.from_name("PN M 1-16")
	ra = [c.ra.value, c2.ra.value]
	dec = [c.dec.value, c2.dec.value]
	names = ["NGC 3242", "PN M 1-16"]
	ov = ObjVisibility(date_obs="2021-4-10", location="spm", ra=ra, dec=dec, names=names)
	ov.staralt()