import numpy as np
import matplotlib.pyplot as plt
from labellines import labelLines
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_moon, get_sun
import pytz
from datetime import datetime
import logging as _log
from typing import Optional, Union, List


class ObjVisibility(object):
    """
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
    """

    def __init__(self, date_obs: str, location: str) -> None:
        """
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
        """
        self._i = 0
        self.list_names = []
        self.plot_handler = []

        self.name_location = location
        self.location = EarthLocation.of_site(location)
        self.date_obs = Time(date_obs + " 00:00:00")
        self.time_zone = self.location.info.meta["timezone"]
        self.fig, self.ax = plt.subplots(tight_layout=True)

    def __airmass_secz(self, degrees: np.ndarray) -> np.ndarray:
        """
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

        """
        return 1.0 / np.cos(np.deg2rad(degrees))

    def _add_moon(self, delta_midnight, day_range, frame_date):
        """
        Retrieves the moon Location.
        """
        ax = self.ax

        moon = get_moon(day_range)
        moonaltz = moon.transform_to(frame_date)

        ax.plot(delta_midnight, moonaltz.alt, color="k", ls="--", label="Moon")

        return moon

    def add_object(
        self,
        ra: float,
        dec: float,
        name: str,
        color: Optional[str] = None,
        moon_dist: bool = False,
    ) -> None:
        """
        Adds an object to the plot.
        """
        coords = SkyCoord(ra, dec, unit="deg", frame="icrs")
        self._i += 1
        self.list_names.append(f"{self._i} {name}")
        ax = self.ax

        tz = pytz.timezone(self.time_zone)
        utcoffset = tz.utcoffset(datetime.now()).total_seconds() / 60.0 / 60.0 * u.hour

        # Plotting interval
        midnight = self.date_obs - utcoffset
        delta_midnight = np.linspace(-8.0, 8.0, 100) * u.hour
        day_range = midnight + delta_midnight
        frame_date = AltAz(obstime=day_range, location=self.location)

        objaltz = coords.transform_to(frame_date)

        if moon_dist:
            moon = get_moon(day_range)
            moon_dist = moon.separation(objaltz)
            objplot = ax.scatter(
                delta_midnight,
                objaltz.alt,
                c=moon_dist,
                label=name,
                s=8,
                cmap="viridis",
            )
            self.fig.colorbar(
                objplot, ax=ax, label="Moon distance [deg]", pad=0.07
            )  # Only if there is one object
        else:
            (objplot,) = ax.plot(
                delta_midnight,
                objaltz.alt,
                linestyle="-.",
                label=f"{self._i}",  # {self.names[i]}",
                ms=8,
                color=color,
                lw=0.5,
            )
            self.plot_handler.append(objplot)

    def staralt(
        self,
        ra: Union[List[float], float],
        dec: Union[List[float], float],
        names: Union[List[str], str],
        color: Optional[str] = None,
        linelabels: Optional[str] = False,
        add_moon: Optional[str] = False,
    ) -> None:
        """
        Function to plot the staralt visibility map with moon position.

        Parameters
        ----------
                color:
                        List[str] or str. List same size as object len.

        Returns
        -------
                None.
        """
        if not isinstance(ra, list):
            ra = [ra]
            dec = [dec]
        if not isinstance(names, list):
            names = [names]

        tz = pytz.timezone(self.time_zone)
        utcoffset = tz.utcoffset(datetime.now()).total_seconds() / 60.0 / 60.0 * u.hour

        # Plotting interval
        midnight = self.date_obs - utcoffset
        delta_midnight = np.linspace(-8.0, 8.0, 100) * u.hour
        day_range = midnight + delta_midnight
        frame_date = AltAz(obstime=day_range, location=self.location)

        if add_moon:
            _ = self._add_moon(delta_midnight, day_range, frame_date)

        # Get Sun
        sun = get_sun(day_range)
        sunaltz = sun.transform_to(frame_date)

        # # Plotting
        # plt.plot(delta_midnight, objairmass)
        # plt.show()
        ax = self.ax
        ax.fill_between(
            delta_midnight.value,
            0,
            90,
            sunaltz.alt.value < -0,
            color="0.8",
            zorder=0,
        )
        ax.fill_between(
            delta_midnight.value,
            0.0,
            90.0,
            sunaltz.alt.value < -18,
            color="0.5",
            zorder=0,
        )
        ax.set_ylim(0.0, 90.0)
        ax.set_xlim(-8.0, 8.0)
        ax.set_xlabel("Local Time [h]")
        ax.set_ylabel("Altitude [deg]")
        ax.set_title(
            f"Location: {self.name_location}. UTC offset: {utcoffset}. Date: {self.date_obs.fits.split('T')[0]}"
        )
        ax.grid()

        ## Check if multiple objects
        if len(names) == 1:
            self.add_object(
                ra=ra[0], dec=dec[0], name=names[0], color=color, moon_dist=add_moon
            )
        else:
            for i in range(len(names)):
                c = color[i] if isinstance(color, list) else color
                self.add_object(
                    ra=ra[i], dec=dec[i], name=names[i], color=c, moon_dist=False
                )

            if linelabels:
                print(self.list_names, self._i, len(ax.get_lines()))
                labelLines(ax.get_lines(), fontsize=5, drop_label=True)

    def savePlot(self, show: bool = False, formatFig: str = "pdf") -> None:
        """
        Saves the plot to a file.
        """
        ax = self.ax
        tz = pytz.timezone(self.time_zone)
        utcoffset = tz.utcoffset(datetime.now()).total_seconds() / 60.0 / 60.0 * u.hour

        # Legend
        ax.legend(handles=self.plot_handler, labels=self.list_names, fontsize=5)

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

        if not show:
            self.fig.savefig(
                self.date_obs.fits.split("T")[0]
                + "_"
                + self.name_location
                + "_visibility."
                + formatFig
            )
        else:
            plt.show()


if __name__ == "__main__":
    c = SkyCoord.from_name("PN M 1-16")
    c2 = SkyCoord.from_name("NGC 3242")
    ra = [c.ra.value, c2.ra.value]
    dec = [c.dec.value, c2.dec.value]
    names = ["M 1-16", "N3242"]
    ov = ObjVisibility(
        date_obs="2021-4-10",
        location="spm",
    )
    ov.staralt(ra=ra, dec=dec, names=names, color="C0", linelabels=True, add_moon=True)
    ov.savePlot(show=False, formatFig="pdf")
