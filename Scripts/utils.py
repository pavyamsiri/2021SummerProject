"""
Library of utility functions to process the star data.
"""
# Standard libraries
from typing import Tuple

# External libraries
import lightkurve as lk
import matplotlib
import numpy as np
from astropy import units as u
from astropy.timeseries import LombScargle
from lightkurve import MPLSTYLE
from matplotlib import gridspec
from matplotlib import pyplot as plt

# Internal libraries
from star import StarData


def grab_lightcurve(star: str, download_dir: str, fast: bool, mission: str = "TESS") -> lk.lightcurve.LightCurve:
    """Searches for a light curve file using the given star ID. Then downloads all TESS data available,
    stitches the curves together and then cleans the, data returning the result. Returns None in case of errors.

    Parameters
    ----------
    star : str
        a valid star ID i.e. HD, GJ, TIC etc.
    download_dir : str
        a valid directory path to the cache used to store light curve data.
    fast : bool
        if true will instead use only sectors with fast 20 second cadence, else will avoid such sectors
    mission : str
        the specific mission to download from i.e. Kepler, K2, TESS. Default value is "TESS".

    Returns
    -------
    light_curve : lightkurve.lightcurve.LightCurve
       the light curve of the target star
    """

    file_collection = []
    # Search with ID
    search = lk.search_lightcurvefile(star, mission=mission)
    # Abort process if search is empty
    if len(search) == 0:
        return None
    # Create collection of light curve files of the same target
    current_target = None
    for idx, (target, filename) in enumerate(
            search.table[["target_name", "productFilename"]]
        ):
        # Set current target name
        if current_target is None:
            current_target = target
        # If the target is not the current target continue to next observation
        if current_target != target:
            continue
        # Download the lightcurve file
        if fast and "fast" in filename:
            file_collection.append(search[idx].download(download_dir=download_dir))
        elif not fast and "fast" not in filename:
            file_collection.append(search[idx].download(download_dir=download_dir))
        else:
            continue

    # Return light curve if valid light curve files exist
    if len(file_collection) > 0:
        lc_collection = lk.LightCurveFileCollection(file_collection)
        # Stitch the lightcurves together and remove outliers and nans
        light_curve = lc_collection.PDCSAP_FLUX.stitch().remove_outliers().remove_nans()
        return light_curve
    # Otherwise return None
    return None


def grab_periodogram_psd(light_curve: lk.lightcurve.LightCurve, custom: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """Computes the power spectral density spectrum of a given light curve.

    Parameters
    ----------
    light_curve : lightkurve.lightcurve.LightCurve
        the light curve
    custom : bool
        if true will compute the psd using astropy's Lombscargle with custom parameters otherwise it will use lightkurve's
        `to_periodogram` method instead with `"psd"` normalization.

    Returns
    -------
    frequency : np.ndarray
        frequency of the psd
    power : np.ndarray
        power of the psd
    """

    if custom:
        # Calculate power spectral density
        cadence = np.nanmedian(np.diff(light_curve.time) * 24.0 * 60.0 * 60.0)
        # Convert cadence to minutes
        cadence /= 60.0
        # Calculate Nyquist frequency in cycles per day
        nyquist = 1 / (2.0*cadence/60.0/24.0)
        frequency, amplitude = LombScargle(
            light_curve.time, light_curve.flux
        ).autopower(method="fast", samples_per_peak=10, maximum_frequency=nyquist)
        # Convert c/d to muHz
        frequency = 1000.0*frequency/86.4
        # Calculate power density
        bin_width = frequency[1] - frequency[0]
        power = 2.0 * amplitude * np.var(light_curve.flux * 1e6)/(np.sum(amplitude) * bin_width)
    else:
        periodogram = light_curve.to_periodogram(
            method="lombscargle",
            normalization="psd",
            freq_unit=u.microhertz
        )
        frequency = periodogram.frequency.to_value()
        power = periodogram.power.to_value()
    return frequency, power


def create_periodogram(
        light_curve: lk.lightcurve.LightCurve,
        numax: dict = {"ATL": None, "TIC": None},
        min_freq: float = 0.01,
        max_freq: float = 2000.0
) -> Tuple[lk.periodogram.Periodogram, lk.periodogram.Periodogram, lk.periodogram.Periodogram]:
    """Creates three periodograms from a light curve using amplitude normalization. Used specfically for `plotter.py`.

    Parameters
    ----------
    light_curve : lightkurve.lightcurve.LightCurve
        the light curve
    numax : dict
        a dictionary containing an entry for `"ATL"` and `"TIC"`. Add ATL numax and TIC numax here for cutout purposes.
    min_freq : float
        mininum frequency of the medium zoom periodogram assumed to be in muHz. Default value is 0.01.
    max_freq : float
        maximum frequency of the medium zoom periodogram assumed to be in muHz. Default value is 2000.

    Returns
    -------
    normal_pg : lightkurve.periodogram.Periodogram
        an oversampled periodogram from 0.01 muHz to the Nyquist frequency.
    medium_pg : lightkurve.periodogram.Periodogram
        an oversampled periodogram from `0.3*numax` to `1.7*numax` or from `min_freq` to `max_freq` if `numax` is empty.
    small_pg : lightkurve.periodogram.Periodogram
        an oversampled periodogram from 0.01 muHz to 100 muHz.
    """

    # Create periodogram from lightcurve
    pg_method = "lombscargle"
    norm_method = "amplitude"

    # Periodogram up to Nyquist frequency
    normal_pg = light_curve.to_periodogram(
        method=pg_method,
        normalization=norm_method,
        minimum_frequency=0.01,
        oversample_factor=10,
        freq_unit=u.microhertz
    )

    atl_numax = numax["ATL"]
    tic_numax = numax["TIC"]
    # ATL numax estimate exists but is the only one
    if atl_numax is not None and (tic_numax is None or len(tic_numax) == 0):
        min_freq = max(0.3*atl_numax, 0.01)
        max_freq = min(1.7*atl_numax, normal_pg.frequency[-1].value)
    # ATL numax does not exist but TIC estimates do
    elif atl_numax is None and tic_numax is not None and len(tic_numax) >= 3:
        min_freq = max(0.3*tic_numax[1], 0.01)
        max_freq = min(1.7*tic_numax[1], normal_pg.frequency[-1].value)
    elif atl_numax is not None and tic_numax is not None:
        min_freq = max(max(0.3*tic_numax[1], 0.01), max(0.3*atl_numax, 0.01))
        max_freq = min(min(1.7*tic_numax[1], normal_pg.frequency[-1].value), min(1.7*atl_numax, normal_pg.frequency[-1].value))
    # No estimates for numax are available so just use given values
    
    # Periodogram from 0.1 microHertz to 2000 microHertz
    # or 0.3*numax to 1.7*numax microHertz around ATL's estimated numax
    medium_pg = light_curve.to_periodogram(
        method=pg_method,
        normalization=norm_method,
        minimum_frequency=min_freq,
        maximum_frequency=max_freq,
        oversample_factor=10,
        freq_unit=u.microhertz
    )

    # Periodogram from 0.1 microHertz to 100 microHertz
    small_pg = light_curve.to_periodogram(
        method=pg_method,
        normalization=norm_method,
        minimum_frequency=0.01,
        maximum_frequency=100,
        oversample_factor=10,
        freq_unit=u.microhertz
    )

    return (normal_pg, medium_pg, small_pg)


def create_plots(star_data: StarData,
        light_curve: lk.lightcurve.LightCurve,
        numax: dict = {"ATL": None, "TIC": None},
        min_freq: float = 0.01,
        max_freq: float = 2000.0
) -> matplotlib.figure.Figure:
    """Plots the lightcurve and periodograms as three subplots, and returns the figure.

    Parameters
    ----------
    star_data : star.StarData
        data pertaining to the target star
    light_curve : lightkurve.lightcurve.LightCurve
        the light curve of the target
    numax : dict
        a dictionary containing an entry for `"ATL"` and `"TIC"`. Add ATL numax and TIC numax here for cutout purposes.
    min_freq : float
        mininum frequency of the medium zoom periodogram assumed to be in muHz. Default value is 0.01.
    max_freq : float
        maximum frequency of the medium zoom periodogram assumed to be in muHz. Default value is 2000.

    Returns
    -------
    plot : matplotlib.figure.Figure
        the final plot figure
    """

    # Create plots
    (normal_pg, medium_pg, small_pg) = create_periodogram(
        light_curve,
        numax=numax,
        min_freq=min_freq,
        max_freq=max_freq
    )

    with plt.style.context(MPLSTYLE):
        fig = plt.figure(constrained_layout=True, figsize=(20, 30))
        spec = gridspec.GridSpec(ncols=2, nrows=3, figure=fig)
        fig_ax1 = fig.add_subplot(spec[0, 0])
        fig_ax2 = fig.add_subplot(spec[0, 1])
        fig_ax3 = fig.add_subplot(spec[1, :])
        fig_ax4 = fig.add_subplot(spec[2, :])

        # Light curve
        fig_ax1.scatter(light_curve.time, light_curve.flux)
        fig_ax1.set_title(f"{star_data.star} light curve")
        fig_ax1.set_xlabel("Time - 2457000 [BTJD days]", fontsize=21)
        fig_ax1.set_ylabel("Normalized Flux", fontsize=21)
        fig_ax1.tick_params(axis="both", which="major", labelsize=21)

        # Periodogram labels
        xlabel = "Frequency ({})".format(normal_pg.frequency.unit.to_string('latex'))

        # log-log periodogram
        smooth_pg = normal_pg.smooth(method="boxkernel", filter_width=1.0)
        fig_ax2.plot(normal_pg.frequency, normal_pg.power)
        fig_ax2.plot(smooth_pg.frequency, smooth_pg.power, "r")
        fig_ax2.set_title(f"{star_data.star} periodogram (log-log scale)")
        fig_ax2.set_xlabel(xlabel, fontsize=21)
        fig_ax2.set_xscale("log")
        fig_ax2.set_xlim(0.1)
        fig_ax2.set_ylabel("Amplitude", fontsize=21)
        fig_ax2.set_yscale("log")
        fig_ax2.tick_params(axis="both", which="major", labelsize=21)

        # Smallest periodogram
        fig_ax3.plot(small_pg.frequency, small_pg.power, linewidth=3)
        fig_ax3.set_title(f"{star_data.star} periodogram")
        fig_ax3.set_xlabel(xlabel, fontsize=21)
        fig_ax3.set_xlim((small_pg.frequency[0].value, small_pg.frequency[-1].value))
        fig_ax3.set_ylabel("Amplitude", fontsize=21)
        fig_ax3.set_ylim(bottom=0)
        fig_ax3.tick_params(axis="both", which="major", labelsize=21)

        # Medium size periodogram
        fig_ax4.plot(medium_pg.frequency, medium_pg.power, zorder=1)
        fig_ax4.set_title(f"{star_data.star} periodogram")
        fig_ax4.set_xlabel(xlabel, fontsize=21)
        fig_ax4.set_xlim((medium_pg.frequency[0].value, medium_pg.frequency[-1].value))
        fig_ax4.set_ylabel("Amplitude", fontsize=21)
        fig_ax4.set_ylim(bottom=0)
        fig_ax4.tick_params(axis="both", which="major", labelsize=21)
        
        estimate_colours = ["xkcd:green", "xkcd:red", "xkcd:orange", "xkcd:blue"]
        
        atl_numax = numax["ATL"]
        tic_numax = numax["TIC"]
        if atl_numax is not None:
            fig_ax4.axvline(atl_numax, c=estimate_colours[0], lw=5, ls="dashed")
        if tic_numax is not None:
            for idx in range(len(tic_numax)):
                fig_ax4.axvline(tic_numax[idx], c=estimate_colours[idx+1], lw=5, ls="dashed")
                    

        fig_ax4.text(
            0.8,
            0.65,
            star_data,
            fontsize=20,
            color="black",
            transform=fig_ax4.transAxes
        )

        return fig


def calculate_angular_diameter(
        numax: float,
        numax_error: float,
        dnu: float,
        dnu_error: float,
        teff: float,
        teff_error: float,
        distance: float,
        distance_error: float
    ) -> Tuple[float, float]:
    """Calculates angular diameter and angular diameter error using asteroseismic parameters.

    Parameters
    ----------
    numax : float
        the numax frequency in muHz
    numax_error : float
        the error in numax frequency in muHz
    dnu : float
        the large frequency separation in muHz
    dnu_error : float
        the error in the large frequency separation in muHz
    teff : float
        the effective temperature in K
    teff_error : float
        the error in effective temperature in K
    distance : float
        the distance to the star in parsecs
    distance_error : float
        the error in distance to the star in parsecs

    Returns
    -------
    ang_diam : float
        the angular diameter of the star in mas
    ang_diam_error : float
        the error in the angular diameter of the star in mas
    """

    # Solar values
    NUMAX_SOL = (3090, 30) # microhertz | Huber et al. 2011
    DELTANU_SOL = (135.1, 0.1) # microhertz | Huber et al. 2011
    TEFF_SOL = (5772.0, 0.8) # Kelvin    | Prsa et al. 2016

    # Asteroseismic stellar radius
    radius = numax/NUMAX_SOL[0] * np.power(dnu/DELTANU_SOL[0], -2) * np.power(teff/TEFF_SOL[0], 0.5)
    # Calculate angular diameter
    ang_diam = 206265*(2*696340*radius/(30856775812800*distance))

    # Calculate error
    radius_error = radius*np.sqrt(
        np.power(numax_error/numax + NUMAX_SOL[1]/NUMAX_SOL[0], 2)
        + np.power(2*(dnu_error/dnu + DELTANU_SOL[1]/DELTANU_SOL[0]), 2)
        + np.power(0.5*(teff_error/teff + TEFF_SOL[1]/TEFF_SOL[0]), 2)
    )
    ang_diam_error = ang_diam*(radius_error/radius + distance_error/distance)

    # Convert from as to mas
    return 1000*ang_diam, 1000*ang_diam_error
