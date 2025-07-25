"""
Utilities for velocity aberration correction.

Script to add velocity aberration correction information to the FITS
files provided to it on the command line (one or more).

It assumes the following keywords are present in the file header:

* JWST_DX (km/sec)
* JWST_DY (km/sec)
* JWST_DZ (km/sec)
* RA_REF (deg)
* DEC_REF (deg)

The keywords added are:

* VA_SCALE (dimensionless scale factor)

It does not currently place the new keywords in any particular location
in the header other than what is required by the standard.
"""

import logging

import numpy as np
from gwcs.geometry import CartesianToSpherical, SphericalToCartesian
from scipy.constants import speed_of_light

import jwst.datamodels as dm
from jwst.datamodels import Level1bModel  # type: ignore[attr-defined]

# Configure logging
logger = logging.getLogger(__name__)

SPEED_OF_LIGHT = speed_of_light / 1000  # km / s

__all__ = ["compute_va_effects_vector", "compute_va_effects", "add_dva"]


def compute_va_effects_vector(velocity_x, velocity_y, velocity_z, u):
    """
    Compute velocity aberration effects scale factor.

    Computes constant scale factor due to velocity aberration as well as
    corrected ``RA`` and ``DEC`` values, in vector form.

    Parameters
    ----------
    velocity_x, velocity_y, velocity_z : float
        The components of the velocity of JWST, in km / s with respect to
        the Sun.  These are celestial coordinates, with x toward the
        vernal equinox, y toward right ascension 90 degrees and declination
        0, z toward the north celestial pole.

    u : numpy.array([u0, u1, u2])
        The vector form of right ascension and declination of the target (or some other
        point, such as the center of a detector) in the barycentric coordinate
        system.  The equator and equinox should be the same as the coordinate
        system for the velocity.

    Returns
    -------
    scale_factor : float
        Multiply the nominal image scale (e.g., in degrees per pixel) by
        this value to obtain the image scale corrected for the "aberration
        of starlight" due to the velocity of JWST with respect to the Sun.

    u_corr : numpy.array([ua0, ua1, ua2])
        Apparent position vector in the moving telescope frame.
    """
    beta = np.array([velocity_x, velocity_y, velocity_z]) / SPEED_OF_LIGHT
    beta2 = np.dot(beta, beta)  # |beta|^2
    if beta2 == 0.0:
        logger.warning("Observatory speed is zero. Setting VA scale to 1.0")
        return 1.0, u

    u_beta = np.dot(u, beta)
    igamma = np.sqrt(1.0 - beta2)  # inverse of usual gamma
    scale_factor = (1.0 + u_beta) / igamma

    # Algorithm below is from Colin Cox notebook.
    # Also see: Instrument Science Report OSG-CAL-97-06 by Colin Cox (1997).
    u_corr = (igamma * u + beta * (1.0 + (1.0 - igamma) * u_beta / beta2)) / (1.0 + u_beta)

    return scale_factor, u_corr


def compute_va_effects(velocity_x, velocity_y, velocity_z, ra, dec):
    """
    Compute velocity aberration effects.

    Computes constant scale factor due to velocity aberration as well as
    corrected ``RA`` and ``DEC`` values.

    Parameters
    ----------
    velocity_x, velocity_y, velocity_z : float
        The components of the velocity of JWST, in km / s with respect to
        the Sun.  These are celestial coordinates, with x toward the
        vernal equinox, y toward right ascension 90 degrees and declination
        0, z toward the north celestial pole.

    ra, dec : float
        The right ascension and declination of the target (or some other
        point, such as the center of a detector) in the barycentric coordinate
        system.  The equator and equinox should be the same as the coordinate
        system for the velocity.

    Returns
    -------
    scale_factor : float
        Multiply the nominal image scale (e.g., in degrees per pixel) by
        this value to obtain the image scale corrected for the "aberration
        of starlight" due to the velocity of JWST with respect to the Sun.

    apparent_ra : float
        Apparent star position in the moving telescope frame.

    apparent_dec : float
        Apparent star position in the moving telescope frame.
    """
    u = np.asanyarray(SphericalToCartesian()(ra, dec))
    scale_factor, u_corr = compute_va_effects_vector(velocity_x, velocity_y, velocity_z, u)
    apparent_ra, apparent_dec = CartesianToSpherical()(*u_corr)
    return scale_factor, apparent_ra, apparent_dec


def add_dva(filename, force_level1bmodel=True):
    """
    Determine velocity aberration.

    Given the name of a valid partially populated level 1b JWST file,
    determine the velocity aberration scale factor and apparent target position
    in the moving (telescope) frame.

    It presumes all the accessed keywords are present (see first block).

    Parameters
    ----------
    filename : str
        The name of the file to be updated.
    force_level1bmodel : bool, optional
        If True, the input file will be force-opened as a Level1bModel.  If False,
        the file will be opened using the generic DataModel.  The default is True.
    """
    if force_level1bmodel:
        model = Level1bModel(filename)
    else:
        model = dm.open(filename)
    scale_factor, apparent_ra, apparent_dec = compute_va_effects(
        velocity_x=model.meta.ephemeris.velocity_x_bary,
        velocity_y=model.meta.ephemeris.velocity_y_bary,
        velocity_z=model.meta.ephemeris.velocity_z_bary,
        ra=model.meta.wcsinfo.ra_ref,
        dec=model.meta.wcsinfo.dec_ref,
    )

    # update header
    model.meta.velocity_aberration.scale_factor = scale_factor
    model.meta.velocity_aberration.va_ra_ref = apparent_ra
    model.meta.velocity_aberration.va_dec_ref = apparent_dec
    model.save(filename)
