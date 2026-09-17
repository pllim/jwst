import logging

from gwcs.wcs import WCS

log = logging.getLogger(__name__)

__all__ = [
    "check_exptime",
    "combine_1d_spectra",
]


class _InputSpectrumModel:
    """
    Class to hold common attributes from different models.

    Attributes
    ----------
    wavelength : ndarray
        Input wavelength.
    flux : ndarray
        Input flux.
    flux_error : ndarray
        Input error on the flux.
    surf_bright : ndarray
        Input surface brightness.
    sb_error : ndarray
        Input error on the surface brightness.
    dq : ndarray
        Input DQ array.
    weight_factor : float
        Weight factor to apply to all the spectral elements.
    right_ascension : float
        RA value for the spectrum.
    declination : float
        Dec value for the spectrum.
    name : str
        Slit name for the spectrum.
    """

    def __init__(
        self,
        wavelength,
        flux,
        flux_error,
        surf_bright,
        sb_error,
        dq,
        weight_factor,
        right_ascension,
        declination,
        name,
    ):
        self.wavelength = wavelength
        self.flux = flux
        self.flux_error = flux_error
        self.surf_bright = surf_bright
        self.sb_error = sb_error
        self.dq = dq
        self.weight_factor = weight_factor  # FIXME: Don't want scalar, just remove???
        self.right_ascension = right_ascension
        self.declination = declination
        self.name = name


def _inputlist_from_multispec(ms, spec, exptime_key):
    """Create list of _InputSpectrumModel from MultiSpecModel or MRSMultiSpecModel."""
    wavelength = spec.spec_table.field("wavelength")
    flux = spec.spec_table.field("flux")
    flux_error = spec.spec_table.field("flux_error")
    surf_bright = spec.spec_table.field("surf_bright")
    sb_error = spec.spec_table.field("sb_error")
    dq = spec.spec_table.field("dq")
    weight_factor = _exptime_to_fac(exptime_key, ms)
    right_ascension, declination = _ms_spec_to_radec(ms, spec)
    name = spec.name
    return [
        _InputSpectrumModel(
            wavelength,
            flux,
            flux_error,
            surf_bright,
            sb_error,
            dq,
            weight_factor,
            right_ascension,
            declination,
            name,
        )
    ]


def _inputlist_from_tso(ms, spec, exptime_key):
    """Create list of _InputSpectrumModel from TSOMultiSpecModel."""
    raise NotImplementedError


def _inputlist_from_wfss(ms, spec, exptime_key):
    """Create list of _InputSpectrumModel from WFSSMultiSpecModel."""
    raise NotImplementedError


def _ms_spec_to_radec(ms, spec):
    """Find RA and Dec for given spectrum/model."""
    if hasattr(spec.meta, "wcs") and isinstance(spec.meta.wcs, WCS):
        ra, dec = spec.meta.wcs(0.0)[:2]
    else:
        ra = ms.meta.target.ra
        dec = ms.meta.target.dec
        # This exception is hit for NIRISS and NIRCam WFSS data,
        # for which it doesn't matter anyway, since the RA and Dec are not used
        # in any meaningful way to combine the spectra.
        log.debug("There is no WCS in the input. Getting RA, Dec from target metadata.")
    return ra, dec


def _exptime_to_fac(exptime_key, ms):
    """Compute weight factor for given model and exptime_key."""
    if exptime_key == "integration_time":
        fac = ms.meta.exposure.integration_time
    elif exptime_key == "exposure_time":
        fac = ms.meta.exposure.exposure_time
    else:  # "unit_weight" (already validated in check_exptime)
        fac = 1.0
    return fac


def check_exptime(exptime_key):
    """
    Check exptime_key for validity.

    This function checks ``exptime_key``.  If it is valid, the corresponding
    value used by the metadata interface will be returned.  This will be
    either "integration_time" or "exposure_time".  If it is invalid,
    "unit weight" will be returned (meaning that a weight of 1 will be
    used when averaging spectra), and a warning will be logged.

    Parameters
    ----------
    exptime_key : str
        A keyword or string indicating what value (integration time or
        exposure time) should be used as a weight when combing spectra.

    Returns
    -------
    exptime_key : str
        The value will be either "integration_time", "exposure_time",
        or "unit_weight".
    """
    exptime_lwr = exptime_key.lower()
    if exptime_lwr.startswith("integration") or exptime_lwr == "effinttm":
        exptime_key = "integration_time"
        log.info("Using integration time as the weight.")
    elif exptime_lwr.startswith("exposure") or exptime_lwr == "effexptm":
        exptime_key = "exposure_time"
        log.info("Using exposure time as the weight.")
    elif exptime_lwr == "unit_weight" or exptime_lwr == "unit weight":
        exptime_key = "unit_weight"
        log.info("Using weight = 1.")
    else:
        log.warning(f"Don't understand exptime_key = '{exptime_key}'; using unit weight.")
        log.info("The options for exptime_key are:")
        log.info("  integration_time, effinttm, exposure_time, effexptm, unit_weight, unit weight")
        exptime_key = "unit_weight"

    return exptime_key


def _read_input_spectra(input_model, exptime_key, input_spectra):
    """
    Read input spectra from a datamodel.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.MultiSpecModel`, \
                  `~stdatamodels.jwst.datamodels.MRSSpecModel`, \
                  `~stdatamodels.jwst.datamodels.TSOMultiSpecModel`, or \
                  `~stdatamodels.jwst.datamodels.WFSSMultiSpecModel`
        A datamodel with a ``spec`` attribute, containing spectra.
        For TSO and WFSS, integrations in the spectral table rows
        are expanded into separate spectra. For WFSS, each source ID
        is processed separately.
    exptime_key : str
        Exposure time key to use for weighting.
    input_spectra : dict
        Dictionary to hold input spectra, keyed by spectral order;
        Updated in place.


    """


def combine_1d_spectra(input_model, exptime_key, sigma_clip=None):
    """
    Combine the input spectra.

    Parameters
    ----------
    input_model : `~jwst.datamodels.container.ModelContainer`, \
                  `~stdatamodels.jwst.datamodels.MultiSpecModel`, \
                  `~stdatamodels.jwst.datamodels.MRSMultiSpecModel`, \
                  `~stdatamodels.jwst.datamodels.TSOMultiSpecModel`, or \
                  `~stdatamodels.jwst.datamodels.WFSSMultiSpecModel`
        The input spectra.
        They may have different spectral orders
        or wavelengths but should all share the same target.
        Each datamodel must have a ``spec`` attribute, containing spectra.
        For TSO and WFSS, integrations in the spectral table rows
        are expanded into separate spectra. For WFSS, each source ID
        is processed separately.
        May be updated in-place if processing is skipped.
    exptime_key : str
        A string identifying which keyword to use to get the exposure time,
        which is used as a weight when combining spectra.  The value should
        be one of:  "exposure_time", "integration_time",
        or "unit_weight".
    sigma_clip : float or None, optional
        Factor for clipping outliers in spectral combination.
        Default is None (no clipping).

    Returns
    -------
    output_model : `~stdatamodels.jwst.datamodels.MultiCombinedSpecModel` or \
                   `~stdatamodels.jwst.datamodels.WFSSMultiCombinedSpecModel`
        A combined spectra datamodel.
    """
    log.debug(f"Using exptime_key = {exptime_key}.")

    exptime_key = check_exptime(exptime_key)
