"""Create High-Level Science Products for a KLIP-processed coronagraphic exposure."""

import logging
import math

import numpy as np
from stdatamodels.jwst import datamodels

log = logging.getLogger(__name__)

__all__ = ["snr_image", "contrast_curve"]


def snr_image(target_model):
    """
    Compute the signal-to-noise ratio map for a KLIP-processed image.

    The SNR image is computed as the ratio of the input image DATA and ERR arrays.

    Parameters
    ----------
    target_model : ImageModel
        Science target ImageModel

    Returns
    -------
    output_model : ImageModel
        ImageModel with signal-to_noise data attribute
    """
    # Initialize the output model as a copy of the input model
    output_model = target_model.copy()

    # Compute the SNR image, which is simply the ratio of the target
    # data and err images
    output_model.data = target_model.data / target_model.err

    return output_model


def contrast_curve(target_model, width):
    """
    Compute "contrast curve" data for a PSF-subtracted image.

    What's computed is the 1-sigma noise in a series of concentric annuli
    centered at the geometric center of the image. The width of the annuli
    is set by the input argument "width". A mask array is used iteratively to
    exclude all pixels outside of each annulus when the standard deviation
    is computed. Annuli are defined to the nearest pixel. No partial-pixel
    computations are performed.

    Parameters
    ----------
    target_model : ImageModel
        The science target ImageModel
    width : float
        Width of the annuli, in pixels

    Returns
    -------
    output_model : ContrastModel
        ContrastModel with a table of radius and standard deviation of non-masked pixels
    """
    # Get the target array size and center
    nrows = target_model.data.shape[0]
    ncols = target_model.data.shape[1]
    crow = int(round(nrows / 2)) - 1
    ccol = int(round(ncols / 2)) - 1

    # Create a series of annuli inner radii that starts at zero, extends
    # to the nearest edge of the image, and goes in increments of "width"
    limits = list(range(0, min(crow, ccol), width))

    # Initialize the mask array and the output lists of radii and sigma values
    mask = target_model.data * 0.0
    radius = []
    sigma = []

    # Loop over the list of annuli radii
    for rmin in limits:
        rmax = rmin + width
        rmid = (rmin + rmax) / 2.0

        # Initialize all pixels in the mask to NaN
        mask[:] = np.nan

        # Loop over all pixels in the data array, determining which ones are
        # within the limits of the current annulus
        for r in range(nrows):
            for c in range(ncols):
                # Compute the distance of this pixel from the image center
                d = math.sqrt((r - crow) ** 2 + (c - ccol) ** 2)

                # If this pixel is within the current annulus,
                # reset its mask value to 1
                if d >= rmin and d < rmax:
                    mask[r, c] = 1

        # Compute the standard deviation of all unmasked pixels
        std = np.nanstd(mask * target_model.data)
        log.debug(" rmin=%d, rmax=%d, rmid=%g, sigma=%g", rmin, rmax, rmid, std)

        # Append the radius and standard deviation values for this annulus
        # to the output lists
        radius.append(rmid)
        sigma.append(std)

    # Convert the lists to arrays
    radii = np.asarray(radius)
    sigmas = np.asarray(sigma)

    # Create the output contrast curve data model
    output_model = datamodels.ContrastModel(contrast_table=list(zip(radii, sigmas, strict=True)))
    output_model.update(target_model)

    return output_model
