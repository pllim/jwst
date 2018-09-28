"""
Test NIRCAM grism WCS transformations.

Notes:
These test the stability of the WFSS and TSO transformations based on the
the reference file that is returned from CRDS. The absolute validity
of the results is verified by the team based on the specwcs reference
file.

"""
from numpy.testing.utils import assert_allclose
import pytest


from astropy.io import fits
from gwcs import wcs

from ...datamodels.image import ImageModel
from ...datamodels import CubeModel
from ..assign_wcs_step import AssignWcsStep
from .. import nircam


# Allowed settings for nircam
tsgrism_filters = ['F277W', 'F444W', 'F322W2', 'F356W']

nircam_wfss_frames = ['grism_detector', 'detector', 'v2v3', 'world']

nircam_tsgrism_frames = ['grism_detector', 'full_detector', 'v2v3', 'world']

nircam_imaging_frames = ['detector', 'v2v3', 'world']


# Default wcs information
wcs_wfss_kw = {'wcsaxes': 2, 'ra_ref': 53.1423683802, 'dec_ref': -27.8171119969,
               'v2_ref': 86.103458, 'v3_ref': -493.227512, 'roll_ref': 45.04234459270135,
               'crpix1': 1024.5, 'crpix2': 1024.5,
               'crval1': 53.1423683802, 'crval2': -27.8171119969,
               'cdelt1': 1.74460027777777e-05, 'cdelt2': 1.75306861111111e-05,
               'ctype1': 'RA---TAN', 'ctype2': 'DEC--TAN',
               'pc1_1': -1, 'pc1_2': 0,
               'pc2_1': 0, 'pc2_2': 1,
               'cunit1': 'deg', 'cunit2': 'deg',
               }

wcs_tso_kw = {'wcsaxes': 2, 'ra_ref': 86.9875, 'dec_ref': 23.423,
              'v2_ref': 95.043034, 'v3_ref': -556.150466, 'roll_ref': 359.9521,
              'crpix1': 887.0, 'crpix2': 35.0,
              'cdelt1': 1.76686111111111e-05, 'cdelt2': 1.78527777777777e-05,
              'ctype1': 'RA---TAN', 'ctype2': 'DEC--TAN',
              'pc1_1': -1, 'pc1_2': 0,
              'pc2_1': 0, 'pc2_2': 1,
              'cunit1': 'deg', 'cunit2': 'deg',
              }


def create_hdul(detector='NRCALONG', channel='LONG', module='A',
                filtername='F444W', exptype='NRC_IMAGE', pupil='GRISMR',
                subarray='FULL', wcskeys=wcs_wfss_kw):
    hdul = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header['telescop'] = "JWST"
    phdu.header['filename'] = "test+" + filtername
    phdu.header['instrume'] = 'NIRCAM'
    phdu.header['channel'] = channel
    phdu.header['detector'] = detector
    phdu.header['FILTER'] = filtername
    phdu.header['PUPIL'] = pupil
    phdu.header['MODULE'] = module
    phdu.header['time-obs'] = '8:59:37'
    phdu.header['date-obs'] = '2017-09-05'
    phdu.header['exp_type'] = exptype
    scihdu = fits.ImageHDU()
    scihdu.header['EXTNAME'] = "SCI"
    scihdu.header['SUBARRAY'] = subarray
    scihdu.header.update(wcskeys)
    hdul.append(phdu)
    hdul.append(scihdu)
    return hdul


def create_wfss_wcs(pupil, filtername='F444W'):
    """Help create WFSS GWCS object."""
    hdul = create_hdul(exptype='NRC_WFSS', filtername=filtername, pupil=pupil)
    im = ImageModel(hdul)
    ref = get_reference_files(im)
    pipeline = nircam.create_pipeline(im, ref)
    wcsobj = wcs.WCS(pipeline)
    return wcsobj


def create_imaging_wcs():
    hdul = create_hdul()
    image = ImageModel(hdul)
    ref = get_reference_files(image)
    pipeline = nircam.create_pipeline(image, ref)
    wcsobj = wcs.WCS(pipeline)
    return wcsobj


def create_tso_wcs(filtername=tsgrism_filters[0], subarray="SUBGRISM256"):
    """Help create tsgrism GWCS object."""
    hdul = create_hdul(exptype='NRC_TSGRISM', pupil='GRISMR',
                       filtername=filtername, detector='NRCALONG',
                       subarray=subarray, wcskeys=wcs_tso_kw)
    im = CubeModel(hdul)
    ref = get_reference_files(im)
    pipeline = nircam.create_pipeline(im, ref)
    wcsobj = wcs.WCS(pipeline)
    return wcsobj


def get_reference_files(datamodel):
    """Get the reference files associated with the step."""
    refs = {}
    step = AssignWcsStep()
    for reftype in AssignWcsStep.reference_file_types:
        refs[reftype] = step.get_reference_file(datamodel, reftype)
    return refs


def test_nircam_wfss_available_frames():
    """Make sure that the expected GWCS reference frames are created."""
    for p in ['GRISMR', 'GRISMC']:
        wcsobj = create_wfss_wcs(p)
        available_frames = wcsobj.available_frames
        assert all([a == b for a, b in zip(nircam_wfss_frames, available_frames)]), "available frame mismatch"


def test_nircam_tso_available_frames():
    """Make sure that the expected GWCS reference frames for TSO are created."""
    wcsobj = create_tso_wcs()
    available_frames = wcsobj.available_frames
    assert all([a == b for a, b in zip(nircam_tsgrism_frames, available_frames)]), "available frame mismatch"


def traverse_wfss_trace(pupil):
    """Make sure that the WFSS dispersion polynomials are reversable.""" 
    wcsobj = create_wfss_wcs(pupil)
    detector_to_grism = wcsobj.get_transform('detector', 'grism_detector')
    grism_to_detector = wcsobj.get_transform('grism_detector', 'detector')

    # check the round trip, grism pixel 100,100, source at 110,110,order 1
    xgrism, ygrism, xsource, ysource, orderin = (100, 100, 110, 110, 1)
    x0, y0, lam, order = grism_to_detector(xgrism, ygrism, xsource, ysource, orderin)
    x, y, xdet, ydet, orderdet = detector_to_grism(x0, y0, lam, order)

    assert x0 == xsource, "Output grism x-source location changed"
    assert y0 == ysource, "Output grism y-source location changed"
    assert order == orderin, "Output order differs from input grism order"
    assert xdet == xsource, "Grism x-detector pixel changed"
    assert ydet == ysource, "Grism y-detector pixel changed"
    assert orderdet == orderin, "Grism order out doesn't match order in"


def test_traverse_wfss_grisms():
    """Check the dispersion polynomials for each grism."""
    for pupil in ['GRISMR', 'GRISMC']:
        traverse_wfss_trace(pupil)


def test_traverse_tso_grism():
    """Make sure that the TSO dispersion polynomials are reversable."""
    wcsobj = create_tso_wcs()
    detector_to_grism = wcsobj.get_transform('full_detector', 'grism_detector')
    grism_to_detector = wcsobj.get_transform('grism_detector', 'full_detector')

    # TSGRISM always has same source locations
    # takes x,y,order -> ra, dec, wave, order
    xin, yin, order = (100, 100, 1)

    x0, y0, lam, orderdet = grism_to_detector(xin, yin, order)
    x, y, orderdet = detector_to_grism(x0, y0, lam, order)

    assert x0 == wcs_tso_kw['crpix1'], "x-source location not equal to crpix1"
    assert y0 == wcs_tso_kw['crpix2'], "y-source location not equal to crpix2"
    assert order == orderdet, "grism order not equal to detector order in"
    assert_allclose(x, xin), "x grism location changed"
    assert y == wcs_tso_kw['crpix2'], "y grism location not equal to crpix2"


def test_imaging_frames():
    """Verify the available imaging mode reference frames."""
    wcsobj = create_imaging_wcs()
    available_frames = wcsobj.available_frames
    assert all([a == b for a, b in zip(nircam_imaging_frames, available_frames)]), "available frame mismatch"


@pytest.mark.xfail
def test_imaging_distortion():
    """Verify that the distortion correction round trips."""
    wcsobj = create_imaging_wcs()
    sky_to_detector = wcsobj.get_transform('world', 'detector')
    detector_to_sky = wcsobj.get_transform('detector', 'sky')

    # we'll use the crpix as the simplest reference point
    ra = wcs_wfss_kw['crval1']
    dec = wcs_wfss_kw['crval2']

    x, y = sky_to_detector(ra, dec)
    raout, decout = detector_to_sky(x, y)

    assert_allclose(x, wcs_wfss_kw['crpix1'])
    assert_allclose(y, wcs_wfss_kw['crpix2'])
    assert_allclose(raout, ra), "RA position did not roundtrip distortion"
    assert_allclose(decout, dec), "DEC position did not roundtrip distortion"