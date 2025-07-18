"""Unit tests for Residual Fringe Correction step interface."""

import numpy as np
import pytest
from stdatamodels.jwst import datamodels

from jwst.residual_fringe import ResidualFringeStep, residual_fringe


@pytest.fixture(scope="function")
def miri_image():
    image = datamodels.IFUImageModel((20, 20))
    image.meta.instrument.name = "MIRI"
    image.meta.instrument.detector = "MIRIFULONG"
    image.meta.exposure.type = "MIR_MRS"
    image.meta.instrument.channel = "12"
    image.meta.instrument.band = "SHORT"
    image.meta.filename = "test_miri.fits"

    rng = np.random.default_rng(42)
    image.data = rng.random((20, 20))

    return image


def test_bad_ignore_regions(tmp_cwd, miri_image):
    # testing the ignore_regions_min
    # There has to be an equal number of min and max ignore region values
    # --ignore_region_min="4.9,"  --ignore_region_max='5.5,"

    step = ResidualFringeStep()
    step.ignore_region_min = [4.9, 5.7]
    step.ignore_region_max = [5.6, 6.5, 9.0]
    step.skip = False

    # If the number ignore min and max regions is not the same a value error is returned
    with pytest.raises(ValueError):
        step.run(miri_image)


def test_ignore_regions(tmp_cwd, monkeypatch, miri_image, log_watcher):
    # Set some reasonable wavelength regions - these should be read in properly
    watcher = log_watcher("stpipe.ResidualFringeStep", message="Ignoring 2 wavelength regions")

    step = ResidualFringeStep()
    step.ignore_region_min = [4.9, 5.7]
    step.ignore_region_max = [5.6, 6.5]
    step.skip = False

    # monkeypatch the reference file retrieval so step aborts but does
    # not error out for this incomplete input
    monkeypatch.setattr(step, "get_reference_file", lambda *args: "N/A")

    # check for ignore regions log message
    step.run(miri_image)
    watcher.assert_seen()


def test_fringe_flat_applied(tmp_cwd, miri_image):
    miri_image.meta.cal_step.fringe = "SKIPPED"
    residual_fringe_reference_file = None
    regions_reference_file = None
    save_intermediate_results = False
    transmission_level = 2
    ignore_regions = {}
    pars = {
        "save_intermediate_results": save_intermediate_results,
        "transmission_level": transmission_level,
    }

    rfc = residual_fringe.ResidualFringeCorrection(
        miri_image, residual_fringe_reference_file, regions_reference_file, ignore_regions, **pars
    )

    # test that the fringe flat step has to be already run
    # on the data before running residual fringe step
    with pytest.raises(residual_fringe.NoFringeFlatError):
        rfc.do_correction()


def test_rf_step_wrong_input_type():
    model = datamodels.ImageModel()
    with pytest.raises(TypeError, match="Failed to process input type"):
        ResidualFringeStep.call(model, skip=False)


def test_rf_step_wrong_exptype(miri_image, log_watcher):
    model = miri_image
    model.meta.exposure.type = "NRS_IFU"

    watcher = log_watcher("stpipe.ResidualFringeStep", message="only for MIRI MRS")
    result = ResidualFringeStep.call(model, skip=False)
    assert result.meta.cal_step.residual_fringe == "SKIPPED"
    watcher.assert_seen()
