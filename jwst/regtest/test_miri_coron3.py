import pytest

from jwst.regtest.st_fitsdiff import STFITSDiff as FITSDiff
from jwst.stpipe import Step

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


@pytest.fixture(scope="module")
def run_pipeline(rtdata_module, resource_tracker):
    """Run calwebb_coron3 on MIRI 4QPM coronographic data."""
    rtdata = rtdata_module
    rtdata.get_asn("miri/coron/jw01386-c1002_20230109t015044_coron3_00001_asn.json")

    # Run the calwebb_coron3 pipeline on the association
    args = [
        "calwebb_coron3",
        rtdata.input,
    ]
    with resource_tracker.track():
        Step.from_cmdline(args)

    return rtdata


def test_log_tracked_resources_coron3(log_tracked_resources, run_pipeline):
    log_tracked_resources()


@pytest.mark.parametrize("suffix", ["crfints", "psfalign", "psfsub"])
@pytest.mark.parametrize("exposure", ["4", "5"])
def test_miri_coron3_sci_exp(run_pipeline, suffix, exposure, fitsdiff_default_kwargs):
    """Check intermediate results of calwebb_coron3"""
    rtdata = run_pipeline

    output = "jw0138600" + exposure + "001_04101_00001_mirimage_c1002_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_miri_coron3/" + output)

    fitsdiff_default_kwargs["atol"] = 1e-2
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()


@pytest.mark.parametrize("suffix", ["psfstack", "i2d"])
def test_miri_coron3_product(run_pipeline, suffix, fitsdiff_default_kwargs):
    """Check final products of calwebb_coron3"""
    rtdata = run_pipeline

    output = "jw01386-c1002_t001_miri_f1140c-mask1140_" + suffix + ".fits"
    rtdata.output = output
    rtdata.get_truth("truth/test_miri_coron3/" + output)

    fitsdiff_default_kwargs["atol"] = 1e-4
    fitsdiff_default_kwargs["rtol"] = 1e-4
    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
