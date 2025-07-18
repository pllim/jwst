import sys

import numpy as np
import pytest
from astropy.utils.data import get_pkg_data_filename
from stdatamodels.jwst import datamodels
from stpipe import crds_client

from jwst.stpipe import Pipeline, Step
from jwst.stpipe.tests.steps import PipeWithReference, StepWithReference


def library_function():
    import logging

    log = logging.getLogger(__name__)
    log.info("This is a library function log")


class FlatField(Step):
    """
    An example flat-fielding Step.
    """

    spec = """
        threshold = float(default=0.0)  # The threshold below which to remove
        multiplier = float(default=1.0) # Multiply by this number
    """

    # Load the spec from a file

    def process(self, science, flat):
        self.log.info("Removing flat field")
        self.log.info("Threshold: {0}".format(self.threshold))
        library_function()

        output = datamodels.ImageModel(data=science.data - flat.data)
        return output


class Combine(Step):
    """
    A Step that combines a list of images.
    """

    def process(self, images):
        combined = np.zeros((50, 50))
        for image in images:
            combined += image.data
        return datamodels.ImageModel(data=combined)


class Display(Step):
    """
    A Step to display an image.
    """

    def process(self, image):
        pass


class MultiplyBy2(Step):
    """
    A Step that does the incredibly complex thing of multiplying by 2.
    """

    def process(self, image):
        with datamodels.ImageModel(image) as dm:
            dm2 = datamodels.ImageModel()
            dm2.data = dm.data * 2
            return dm2


class MyPipeline(Pipeline):
    """
    A test pipeline.
    """

    step_defs = {"flat_field": FlatField, "combine": Combine, "display": Display}

    spec = """
    science_filename = input_file()  # The input science filename
    flat_filename = input_file(default=None)     # The input flat filename
    """

    def process(self, *args):
        science = datamodels.open(self.science_filename)
        if self.flat_filename is None:
            self.flat_filename = get_pkg_data_filename(
                "data/flat.fits", package="jwst.stpipe.tests"
            )
        flat = datamodels.open(self.flat_filename)
        calibrated = []
        calibrated.append(self.flat_field.run(science, flat))
        combined = self.combine.run(calibrated)
        self.display.run(combined)
        dm = datamodels.ImageModel(combined)
        self.save_model(dm)
        science.close()
        flat.close()
        return dm


def test_pipeline_from_config_file(tmp_cwd):
    config_file_path = get_pkg_data_filename(
        "steps/python_pipeline.cfg", package="jwst.stpipe.tests"
    )
    pipe = Pipeline.from_config_file(config_file_path)

    assert pipe.flat_field.threshold == 42.0
    assert pipe.flat_field.multiplier == 2.0

    pipe.run()


def test_pipeline_python(tmp_cwd):
    steps = {"flat_field": {"threshold": 42.0}}

    pipe = MyPipeline(
        "MyPipeline",
        steps=steps,
        science_filename=get_pkg_data_filename("data/science.fits", package="jwst.stpipe.tests"),
        flat_filename=get_pkg_data_filename("data/flat.fits", package="jwst.stpipe.tests"),
        output_file="python.fits",
    )

    assert pipe.flat_field.threshold == 42.0
    assert pipe.flat_field.multiplier == 1.0

    pipe.run()


def test_prefetch(tmp_cwd, monkeypatch):
    """Test prefetching"""

    # Setup mock to crds to flag if the call was made.
    class MockGetRef:
        called = False

        def mock(self, parameters, reference_file_types, observatory):
            if "flat" in reference_file_types:
                self.called = True
            result = {reftype: "N/A" for reftype in reference_file_types}
            return result

        __call__ = mock

    mock_get_ref = MockGetRef()
    monkeypatch.setattr(crds_client, "get_multiple_reference_paths", mock_get_ref)

    # Create some data
    model = datamodels.ImageModel((19, 19))
    model.meta.instrument.name = "NIRCAM"
    model.meta.instrument.detector = "NRCA1"
    model.meta.instrument.filter = "F070W"
    model.meta.instrument.pupil = "CLEAR"
    model.meta.observation.date = "2019-01-01"
    model.meta.observation.time = "00:00:00"

    # Run the pipeline with prefetch set.
    StepWithReference.prefetch_references = True
    PipeWithReference.call(model)
    assert mock_get_ref.called

    # Now run with prefetch unset.
    mock_get_ref.called = False
    StepWithReference.prefetch_references = False
    PipeWithReference.call(model)
    assert not mock_get_ref.called


def test_pipeline_from_cmdline_cfg(tmp_cwd):
    args = [
        get_pkg_data_filename("steps/python_pipeline.cfg", package="jwst.stpipe.tests"),
        "--steps.flat_field.threshold=47",
    ]

    pipe = Step.from_cmdline(args)

    assert pipe.flat_field.threshold == 47.0
    assert pipe.flat_field.multiplier == 2.0

    pipe.run()


def test_pipeline_from_cmdline_class(tmp_cwd):
    science_filename = get_pkg_data_filename("data/science.fits", package="jwst.stpipe.tests")
    args = [
        "jwst.stpipe.tests.test_pipeline.MyPipeline",
        f"--science_filename={science_filename}",
        "--output_file=output.fits",
        "--steps.flat_field.threshold=47",
    ]

    pipe = Step.from_cmdline(args)

    assert pipe.flat_field.threshold == 47.0
    assert pipe.flat_field.multiplier == 1.0

    pipe.run()


def test_pipeline_commandline_invalid_args():
    from io import StringIO

    flat_filename = get_pkg_data_filename("data/flat.fits", package="jwst.stpipe.tests")
    args = [
        "jwst.stpipe.tests.test_pipeline.MyPipeline",
        # The file_name parameters are *required*, and one of them
        # is missing, so we should get a message to that effect
        # followed by the commandline usage message.
        f"--flat_filename={flat_filename}",
        "--steps.flat_field.threshold=47",
    ]

    sys.stdout = buffer = StringIO()

    with pytest.raises(ValueError):
        Step.from_cmdline(args)

    help = buffer.getvalue()
    assert "Multiply by this number" in help
