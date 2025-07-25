#!/usr/bin/env python
from collections import defaultdict
from pathlib import Path

from stdatamodels.jwst import datamodels

# step imports
from jwst.coron import align_refs_step, klip_step, stack_refs_step
from jwst.datamodels import ModelContainer, ModelLibrary
from jwst.model_blender import ModelBlender
from jwst.outlier_detection import outlier_detection_step
from jwst.resample import resample_step
from jwst.stpipe import Pipeline

__all__ = ["Coron3Pipeline"]


def to_container(model):
    """
    Convert a CubeModel to a ModelContainer of ImageModels for each plane.

    Parameters
    ----------
    model : CubeModel
        The input model to convert

    Returns
    -------
    container : ModelContainer
        The container of ImageModels
    """
    container = ModelContainer()
    for plane in range(model.shape[0]):
        image = datamodels.ImageModel()
        for attribute in [
            "data",
            "dq",
            "err",
            "zeroframe",
            "area",
            "var_poisson",
            "var_rnoise",
            "var_flat",
        ]:
            try:
                arr = model.getarray_noinit(attribute)
                if arr.ndim == 3:
                    setattr(image, attribute, arr[plane])
                elif arr.ndim == 2:
                    setattr(image, attribute, arr)
                else:
                    msg = f"Unexpected array shape for extension {attribute}: {arr.shape}"
                    raise ValueError(msg)
            except AttributeError:
                pass
        image.update(model)
        try:
            image.meta.wcs = model.meta.wcs
        except AttributeError:
            pass
        container.append(image)
    return container


class Coron3Pipeline(Pipeline):
    """
    Class for defining Coron3Pipeline.

    Coron3Pipeline: Apply all level-3 calibration steps to a
    coronagraphic association of exposures. Included steps are:

    #. stack_refs (assemble reference PSF inputs)
    #. align_refs (align reference PSFs to target images)
    #. klip (PSF subtraction using the KLIP algorithm)
    #. outlier_detection (flag outliers)
    #. resample (image combination and resampling)
    """

    class_alias = "calwebb_coron3"

    spec = """
        suffix = string(default='i2d')
    """  # noqa: E501

    # Define aliases to steps
    step_defs = {
        "stack_refs": stack_refs_step.StackRefsStep,
        "align_refs": align_refs_step.AlignRefsStep,
        "klip": klip_step.KlipStep,
        "outlier_detection": outlier_detection_step.OutlierDetectionStep,
        "resample": resample_step.ResampleStep,
    }

    prefetch_references = False

    # Main processing
    def process(self, user_input):
        """
        Run the coron3 pipeline on the input data.

        Parameters
        ----------
        user_input : str, Level3 Association, or ~jwst.datamodels.JwstDataModel
            The exposure or association of exposures to process
        """
        self.log.info("Starting calwebb_coron3 ...")
        asn_exptypes = ["science", "psf"]

        # Create a DM object using the association table
        input_models = datamodels.open(user_input, asn_exptypes=asn_exptypes)

        # This asn_id assignment is important as it allows outlier detection
        # to know the asn_id since that step receives the cube as input.
        self.asn_id = input_models.asn_table["asn_id"]

        # Store the output file for future use
        self.output_file = input_models.asn_table["products"][0]["name"]

        # Find all the member types in the product
        members_by_type = defaultdict(list)
        prod = input_models.asn_table["products"][0]

        for member in prod["members"]:
            members_by_type[member["exptype"].lower()].append(member["expname"])

        # Set up required output products and formats
        self.outlier_detection.suffix = "crfints"
        self.outlier_detection.mode = "coron"
        self.outlier_detection.save_results = self.save_results
        self.resample.blendheaders = False

        # Save the original outlier_detection.skip setting from the
        # input, because it may get toggled off within loops for
        # processing individual inputs
        skip_outlier_detection = self.outlier_detection.skip

        # Extract lists of all the PSF and science target members
        psf_files = members_by_type["psf"]
        targ_files = members_by_type["science"]

        # Make sure we found some PSF and target members
        if len(psf_files) == 0:
            err_str1 = "No reference PSF members found in association table."
            self.log.error(err_str1)
            self.log.error("Calwebb_coron3 processing will be aborted")
            return

        if len(targ_files) == 0:
            err_str1 = "No science target members found in association table"
            self.log.error(err_str1)
            self.log.error("Calwebb_coron3 processing will be aborted")
            return

        for member in psf_files + targ_files:
            self.prefetch(member)

        # Assemble all the input psf files into a single ModelContainer
        psf_models = ModelContainer()
        for i in range(len(psf_files)):
            psf_input = datamodels.CubeModel(psf_files[i])
            psf_models.append(psf_input)

            psf_input.close()

        # Perform outlier detection on the PSFs.
        if not skip_outlier_detection:
            for model in psf_models:
                self.outlier_detection.run(model)
                # step may have been skipped for this model;
                # turn back on for next model
                self.outlier_detection.skip = False
        else:
            self.log.info("Outlier detection skipped for PSFs")

        # Stack all the PSF images into a single CubeModel
        psf_stack = self.stack_refs.run(psf_models)
        psf_models.close()

        # Save the resulting PSF stack
        self.save_model(psf_stack, suffix="psfstack")

        # Call the sequence of steps: outlier_detection, align_refs, and klip
        # once for each input target exposure
        resample_input = ModelContainer()
        model_blender = ModelBlender()
        for target_file in targ_files:
            with datamodels.open(target_file) as target:
                model_blender.accumulate(target)

                # Remove outliers from the target
                if not skip_outlier_detection:
                    target = self.outlier_detection.run(target)
                    # step may have been skipped for this model;
                    # turn back on for next model
                    self.outlier_detection.skip = False

                # Call align_refs
                psf_aligned = self.align_refs.run(target, psf_stack)

                # Save the alignment results
                self.save_model(
                    psf_aligned, output_file=target_file, suffix="psfalign", acid=self.asn_id
                )

                # Call KLIP
                psf_sub = self.klip.run(target, psf_aligned)
                psf_aligned.close()

                # Save the psf subtraction results
                self.save_model(psf_sub, output_file=target_file, suffix="psfsub", acid=self.asn_id)

                # Split out the integrations into separate models
                # in a ModelContainer to pass to `resample`
                for model in to_container(psf_sub):
                    resample_input.append(model)

        # Call the resample step to combine all psf-subtracted target images
        # for compatibility with image3 pipeline use of ModelLibrary,
        # convert ModelContainer to ModelLibrary
        resample_library = ModelLibrary(resample_input, on_disk=False)

        # Output is a single datamodel
        result = self.resample.run(resample_library)

        # Blend the science headers
        try:
            completed = result.meta.cal_step.resample
        except AttributeError:
            self.log.debug("Could not determine if resample was completed.")
            self.log.debug("Presuming not.")

            completed = "SKIPPED"
        if completed == "COMPLETE":
            self.log.debug(f"Blending metadata for {result}")
            model_blender.finalize_model(result)

        try:
            result.meta.asn.pool_name = input_models.asn_pool_name
            result.meta.asn.table_name = Path(user_input).name
        except AttributeError:
            self.log.debug("Cannot set association information on final")
            self.log.debug(f"result {result}")

        # Save the final result
        self.save_model(result, suffix=self.suffix)

        # We're done
        self.log.info("...ending calwebb_coron3")

        return
