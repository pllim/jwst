from pathlib import Path

from jwst.datamodels import ModelContainer
from jwst.stpipe import Step
from jwst.wfs_combine import wfs_combine

__all__ = ["WfsCombineStep"]


class WfsCombineStep(Step):
    """Combine pairs of dithered PSF images."""

    class_alias = "calwebb_wfs-image3"

    spec = """
        do_refine = boolean(default=False)
        flip_dithers = boolean(default=True) # change the sign and switch order of images when x offset is negative
        psf_size = integer(default=100)
        blur_size = integer(default=10)
        n_size = integer(default=2)
        suffix = string(default="wfscmb")
    """  # noqa: E501

    def make_output_path(self, basepath, *args, **kwargs):  # noqa:  ARG002
        """
        Filename formatting bypass in stpipe.

        Returns
        -------
        basepath : str
            Basename of the file
        """
        return basepath

    def process(self, input_table):
        """
        Combine image pairs.

        Returns
        -------
        output_container : ModelContainer
            Model container with all combined data
        """
        self.suffix = "wfscmb"
        self.output_use_model = True

        # Load the input ASN table
        asn_table = self.load_as_level3_asn(input_table)

        self.log.info("Using input table: %s", input_table)
        self.log.info("The number of pairs of input files: %g", len(asn_table["products"]))

        output_container = ModelContainer()

        # Process each pair of input images listed in the association table
        for which_set in asn_table["products"]:
            # Get the list of science members in this pair
            science_members = [
                member for member in which_set["members"] if member["exptype"].lower() == "science"
            ]
            infile_1 = science_members[0]["expname"]
            infile_2 = science_members[1]["expname"]
            outfile = which_set["name"]

            # Create the step instance
            wfs = wfs_combine.DataSet(
                infile_1,
                infile_2,
                outfile,
                self.do_refine,
                self.flip_dithers,
                self.psf_size,
                self.blur_size,
                self.n_size,
            )

            # Do the processing
            output_model = wfs.do_all()

            # The DataSet class does not close its resources.  Do that here.
            wfs.input_1.close()
            wfs.input_2.close()

            # Update necessary meta info in the output
            output_model.meta.cal_step.wfs_combine = "COMPLETE"
            output_model.meta.asn.pool_name = asn_table["asn_pool"]
            output_model.meta.asn.table_name = Path(input_table).name
            # format the filename here
            output_model.meta.filename = (
                which_set["name"].format(suffix=self.suffix) + self.output_ext
            )

            output_container.append(output_model)

        # Return the output so it can be tested.
        return output_container
