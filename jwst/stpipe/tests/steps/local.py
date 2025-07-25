import logging

from jwst.stpipe import Step

log = logging.getLogger("FOO")


class DummyStep(Step):
    """Placeholder step."""

    spec = """
    foo = string()
    """

    def process(self, *args):  # noqa: D102
        from stdatamodels.jwst.datamodels import ImageModel

        log.info("Default logger")
        log.debug("Default logger")

        self.log.info(f"Foo: {self.foo}")

        self.log.debug("Debug!!!")

        return ImageModel(args[0])
