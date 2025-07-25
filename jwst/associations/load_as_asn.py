"""Treat various objects as Associations."""

from functools import partial
from pathlib import Path

from jwst.associations import Association, AssociationRegistry, libpath
from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase
from jwst.associations.load_asn import load_asn

__all__ = [
    "LoadAsAssociation",
    "LoadAsLevel2Asn",
]


DEFAULT_NAME = "singleton"
DEFAULT_ASN_META = {"program": DEFAULT_NAME, "target": DEFAULT_NAME, "asn_pool": DEFAULT_NAME}


class LoadAsAssociation(dict):
    """
    Read in or create an association.

    Notes
    -----
    This class is normally not instantiated.
    the `load` method should be used as the factory
    method to read an association or create one from
    a string or `Datamodel` object, or a list of such
    objects.
    """

    @classmethod
    def load(
        cls,
        obj,
        meta=DEFAULT_ASN_META,
        registry=AssociationRegistry,
        rule=Association,
        product_name_func=None,
    ):
        """
        Load object and return an association of it.

        Parameters
        ----------
        obj : Association, str, Datamodel, [str[,...]], [Datamodel[,...]]
            The obj to return as an association.

        registry : AssociationRegistry
            The registry to use to load an association file with.

        rule : Association
            The rule to use if an association needs to be created.

        product_name_func : func
            A function, when given the argument of `obj`, or
            if `obj` is a list, each item in `obj`, returns
            a string that will be used as the product name in
            the association.

        Returns
        -------
        association : Association
            An association created using given obj.

        Notes
        -----
        Along with the attributes belonging to a Level2 association, the
        filename is added here, if such a file was passed in. Otherwise
        a default value is given.
        """
        try:
            with Path(obj).open() as fp:
                pure_asn = load_asn(fp, registry=registry)
        except Exception:
            if not isinstance(obj, list):
                obj = [obj]
            asn = asn_from_list(obj, rule=rule, meta=meta, product_name_func=product_name_func)
            asn.filename = DEFAULT_NAME
        else:
            asn = rule()
            asn.update(pure_asn)
            asn.filename = obj

        return asn


class LoadAsLevel2Asn(LoadAsAssociation):
    """Read in or create a Level2 association."""

    @classmethod
    def load(cls, obj, basename=None):
        """
        Open object and return a Level2 association of it.

        Parameters
        ----------
        obj : Association, str, Datamodel, [str[,...]], [Datamodel[,...]]
            The obj to return as an association.

        basename : str
            If specified, use as the basename, with an index appended.

        Returns
        -------
        association : DMSLevel2bBase
            An association created using given obj.

        Notes
        -----
        Along with the attributes belonging to a Level2 association, the
        filename is added here, if such a file was passed in. Otherwise
        a default value is given.
        """
        product_name_func = cls.model_product_name
        if basename is not None:
            product_name_func = partial(cls.name_with_index, basename)

        # if the input string is a FITS file create an asn and return
        if isinstance(obj, str):
            file_name = Path(obj).name
            file_ext = Path(obj).suffix

            if file_ext == ".fits":
                items = [(obj, "science")]
                asn = asn_from_list(
                    items,
                    product_name=file_name,
                    rule=DMSLevel2bBase,
                    with_exptype=True,
                    meta={"asn_pool": "singleton"},
                )
                return asn

        asn = super(LoadAsLevel2Asn, cls).load(
            obj,
            registry=AssociationRegistry(
                definition_files=[libpath() / "rules_level2b.py"], include_default=False
            ),
            rule=DMSLevel2bBase,
            product_name_func=product_name_func,
        )
        return asn

    @staticmethod
    def model_product_name(model, _idx):
        """
        Produce a model product name based on the model.

        Parameters
        ----------
        model : DataModel
            The model to get the name from
        _idx : int
            The parent method is sometimes passed an index,
            which this method ignores.

        Returns
        -------
        product_name : str
            The basename of filename from the model
        """
        return Path(model.meta.filename).stem

    @staticmethod
    def name_with_index(basename, idx):
        """
        Produce a name with the basename and index appended.

        Parameters
        ----------
        basename : str
            The base of the file name
        idx : int
            The current index of the added item.

        Returns
        -------
        product_name : str
            The concatenation of basename, '_', idx

        Notes
        -----
        If the index is less than or equal to 1, no appending is done.
        """
        basename = Path(basename).stem

        if idx > 1:
            basename = basename + "_" + str(idx)
        return basename
