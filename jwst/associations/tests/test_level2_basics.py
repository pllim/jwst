"""Test basic usage of Level2 associations"""

import re

from astropy.utils.data import get_pkg_data_filename

from jwst.associations import (
    generate,
    load_asn,
)
from jwst.associations.main import Main
from jwst.associations.tests.helpers import (
    combine_pools,
    registry_level2_only,
)

REGEX_LEVEL2A = r"(?P<path>.+)(?P<type>_rate(ints)?)(?P<extension>\..+)"


def from_level2_schema():
    with open(
        get_pkg_data_filename("data/asn_level2.json", package="jwst.associations.tests")
    ) as asn_file:
        asn = load_asn(asn_file)
    return [asn]


def cmd_from_pool(pool_path, args):
    """Run commandline on pool

    Parameters
    ---------
    pool_path: str
        The pool to run on.

    args: [arg(, ...)]
        Additional command line arguments in the form `sys.argv`
    """
    full_args = [
        "--dry-run",
        "-D",
        "-r",
        get_pkg_data_filename("lib/rules_level2b.py", package="jwst.associations"),
        "--ignore-default",
    ]
    full_args.extend(args)
    result = Main(full_args, pool=pool_path)
    return result


def test_level2_productname():
    rules = registry_level2_only()
    pool = combine_pools(
        get_pkg_data_filename("data/pool_002_image_miri.csv", package="jwst.associations.tests")
    )
    for asn in generate(pool, rules):
        for product in asn["products"]:
            science = [member for member in product["members"] if member["exptype"] == "science"]
            assert len(science) == 1
            match = re.match(REGEX_LEVEL2A, science[0]["expname"])
            assert match.groupdict()["path"] == product["name"]
