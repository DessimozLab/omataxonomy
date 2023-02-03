from .query import Taxonomy, EnvReleaseTaxonomy, is_taxadb_up_to_date
from .build_db import build_combined_tarball
from contextlib import contextmanager
from importlib import metadata
import os

__version__ = metadata.version(__package__)

@contextmanager
def cwd(path):
    oldpwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(oldpwd)


