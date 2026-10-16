import os
from unittest.mock import patch

import pytest

from omataxonomy import Taxonomy
from omataxonomy.query import DB_VERSION

# GTDB release the expected test outputs are based on. Pinning it keeps the
# tests from breaking whenever GTDB publishes a new release. When bumping it,
# the expected outputs in the tests and in docs/examples/session.py need to
# be regenerated.
GTDB_TEST_RELEASE = "232"
TEST_DB = os.path.expanduser(os.path.join(
    "~", ".cache", "omataxonomy", f"test-gtdb-r{GTDB_TEST_RELEASE}-v{DB_VERSION}", "taxonomy.sqlite"))


@pytest.fixture(scope="session", autouse=True)
def pinned_gtdb_release_db():
    """Build a database with the pinned GTDB release (once) and make it the
    default database for all tests and doctests."""
    Taxonomy(dbfile=TEST_DB, gtdb_release=GTDB_TEST_RELEASE)
    with patch("omataxonomy.query.DEFAULT_DB", TEST_DB):
        yield TEST_DB
