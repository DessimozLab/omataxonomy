import doctest
import omataxonomy


def test_core_doctests():
    failures, _ = doctest.testmod(omataxonomy)
    assert failures == 0


def test_session_example():
    failures, _ = doctest.testfile("../docs/examples/session.py")
    assert failures == 0