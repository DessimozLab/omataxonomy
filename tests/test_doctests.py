import doctest

import pytest

import omataxonomy


def test_core_doctests():
    failures, _ = doctest.testmod(omataxonomy)
    assert failures == 0


def test_session_example():
    failures, _ = doctest.testfile("../docs/examples/session.py")
    assert failures == 0


@pytest.mark.parametrize("page", ["background.md", "tutorial.md"])
def test_documentation_examples(page):
    failures, _ = doctest.testfile(f"../docs/{page}")
    assert failures == 0