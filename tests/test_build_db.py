import io
import os
import tarfile
import tempfile
import unittest
from unittest.mock import patch

from omataxonomy.build_db import update_db

NAMES = [("1", "root", "", "scientific name"),
         ("-5", "d__Bacteria", "", "scientific name"),
         ("-7", "s__Foo bar", "", "scientific name")]
NODES = [("1", "1", "no rank"),
         ("-5", "1", "superkingdom"),
         ("-7", "-5", "species")]


def write_taxdump(path, names, nodes):
    def add(tar, name, rows):
        data = "".join("\t|\t".join(row) + "\t|\n" for row in rows).encode()
        info = tarfile.TarInfo(name)
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))

    with tarfile.open(path, "w:gz") as tar:
        add(tar, "names.dmp", names)
        add(tar, "nodes.dmp", nodes)
        add(tar, "merged.dmp", [])


class UpdateDbTest(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.dbfile = os.path.join(self.tmpdir.name, "taxonomy.sqlite")
        self.targz = os.path.join(self.tmpdir.name, "omatax.tar.gz")

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_rejects_archive_with_duplicated_taxids(self):
        # an archive built before hash collisions were resolved: two GTDB
        # entries share taxid -7, the second would silently replace the first
        write_taxdump(self.targz,
                      NAMES + [("-7", "s__Other species", "", "scientific name")],
                      NODES + [("-7", "-5", "species")])
        with self.assertRaisesRegex(ValueError, "1 duplicated taxids"):
            update_db(self.dbfile, self.targz)
        self.assertFalse(os.path.exists(self.dbfile))
        self.assertTrue(os.path.exists(self.targz))

    @patch("omataxonomy.build_db.update_mnemonic_codes")  # avoids downloading uniprot speclist
    def test_accepts_archive_without_duplicates(self, _):
        write_taxdump(self.targz, NAMES, NODES)
        update_db(self.dbfile, self.targz)
        self.assertTrue(os.path.exists(self.dbfile))
