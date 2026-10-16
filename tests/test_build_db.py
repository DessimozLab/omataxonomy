import io
import os
import tarfile
import tempfile
import unittest
import unittest.mock
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

    @patch("omataxonomy.build_db.update_mnemonic_codes")
    def test_accepts_relative_paths(self, _):
        from omataxonomy import cwd
        write_taxdump(self.targz, NAMES, NODES)
        with cwd(self.tmpdir.name):
            update_db("taxonomy.sqlite", "omatax.tar.gz")
        self.assertTrue(os.path.exists(self.dbfile))


class DownloadNcbiTaxonomyTest(unittest.TestCase):
    def test_rejects_corrupted_download(self):
        from omataxonomy import cwd
        from omataxonomy.build_db import download_ncbi_taxonomy

        def fake_download(url, filepath):
            with open(filepath, "wb") as fh:
                fh.write(b"truncated")

        md5_response = unittest.mock.Mock(text="0123456789abcdef0123456789abcdef  taxdump.tar.gz\n")
        with tempfile.TemporaryDirectory() as tmpdir, cwd(tmpdir), \
                patch("omataxonomy.build_db.requests.get", return_value=md5_response), \
                patch("omataxonomy.build_db.download_file", side_effect=fake_download):
            with self.assertRaisesRegex(IOError, "md5 checksum"):
                download_ncbi_taxonomy()


class LocalSpeclistTest(unittest.TestCase):
    def test_reads_local_speclist(self):
        from omataxonomy.mnemonic import iter_mnemonic_species_codes
        content = ("header\n"
                   "_____ _ _______  ____________________________________________________________\n"
                   "HUMAN E    9606: N=Homo sapiens\n"
                   "                 C=Human\n"
                   "MOUSE E   10090: N=Mus musculus\n")
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, "speclist.txt")
            with open(fname, "w") as fh:
                fh.write(content)
            self.assertEqual([("HUMAN", "9606"), ("MOUSE", "10090")],
                             list(iter_mnemonic_species_codes(fname)))
