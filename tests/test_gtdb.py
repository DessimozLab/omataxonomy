import unittest
from unittest.mock import patch

from omataxonomy import Taxonomy
from omataxonomy.gtdb import assign_taxids, hash_to_int


def fake_hash(table):
    return lambda val: table.get(val, hash_to_int(val))


class AssignTaxidsTest(unittest.TestCase):
    def test_without_collisions_taxid_is_hash(self):
        keys = ["s__Nitrosococcus oceani", "000012805", "f__Nitrosococcaceae"]
        self.assertEqual({k: hash_to_int(k) for k in keys}, assign_taxids(keys))

    def test_alphabetically_first_key_keeps_hash(self):
        table = {"b": -5, "a": -5, "c": -5, "b#1": -7, "c#1": -8}
        with patch("omataxonomy.gtdb.hash_to_int", side_effect=fake_hash(table)):
            ids = assign_taxids(["c", "b", "a"])
        self.assertEqual({"a": -5, "b": -7, "c": -8}, ids)

    def test_reassigned_taxid_avoids_used_values(self):
        # "b#1" collides with the hash of "x", "c#1" with the reassigned "b#2"
        table = {"a": -5, "b": -5, "c": -5, "x": -7, "b#1": -7, "b#2": -8, "c#1": -8, "c#2": -9}
        with patch("omataxonomy.gtdb.hash_to_int", side_effect=fake_hash(table)):
            ids = assign_taxids(["a", "b", "c", "x"])
        self.assertEqual({"a": -5, "b": -8, "c": -9, "x": -7}, ids)

    def test_independent_of_input_order(self):
        table = {"a": -5, "b": -5, "c": -5}
        with patch("omataxonomy.gtdb.hash_to_int", side_effect=fake_hash(table)):
            self.assertEqual(assign_taxids(["a", "b", "c"]), assign_taxids(["c", "a", "b", "a"]))


class CollisionsInReleaseTest(unittest.TestCase):
    """Keys that collide in the pinned GTDB release (see conftest.py)."""

    def test_colliding_genomes_are_distinct(self):
        self.assertEqual(hash_to_int("004212075"), hash_to_int("034420595"))
        tax = Taxonomy()
        names = ["GB_GCA_004212075.1", "RS_GCF_034420595.1"]
        taxids = tax.get_name_translator(names)
        self.assertEqual([hash_to_int("004212075")], taxids["GB_GCA_004212075.1"])
        self.assertNotEqual(taxids["GB_GCA_004212075.1"], taxids["RS_GCF_034420595.1"])
        lineages = tax.get_name_lineage(names)
        self.assertIn("c__Methanoliparia", lineages["GB_GCA_004212075.1"])
        self.assertIn("f__Enterobacteriaceae", lineages["RS_GCF_034420595.1"])

    def test_species_colliding_with_genome_keeps_its_genomes(self):
        self.assertEqual(hash_to_int("949080345"), hash_to_int("s__Fodinicurvata halophila"))
        tax = Taxonomy()
        lineage = tax.get_name_lineage(["RS_GCF_042679585.1"])["RS_GCF_042679585.1"]
        self.assertEqual(["s__Fodinicurvata halophila", "RS_GCF_042679585.1"], lineage[-2:])
        self.assertIn("f__Lachnospiraceae", tax.get_name_lineage(["GB_GCA_949080345.1"])["GB_GCA_949080345.1"])
