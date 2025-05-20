import os
import unittest
from unittest.mock import patch

from omataxonomy import Taxonomy, EnvReleaseTaxonomy


class TaxonomyTest(unittest.TestCase):

    def test_named_lineage(self):
        tax = Taxonomy()
        lineages = tax.get_name_lineage(['Homo sapiens', "f__Leptotrichiaceae", "Gallus"])
        self.assertIn("Mammalia", lineages['Homo sapiens'])
        self.assertIn("o__Fusobacteriales", lineages["f__Leptotrichiaceae"])
        self.assertIn("Aves", lineages['Gallus'])
        self.assertIn("d__Bacteria", lineages["f__Leptotrichiaceae"])
        self.assertNotIn("Bacteria", lineages["f__Leptotrichiaceae"])

    def test_stable_gtdb_taxids(self):
        expected = {'f__Leptotrichiaceae': [-769982259],
                    'GB_GCA_001515945.1': [-1326969492],
                    's__Moorella thermoacetica': [-1120457739]}
        tax = Taxonomy()
        taxids = tax.get_name_translator(["f__Leptotrichiaceae", "GB_GCA_001515945.1", "s__Moorella thermoacetica"])
        self.assertEqual(taxids, expected)

    def test_lineage_from_acc(self):
        tax = Taxonomy()
        taxid = tax.get_name_translator(['GCA_006226595.1'])
        res = tax.get_name_lineage(['GCA_006226595.1'])['GCA_006226595.1']
        self.assertEqual("root; d__Bacteria; p__Myxococcota_A; c__UBA796; o__UBA796; f__GCA-2862545; g__M1803; s__M1803 sp006226595; GB_GCA_006226595.1".split('; '),
                         res)

    def test_gtdb_common_name(self):
        tax = Taxonomy()
        spname = "GB_GCA_002731275.1"
        taxid = tax.get_name_translator([spname])[spname][0]
        common = tax.get_common_names([taxid])[taxid]
        self.assertEqual("Deltaproteobacteria bacterium NP119", common)


class EnvBasedTaxonomyTest(unittest.TestCase):
    def test_raises_with_env_set_to_inexisting_value(self):
        with self.assertRaises(KeyError):
            tax = EnvReleaseTaxonomy(env_var="XXX")

    @patch('omataxonomy.query.update_db')
    def test_calls_update_db_with_good_args(self, mocked_updated_db):
        os.environ['TEST_LOC'] = "/tmp/test"
        try:
            tax = EnvReleaseTaxonomy(env_var="TEST_LOC")
        except ValueError:
            pass
        self.assertTrue(mocked_updated_db.called)
        mocked_updated_db.assert_called_with("/tmp/test/taxonomy.sqlite", targz_file="/tmp/test/omatax.tar.gz")
