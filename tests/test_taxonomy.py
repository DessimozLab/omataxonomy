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
                    'GB_GCA_001515945.1': [-1584033363],
                    's__Moorella thermoacetica': [-1120457739]}
        tax = Taxonomy()
        taxids = tax.get_name_translator(["f__Leptotrichiaceae", "GB_GCA_001515945.1", "s__Moorella thermoacetica"])
        self.assertEqual(taxids, expected)


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