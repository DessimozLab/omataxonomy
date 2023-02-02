import unittest
from omataxonomy import Taxonomy


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

