from django.test import override_settings, SimpleTestCase

from cellar3.analysis.pathways.kegg import KEGG
from cellar3.analysis.pathways.parsers.pathway import Pathway


@override_settings(PATHWAYS_FOLDER='tests/unittests/data/pathways')
class KEGGTests(SimpleTestCase):

    @classmethod
    def setUpClass(cls):
        super(KEGGTests, cls).setUpClass()
        cls.kegg = KEGG()

    def test_term2id(self):
        self.assertEquals(self.kegg.term2id('Metabolic pathways'), '01100')
        self.assertEquals(self.kegg.term2id('01100'), '01100')

    def test_get_pathway(self):
        print(Pathway.get_pathway_path(species_code='mmu', pathway_id='01100'))
        self.assertDictEqual(
            self.kegg.get_pathway(species='Mouse', pathway='Metabolic pathways', upregulated=[], downregulated=[]),
            {'id': 'test-pathway-1'})
