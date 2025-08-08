import logging

from rest_framework.response import Response
from rest_framework.views import APIView

from src.cellar3.analysis.functional.enrich import Enrichir
from src.cellar3.analysis.pathways.kegg import KEGG
from src.cellar3.receipt import Receipt

logger = logging.getLogger('cellar.analysis.pathways')


class PathwaysPlotView(APIView):
    analyzer = Enrichir()
    pathways_provider = KEGG()
    gene_sets = {
        'Human': 'KEGG_2021_Human',
        'Mouse': 'KEGG_2019_Mouse'
    }

    def get(self, request, dataset_id):
        return Response({'status': 'This plot is only available through POST'})

    def post(self, request, dataset_id):
        receipt = Receipt.create()

        # Input
        json_data = request.data
        pathway = json_data.get('pathway')
        species = json_data.get('species')
        upregulated = json_data.get('upregulated')
        downregulated = json_data.get('downregulated')

        # 3. Pathways network with highlighted u and down regulated genes
        data = self.pathways_provider.get_pathway(
            species=species,
            pathway=pathway,
            upregulated=upregulated,
            downregulated=downregulated
        )

        return Response({'pathway': data, 'meta': receipt})
