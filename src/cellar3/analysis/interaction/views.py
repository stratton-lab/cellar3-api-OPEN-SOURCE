import logging

from rest_framework.response import Response
from rest_framework.views import APIView

from cellar3.datasets.dataset import Dataset
from cellar3.exceptions import UserException

logger = logging.getLogger('cellar.analysis.interaction')


class InteractionView(APIView):
    """
    Returns the data necessary to draw the Chord diagram for cell-cell interactions.
    """

    def post(self, request, dataset_id):
        """
        Accepts as input the name of source cellType and P-value threshold
        Returns a response of type [{source:'CellType1', target:'CellType2', value: 20}]
        :param request:
        :param dataset_id:
        :return:
        """
        try:
            json_data = request.data
            source = json_data.get('source')
            targets = json_data.get('targets')
            interaction_type = json_data.get('type')
            dataset = Dataset(dataset_id)
            interactions = dataset.get_interactions(source, targets, interaction_type)
            return Response(interactions)
        except Exception as e:
            logger.exception(e)
            raise UserException(f'Could not retrieve cell-cell interaction data: {e}')
