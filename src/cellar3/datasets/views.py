import logging

from rest_framework.response import Response
from rest_framework.views import APIView

from cellar3.datasets.prod_status import get_datasets_prod_status, assign_prod_status
from src.cellar3.datasets.categories import CATEGORIES_READER
from src.cellar3.datasets.dataset import Dataset
from src.cellar3.datasets.meta_reader import META_READER
from src.cellar3.exceptions import UserException

logger = logging.getLogger('cellar.dataset.views')


class DatasetsSearchView(APIView):
    """
    @todo Implement server-side datasets search.
    """
    pass


class DatasetsView(APIView):
    """
    Returns all available public datasets.
    """

    def get(self, request):
        try:
            return Response({
                'datasets': self._get_public_datasets(),
                'categories': CATEGORIES_READER.get_categories()
            })
        except Exception as e:
            logger.exception(e)
            raise UserException(f'Could not retrieve public datasets: {e}')

    @classmethod
    def _get_public_datasets(cls):
        datasets = META_READER.get_public_datasets()
        assign_prod_status(datasets)
        return datasets


class DatasetView(APIView):
    """
    Returns plot data for the dataset.
    """

    # @method_decorator(cache_page(60 * 15))  # Cache page for 15 minutes
    def get(self, request, dataset_id):
        """
        #TODO Cache the dataset object instead, not the meta.
        :param request:
        :param dataset_id:
        :return:
        """
        try:
            group_by = request.query_params.get('group_by', None)
            embedding_name = request.query_params.get('embedding', None)
            dataset = Dataset(dataset_id).as_json(group_by=group_by, embedding_name=embedding_name)
            assign_prod_status([dataset['meta']])
            return Response(dataset)
        except Exception as e:
            logger.exception(e)
            raise UserException(f'Could not retrieve dataset: {e}')

    # @todo Add a post endpoint, accepting an optional list of point ids used for filtering of returned dataset
    def post(self, request, dataset_id):
        raise UserException('POST Not yet implemented.')


class GeneExpressionView(APIView):
    """
    Returns the expression values of a specific gene, in each cell.
    """

    def get(self, request, dataset_id, gene_name):
        """

        :param request:
        :param dataset_id:
        :param gene_name:
        :return:
        """
        try:
            dataset = Dataset(dataset_id)
            return Response({
                'gene': gene_name,
                'cells': dataset.get_gene_expression(gene_name=gene_name)})
        except Exception as e:
            logger.exception(e)
            raise UserException(f'Could not retrieve gene expression values: {e}')


class CellExpressionView(APIView):
    """
    Returns the expression values for each gene, in a specific cell.
    """

    def get(self, request, dataset_id: str, cell_id: str):
        """
        Returns the expression values for each gene, in a specific cell.

        :param request:
        :param dataset_id:
        :param cell_id:
        :return:
        """
        try:
            dataset = Dataset(dataset_id)
            return Response({
                'cell': cell_id,
                'genes': dataset.get_cell_expression(cell_id=cell_id)
            })
        except Exception as e:
            logger.exception(e)
            raise UserException(f'Could not retrieve cell expression values: {e}')


class GeneSampleExpressionView(APIView):
    """
    Returns the expression value of a specific gene, in a specific cell.
    """

    def get(self, request, dataset_id, gene_name: str, cell_id: str):
        """

        :param request:
        :param dataset_id:
        :param gene_name:
        :param cell_id:
        :return:
        """
        try:
            dataset = Dataset(dataset_id)
            return Response({'expression': dataset.get_expression(cell_id=cell_id, gene_name=gene_name)})
        except Exception as e:
            logger.exception(e)
            raise UserException(f'Could not retrieve expression value: {e}')


class UmapView(APIView):
    """
    Returns information to draw UMAP plot.
    """

    def get(self, request, dataset_id):
        try:
            group_by = request.query_params.get('group_by', None)
            embedding_name = request.query_params.get('embedding', None)
            dataset = Dataset(dataset_id)
            umap = dataset.get_umap(group_by=group_by, embedding_name=embedding_name)
            return Response(umap)
        except Exception as e:
            logger.exception(e)
            raise UserException(f'Could not retrieve UMAP: {e}')


class DatasetInfoView(APIView):
    def get(self, request, dataset_id):
        try:
            return Response(Dataset(dataset_id).get_dataset_info())
        except Exception as e:
            logger.exception(e)
            raise UserException(f'Could not retrieve dataset info: {e}')


class DatasetsStatusView(APIView):
    """
    Returns for each submitted dataset id whether it's public or private.
    """

    def post(self, request):
        datasets = META_READER.get_multiple_dataset_info(request.data)
        return Response({dataset['id']: dataset['public'] for dataset in datasets})


class DatasetsProdStatusView(APIView):
    def get(self, request):
        return Response(get_datasets_prod_status(META_READER.get_public_datasets()))
