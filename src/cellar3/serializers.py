from rest_framework.fields import FileField
from rest_framework.serializers import Serializer, ModelSerializer

from src.cellar3.models import DatasetMeta


class UploadSerializer(Serializer):
    file = FileField()

    class Meta:
        fields = ['file']


class DatasetInfoSerializer(ModelSerializer):
    """
    Excludes maintainer fields. Used to retrieve data for the dataset display page.
    """

    class Meta:
        model = DatasetMeta
        fields = ['id', 'public', 'name', 'type', 'species', 'tissue', 'cells', 'categories', 'description', 'file',
                  'image', 'linksPublications', 'linksDatasets', 'keywords', 'info', 'labelsRemap', 'embeddings',
                  'groups', 'pseudotime']


class DatasetInfoPreviewSerializer(ModelSerializer):
    """
    Returns minimal info about a dataset, used to display a list of datasets previews.
    """

    class Meta:
        model = DatasetMeta
        fields = ['id', 'public', 'name', 'type', 'species', 'tissue', 'cells', 'categories', 'description', 'image']


class DatasetInfoFullSerializer(ModelSerializer):
    """
    Serializes into a JSON doc ALL fields of the model.
    """

    class Meta:
        model = DatasetMeta
        fields = '__all__'
