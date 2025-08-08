from django.urls import path

from src.cellar3.datasets.views import DatasetView, GeneExpressionView, GeneSampleExpressionView, CellExpressionView, \
    UmapView, DatasetInfoView, DatasetsStatusView, DatasetsProdStatusView

urlpatterns = [
    path('<str:dataset_id>/', DatasetView.as_view(), name='display'),
    path('<str:dataset_id>/umap', UmapView.as_view(), name='umap'),
    path('<str:dataset_id>/expression/gene/<str:gene_name>', GeneExpressionView.as_view(), name='gene_expression'),
    path('<str:dataset_id>/expression/cell/<str:cell_id>', CellExpressionView.as_view(),
         name='sample_expression'),
    path('<str:dataset_id>/expression/value/<str:gene_name>/<str:cell_id>', GeneSampleExpressionView.as_view(),
         name='gene_sample_expression'),
    path('<str:dataset_id>/info', DatasetInfoView.as_view(), name='dataset_info'),
    path('status/public', DatasetsStatusView.as_view(), name='datasets_status_public'),
    path('status/public/prod', DatasetsProdStatusView.as_view(), name='datasets_status_prod')
]
