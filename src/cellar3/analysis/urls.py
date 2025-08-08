from django.urls import path

from cellar3.analysis.interaction.views import InteractionView
from src.cellar3.analysis.differential.exporters.views import DiffExpressionExportView
from src.cellar3.analysis.differential.plots.views import ViolinPlotView, DotPlotView, UmapPlotView, VolcanoPlotView
from src.cellar3.analysis.functional.exporters.views import FunctionalExportView
from src.cellar3.analysis.functional.plots.views import BarPlotView, NetworkPlotView
from src.cellar3.analysis.pathways.views import PathwaysPlotView

urlpatterns = [
    path('expression/diff/<str:dataset_id>/', VolcanoPlotView.as_view(), name='diff_expressions'),
    path('expression/diff/<str:dataset_id>/export/<str:export_format>', DiffExpressionExportView.as_view(),
         name='dge_export'),
    path('expression/violin/<str:dataset_id>/', ViolinPlotView.as_view(), name='violin_plot'),
    path('expression/dot/<str:dataset_id>/', DotPlotView.as_view(), name='dot_plot'),
    path('expression/umap/<str:dataset_id>/', UmapPlotView.as_view(), name='umap_plot'),

    path('expression/functional/bar/<str:dataset_id>/', BarPlotView.as_view(), name='functional_bar'),
    path('expression/functional/network/<str:dataset_id>/', NetworkPlotView.as_view(), name='functional_network'),
    path('expression/functional/<str:dataset_id>/export/<str:export_format>', FunctionalExportView.as_view(),
         name='functional_export'),
    path('expression/pathways/<str:dataset_id>/', PathwaysPlotView.as_view(), name='pathways'),
    path('interaction/<str:dataset_id>/', InteractionView.as_view(), name='interaction')
]
