from django.urls import path

from src.cellar3.submissions.dataset.views import SubmitDatasetFileView
from src.cellar3.submissions.meta.views import SubmitDatasetMetaView

urlpatterns = [
    path('meta/', SubmitDatasetMetaView.as_view(), name='submit-meta'),
    path('file/<str:submission_id>/', SubmitDatasetFileView.as_view(), name='submit-file')
]
