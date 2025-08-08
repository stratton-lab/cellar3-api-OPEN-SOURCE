from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include

from handlers import handler404json, handler400json, handler500json, forbidden_view
from src.cellar3.admin import admin_site
from src.cellar3.analysis import urls as analysis_urls
from src.cellar3.datasets import urls as dataset_urls
from src.cellar3.datasets.views import DatasetsView
from src.cellar3.submissions import urls as submit_urls
from src.cellar3.tools import is_dev
from src.cellar3.views import IndexView

urlpatterns = [
    path('', IndexView.as_view(), name='index'),
    path('marketplace/', DatasetsView.as_view(), name='marketplace'),
    path('dataset/', include(dataset_urls)),
    path('analysis/', include(analysis_urls)),
    path('submit/', include(submit_urls)),
    path('manage/', admin_site.urls if is_dev() else forbidden_view)
]

if settings.DEBUG:
    # Used to serve media (dataset images) on local deployment.
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

handler404 = handler404json
handler400 = handler400json
handler500 = handler500json
