import logging

from django.contrib.admin import AdminSite
from django.contrib.auth.admin import GroupAdmin, UserAdmin
from django.contrib.auth.models import User, Group

from src.cellar3.manage.admin import DatasetInfoAdmin
from src.cellar3.models import DatasetMeta, Submission
from src.cellar3.submissions.admin import SubmissionAdmin

logger = logging.getLogger('cellar.admin.manage')


class MyAdminSite(AdminSite):
    site_header = 'SingloCell Admin'
    site_title = 'SingloCell Admin'
    index_title = 'SingloCell Admin'


admin_site = MyAdminSite(name='myadmin')

admin_site.register(User, UserAdmin)
admin_site.register(Group, GroupAdmin)
admin_site.register(DatasetMeta, DatasetInfoAdmin)
admin_site.register(Submission, SubmissionAdmin)
