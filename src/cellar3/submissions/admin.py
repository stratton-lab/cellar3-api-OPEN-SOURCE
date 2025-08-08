import json
import logging

from django import forms
from django.conf import settings
from django.contrib import admin
from django.contrib.admin.widgets import AdminTextareaWidget
from django.forms import ModelForm
from django.http import HttpResponseRedirect
from django.urls import path, reverse
from django.utils import timezone
from django.utils.safestring import mark_safe

from src.cellar3.models import Submission, SubmissionStatus
from src.cellar3.tools import raise_if_not_dev

logger = logging.getLogger('cellar.admin.submission')


class SubmissionForm(ModelForm):
    class Meta:
        model = Submission
        fields = []
        widgets = {
            'data': AdminTextareaWidget()
        }


class PrettyJSONWidget(forms.Widget):
    def render(self, name, value, attrs=None, renderer=None):
        try:
            # Pretty-print the JSON content
            pretty_json = json.dumps(value, indent=4, sort_keys=True)
        except (TypeError, ValueError):
            # If value is not valid JSON, fall back to displaying it as-is
            pretty_json = value
        return mark_safe(f'<pre>{pretty_json}</pre>')


class SubmissionAdmin(admin.ModelAdmin):
    form = SubmissionForm
    change_form_template = "admin/cellar3/submission_change_form.html"
    list_display = ('dataset_id', 'submission_date', 'status', 'conversion_status', 'maintainer_name',)
    ordering = ('-submission_date',)
    search_fields = ('maintainer_name', 'status')
    list_filter = ('status', 'conversion_status')

    fields = ['maintainer_name', 'maintainer_email', 'submission_date', 'data', 'decision_date', 'status',
              'conversion_status', 'error', 'rds', 'h5ad']
    readonly_fields = ['maintainer_name', 'maintainer_email', 'submission_date', 'decision_date',
                       'conversion_status', 'error', 'rds', 'h5ad']

    # TODO Only show Accept / Reject buttons if status is RECEIVED

    def get_urls(self):
        urls = super().get_urls()
        custom_urls = [
            path('accept/<uuid:object_id>/', self.admin_site.admin_view(self.accept), name='accept'),
            path('reject/<uuid:object_id>/', self.admin_site.admin_view(self.reject), name='reject'),
            path('delete_files/<uuid:object_id>/', self.admin_site.admin_view(self.delete_files), name='delete_files')
        ]
        return custom_urls + urls

    def accept(self, request, object_id):
        obj = self.get_object(request, object_id)
        raise_if_not_dev()

        obj.status = SubmissionStatus.ACCEPTED
        obj.decision_date = timezone.now()
        obj.save()

        # Redirect back to the detailed view after the action is performed
        # RDS upload page accessed either on local or PROD
        rds_upload_url = f'{settings.FRONTEND_PROD_BASE_URL}submit-rds-file/{obj.id}'
        message = f'Submission ACCEPTED. Please use this URL to upload RDS file: <a target="_blank" href="{rds_upload_url}">Upload RDS</a>'
        self.message_user(request, mark_safe(message))
        return HttpResponseRedirect(reverse('admin:cellar3_submission_change', args=[object_id]))

    def reject(self, request, object_id):
        obj = self.get_object(request, object_id)
        raise_if_not_dev()
        obj.status = SubmissionStatus.REJECTED
        obj.decision_date = timezone.now()
        obj.save()

        # Redirect back to the detailed view after the action is performed
        self.message_user(request, f"Submission {obj} REJECTED.")
        return HttpResponseRedirect(reverse('admin:cellar3_submission_change', args=[object_id]))

    def delete_files(self, request, object_id):
        raise_if_not_dev()
        from src.cellar3.management.commands.delete_submission_file import Command as DeleteSubmissionFileCommand
        DeleteSubmissionFileCommand().handle(submission_id=object_id)
        self.message_user(request, f"Files deleted.")
        return HttpResponseRedirect(reverse('admin:cellar3_submission_change', args=[object_id]))
