from django.contrib import admin
from django.contrib.admin.widgets import AdminTextareaWidget
from django.forms import ModelForm, TextInput, Textarea, Select

from cellar3.manage.fields.groups import GROUPS_FIELD
from src.cellar3.manage.fields.info import INFO_FIELD
from src.cellar3.manage.fields.links import LINKS_FIELD
from src.cellar3.manage.fields.maintainer import MAINTAINER_FIELD
from src.cellar3.manage.fields.pseudotime import PSEUDOTIME_FIELD
from src.cellar3.manage.fields.smart_json import SmartJSONField
from src.cellar3.manage.widgets.dict import JSONDictWidget
from src.cellar3.manage.widgets.list import ListWidget
from src.cellar3.manage.widgets.list_dicts import ListOfDicts
from src.cellar3.manage.widgets.single_list import CategoriesListWidget
from src.cellar3.models import DatasetMeta


class DatasetInfoForm(ModelForm):
    linksPublications = LINKS_FIELD
    linksDatasets = LINKS_FIELD
    embeddings = SmartJSONField(widget=ListOfDicts(fields=['name', 'key', 'image']), required=False)
    groups = GROUPS_FIELD
    info = INFO_FIELD
    infoDefault = SmartJSONField(widget=JSONDictWidget(), required=False)
    pseudotime = PSEUDOTIME_FIELD
    maintainer = MAINTAINER_FIELD

    class Meta:
        model = DatasetMeta
        fields = '__all__'
        widgets = {
            'name': TextInput(attrs={'size': 100}),
            'description': Textarea(attrs={'cols': 100, 'rows': 3}),
            'keywords': ListWidget(attrs={'size': 80}),
            'file': TextInput(attrs={'size': 100}),
            'categories': CategoriesListWidget(),
            'labelsRemap': AdminTextareaWidget(),
            'species': Select(choices=[
                ('Mouse', 'Mouse'),
                ('Human', 'Human')
            ]),
            'type': Select(choices=[
                ('Suspension', 'Suspension'),
                ('Spatial', 'Spatial')
            ])
        }


class DatasetInfoAdmin(admin.ModelAdmin):
    form = DatasetInfoForm
    change_form_template = "admin/cellar3/dataset_change_form.html"
    list_display = ('id', 'public', 'name', 'maintainer_name')
    search_fields = ('name', 'type', 'species', 'tissue')
    list_filter = ('public', 'species', 'type')

    def get_readonly_fields(self, request, obj=None):
        # obj is None when creating a new object
        if obj is None:
            return []  # All fields writable when creating
        return ['id']  # Make 'id' read-only when editing

    def get_fields(self, request, obj=None):
        fields = super().get_fields(request, obj)
        readonly = self.get_readonly_fields(request, obj)
        rest_of_fields = [field for field in fields if field not in readonly]
        return readonly + rest_of_fields
