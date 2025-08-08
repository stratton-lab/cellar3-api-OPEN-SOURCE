import os

from django.conf import settings
from django.core.files.storage import FileSystemStorage
from django.db import models
from django.db.models import FileField


class DatasetMeta(models.Model):
    id = models.CharField(max_length=100, primary_key=True, help_text="Unique identifier of the dataset")
    public = models.BooleanField(default=False, null=False,
                                 help_text="If not set to true, the dataset will only be available if SHOW_PRIVATE_DATASETS is set to True in settings.")
    name = models.CharField(max_length=500, help_text="User friendly name of the dataset")
    type = models.CharField(max_length=100)
    species = models.CharField(max_length=100, help_text="Used for functional analysis.")
    tissue = models.CharField(max_length=100)
    cells = models.IntegerField(help_text="Number of cells in the dataset.")
    categories = models.JSONField(default=list,
                                  help_text="Categories in the explore view of the home page (ex: Transcriptomics)")
    description = models.TextField(max_length=5000)

    file = models.CharField(max_length=100, help_text="Name the h5ad file containing the dataset.")

    image = FileField(
        storage=FileSystemStorage(
            location=os.path.join(settings.MEDIA_ROOT, 'datasets'),
            base_url=os.path.join(settings.MEDIA_URL, 'datasets')
        ), max_length=100)

    linksPublications = models.JSONField(default=list, null=True, blank=True)
    linksDatasets = models.JSONField(default=list, null=True, blank=True)
    keywords = models.JSONField(default=list, null=True, blank=True)

    info = models.JSONField(default=dict,
                            help_text="Maps a standard property name (ex: sample) to the obs column name in the dataset.")
    infoDefault = models.JSONField(default=dict, null=True, blank=True,
                                   help_text="Default values for info columns missing from info mapping.")
    labelsRemap = models.JSONField(default=dict, null=True, blank=True,
                                   help_text="Keys correspond to column names in dataset. Each key contains mapping between original and desired value.")
    embeddings = models.JSONField(default=list, null=True, blank=True,
                                  help_text="Matrix containing 2D cell coordinates (x, y).")
    groups = models.JSONField(default=list, help_text="Columns used to color different points of the data.")
    pseudotime = models.JSONField(default=list, null=True, blank=True)

    maintainer = models.JSONField(default=dict)
    hash = models.CharField(max_length=1000, null=True, blank=True)

    @property
    def maintainer_name(self):
        return self.maintainer.get('name', 'N/A')

    class Meta:
        db_table = 'cellar3_datasets'
        verbose_name = "Dataset"
        verbose_name_plural = 'Datasets'

    def __str__(self):
        return self.name

    def save(self, *args, **kwargs):
        # Using id for image name, with original image extension
        if self.image:
            image_name = f"{self.id}{os.path.splitext(self.image.name)[1]}"
            self.image.name = image_name

        super(DatasetMeta, self).save(*args, **kwargs)


class SubmissionStatus(models.TextChoices):
    RECEIVED = 'received', 'Received'
    ACCEPTED = 'accepted', 'Accepted'
    REJECTED = 'rejected', 'Rejected'


class ConversionStatus(models.TextChoices):
    INITIAL = 'initial', 'Initial'
    JOB_SENT = 'sent', 'Job Submitted'  # Job sent to the conversion script
    JOB_RECEIVED = 'received', 'Job Received'
    RUNNING = 'running', 'Running'
    SUCCESS = 'success', 'Success'
    ERROR = 'error', 'Error'


def get_submission_upload_path(instance, filename):
    # Use the instance's ID as the subfolder name
    return os.path.join(str(instance.id), filename)


class Submission(models.Model):
    """
    Users can submit datasets to the website.
    """
    id = models.CharField(max_length=300, primary_key=True)
    submission_date = models.DateTimeField(auto_now_add=True,
                                           help_text="When the user has filed the submission inquiry.")
    maintainer_email = models.EmailField(max_length=500, help_text="Email of the person who filed the submission.")
    maintainer_name = models.CharField(max_length=500)
    decision_date = models.DateTimeField(null=True, blank=True,
                                         help_text="Date when this submission was Accepted or Rejected.")
    status = models.CharField(max_length=20, choices=SubmissionStatus.choices, default=SubmissionStatus.RECEIVED)
    conversion_status = models.CharField(max_length=20, choices=ConversionStatus.choices,
                                         default=ConversionStatus.INITIAL)
    error = models.CharField(max_length=5000, null=True, blank=True)
    data = models.JSONField(default=dict, null=True, blank=True,
                            help_text="Please review the meta information provided by the user. WARNING: You will not be able to modify the dataset ID ('id' field) after this step.")
    rds = models.FileField(storage=FileSystemStorage(location=settings.SUBMISSIONS_FOLDER),
                           upload_to=get_submission_upload_path, null=True, help_text="RDS file uploaded by the user.")
    h5ad = models.FileField(storage=FileSystemStorage(location=settings.SUBMISSIONS_FOLDER),
                            upload_to=get_submission_upload_path, null=True,
                            help_text="H5AD file generated by the Conversion Service from the RDS file uploaded by the user.")

    @property
    def dataset_id(self):
        return self.data.get('id', 'N/A')

    class Meta:
        db_table = 'cellar3_submissions'
