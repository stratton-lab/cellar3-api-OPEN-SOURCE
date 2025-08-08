from django.core.management import BaseCommand

from cellar3.datasets.dataset import Dataset


class Command(BaseCommand):
    help = 'Displays gene expression in a cell'

    def add_arguments(self, parser):
        parser.add_argument('dataset_id', type=str, help='ID of the dataset.')
        parser.add_argument('gene_id', type=str, help='Gene symbol')
        parser.add_argument('cell_id', type=str, help='Cell ID')

    def handle(self, *args, **kwargs):
        dataset_id = kwargs['dataset_id']
        gene_name = kwargs['gene_id']
        cell_id = kwargs['cell_id']
        dataset = Dataset(dataset_id=dataset_id)
        expression = dataset.get_expression(gene_name=gene_name, cell_id=cell_id)
        print(f'Expression value of {gene_name} in {cell_id}: {expression}')
