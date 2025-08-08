import numpy as np
from django.core.management import BaseCommand

from cellar3.datasets.dataset import Dataset


class Command(BaseCommand):
    help = 'Displays gene expression'

    def add_arguments(self, parser):
        # Optional: Define command-line arguments here
        parser.add_argument('dataset_id', type=str, help='ID of the dataset.')
        parser.add_argument('gene_id', type=str, help='Gene symbol')

    def handle(self, *args, **kwargs):
        dataset_id = kwargs['dataset_id']
        gene_name = kwargs['gene_id']
        dataset = Dataset(dataset_id=dataset_id)
        adata = dataset.data
        # adata[adata.obs['cell_type'] == cell_type].copy()
        gene_expression = adata[:, gene_name].X
        non_zero_cells = np.sum(gene_expression > 0)

        total_cells = adata.n_obs
        print(f'Number of cells with non zero expression for gene {gene_name}: {non_zero_cells} / {total_cells}')
        print(f"Minimum expression value for {gene_name}: {gene_expression.min()}")
        print(f"Maximum expression value for {gene_name}: {gene_expression.max()}")
