import logging
from typing import List, Dict, Union

from src.cellar3.datasets.dataset import Dataset

logger = logging.getLogger('cellar.analysis.expression.umap')


class UMAPPlot:

    def as_plot_data(self, dataset: Dataset, embedding_name: Union[str, None], genes: List[str],
                     meta: Dict[str, any], normalize_to_gene=False) -> Dict:
        pca_x, pca_y = dataset.get_embedding(embedding_name=embedding_name)
        global_max_expression = dataset.get_max_expression()
        logger.info(f'For UMAP plot using embedding {embedding_name} and max expression {global_max_expression}')
        genes_expressions = []
        for gene in genes:
            max_expression = dataset.get_max_gene_expression(gene) if normalize_to_gene else global_max_expression
            genes_expressions.append({
                'gene': gene,
                'expression': dataset.get_normalized_gene_expression(gene, max_expression)
            })

        return {
            'x': Dataset.np2list(pca_x),
            'y': Dataset.np2list(pca_y),
            'genes': genes_expressions
        }
