import tempfile
from typing import List, Dict

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy.sparse import lil_matrix


class AnnDataBuilder:

    def __init__(self, genes: List[str], cells: List[str]):
        self.adata = sc.AnnData(
            X=lil_matrix(np.zeros((len(cells), len(genes))), dtype=np.float64),
            obs=pd.DataFrame(index=cells),
            var=pd.DataFrame(index=genes)
        )

    def set_expression(self, gene: str, cell: str, value: float) -> 'AnnDataBuilder':
        gene_idx = self.adata.var_names.get_loc(gene)
        cell_idx = self.adata.obs_names.get_loc(cell)
        self.adata.X[cell_idx, gene_idx] = value
        return self

    def set_expressions(self, gene: str, cells_expressions: Dict[str, float]) -> 'AnnDataBuilder':
        for cell, expression in cells_expressions.items():
            self.set_expression(gene, cell, expression)
        return self

    def set_cell_obs(self, obs_name: str, cells_observations: Dict[str, str]) -> 'AnnDataBuilder':
        # Check if the column exists, if not create it with default values (e.g., None or np.nan)
        if obs_name not in self.adata.obs:
            self.adata.obs[obs_name] = None

        for cell_name, value in cells_observations.items():
            if cell_name in self.adata.obs.index:
                self.adata.obs.at[cell_name, obs_name] = value
            else:
                print(f"Warning: Cell '{cell_name}' not found in the .obs index.")
        return self

    def set_2d_embedding(self, embedding_name: str, embedding: List[List[float]]) -> 'AnnDataBuilder':
        self.adata.obsm[embedding_name] = np.array(embedding)
        return self

    def build(self):
        self.adata.X = self.adata.X.tocsc()
        return self._as_filebacked()

    def _as_filebacked(self) -> AnnData:
        temp_file = tempfile.NamedTemporaryFile(suffix='.h5ad', delete=False)
        temp_filename = temp_file.name
        temp_file.close()
        try:
            # Save the mock data to the temporary file
            self.adata.write_h5ad(temp_filename)
            return anndata.read_h5ad(temp_filename, backed='r')
        finally:
            pass  # TODO remove tmp file


def get_group(adata, cell: str) -> str:
    return adata.obs.at[cell, 'group']
