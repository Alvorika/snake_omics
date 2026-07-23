"""Shared readers for inspected expression matrices.

The ingestion and QC components use this public helper so that 10x format
selection, Gene Expression feature selection, and integer-count validation
have one implementation.
"""

import warnings
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def read_10x_count_matrix(matrix: dict[str, Any]) -> ad.AnnData:
    """Read one manifest matrix record as a validated Gene Expression matrix."""
    matrix_path_value = matrix.get("selected_path")
    matrix_format = matrix.get("selected_format")
    if not matrix_path_value or not matrix_format:
        raise ValueError("Inspected matrix record has no selected path or format")
    if matrix.get("matrix_semantics") != "raw_counts":
        raise ValueError("Expression matrix must declare matrix_semantics='raw_counts'")

    matrix_path = Path(matrix_path_value)
    if matrix_format == "10x_h5":
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="Variable names are not unique.*",
                category=UserWarning,
            )
            adata = sc.read_10x_h5(matrix_path, gex_only=True)
        gene_symbols = adata.var_names.astype(str).to_numpy()
        if "gene_ids" not in adata.var.columns:
            raise ValueError(f"10x HDF5 matrix has no gene_ids field: {matrix_path}")
        gene_ids = adata.var["gene_ids"].astype(str).to_numpy()
    elif matrix_format == "10x_mtx":
        adata = sc.read_10x_mtx(
            matrix_path,
            var_names="gene_ids",
            make_unique=False,
            gex_only=True,
            cache=False,
        )
        gene_ids = adata.var_names.astype(str).to_numpy()
        if "gene_symbols" in adata.var.columns:
            gene_symbols = adata.var["gene_symbols"].astype(str).to_numpy()
        else:
            gene_symbols = gene_ids.copy()
    else:
        raise ValueError(f"Unsupported inspected matrix format: {matrix_format}")

    adata.obs_names = adata.obs_names.astype(str)
    adata.obs_names.name = "barcode"
    adata.var["gene_symbol"] = gene_symbols
    adata.var_names = pd.Index(gene_ids, name="gene_id")
    if not adata.obs_names.is_unique:
        raise ValueError(f"Expression matrix contains duplicate barcodes: {matrix_path}")
    if adata.obs_names.isna().any() or (adata.obs_names == "").any():
        raise ValueError(f"Expression matrix contains missing barcodes: {matrix_path}")
    if not adata.var_names.is_unique:
        raise ValueError(f"Expression matrix contains duplicate gene IDs: {matrix_path}")

    if sparse.issparse(adata.X):
        adata.X = adata.X.tocsr()
        adata.X.eliminate_zeros()
        values = adata.X.data
    else:
        values = np.asarray(adata.X)
    if not np.isfinite(values).all():
        raise ValueError("Expression matrix contains non-finite values")
    if np.any(values < 0):
        raise ValueError("Expression matrix contains negative values")
    if not np.allclose(values, np.rint(values)):
        raise ValueError("Expression matrix is not an integer-count matrix")
    return adata
