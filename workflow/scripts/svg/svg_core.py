"""Small numerical helpers for per-section ROI SVG analysis."""

from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components


HEX_DELTAS = ((0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1))


def build_visium_hex_graph(array_coordinates: np.ndarray) -> sparse.csr_matrix:
    """Build exact native Visium six-neighbour adjacency."""

    coordinates = np.asarray(array_coordinates)
    if coordinates.ndim != 2 or coordinates.shape[1] != 2:
        raise ValueError("Array coordinates must have shape (n_spots, 2)")
    if not np.isfinite(coordinates).all():
        raise ValueError("Array coordinates must be finite")
    if not np.allclose(coordinates, np.rint(coordinates)):
        raise ValueError("Array coordinates must contain integers")
    coordinates = np.rint(coordinates).astype(np.int64, copy=False)
    coordinate_tuples = [tuple(value) for value in coordinates]
    if len(coordinate_tuples) != len(set(coordinate_tuples)):
        raise ValueError("Array coordinates must be unique within an ROI")

    lookup = {coordinate: index for index, coordinate in enumerate(coordinate_tuples)}
    rows: list[int] = []
    columns: list[int] = []
    for index, (array_row, array_col) in enumerate(coordinate_tuples):
        for delta_row, delta_col in HEX_DELTAS:
            neighbor = lookup.get((array_row + delta_row, array_col + delta_col))
            if neighbor is not None:
                rows.append(index)
                columns.append(neighbor)
    graph = sparse.csr_matrix(
        (np.ones(len(rows), dtype=np.float64), (rows, columns)),
        shape=(len(coordinates), len(coordinates)),
    )
    graph.setdiag(0)
    graph.eliminate_zeros()
    if (graph != graph.T).nnz:
        raise RuntimeError("Native Visium graph is not symmetric")
    return graph


def component_membership(
    graph: sparse.csr_matrix,
    *,
    minimum_spots: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return component id, per-spot size, and retained-component mask."""

    if isinstance(minimum_spots, bool) or int(minimum_spots) < 1:
        raise ValueError("minimum_spots must be a positive integer")
    if graph.shape[0] != graph.shape[1]:
        raise ValueError("Graph must be square")
    if graph.shape[0] == 0:
        empty = np.asarray([], dtype=np.int64)
        return empty, empty, np.asarray([], dtype=bool)
    n_components, component_ids = connected_components(graph, directed=False)
    sizes = np.bincount(component_ids, minlength=n_components).astype(np.int64)
    per_spot_sizes = sizes[component_ids]
    retained = per_spot_sizes >= int(minimum_spots)
    return component_ids.astype(np.int64), per_spot_sizes, retained


def component_center(
    matrix: np.ndarray,
    component_ids: np.ndarray,
) -> np.ndarray:
    """Center each feature within every disconnected graph component.

    Global autocorrelation statistics otherwise treat component-to-component
    mean shifts as spatial structure even though no edge connects the
    components.  Centering removes that unsupported between-component signal
    while leaving within-component gradients unchanged.
    """

    values = np.asarray(matrix, dtype=np.float64)
    was_vector = values.ndim == 1
    if was_vector:
        values = values[:, None]
    if values.ndim != 2:
        raise ValueError("Expression matrix must be one- or two-dimensional")
    membership = np.asarray(component_ids)
    if membership.ndim != 1 or len(membership) != values.shape[0]:
        raise ValueError("Component ids must have one entry per expression row")
    if len(membership) and not np.issubdtype(membership.dtype, np.integer):
        if not np.allclose(membership, np.rint(membership)):
            raise ValueError("Component ids must be integers")
        membership = np.rint(membership).astype(np.int64)

    centered = values.copy()
    for component_id in np.unique(membership):
        mask = membership == component_id
        centered[mask] -= centered[mask].mean(axis=0, keepdims=True)
    return centered[:, 0] if was_vector else centered


def row_normalize_graph(graph: sparse.csr_matrix) -> sparse.csr_matrix:
    """Return row-standardized weights used by both Moran and Geary."""

    graph = graph.astype(np.float64).tocsr(copy=True)
    graph.setdiag(0)
    graph.eliminate_zeros()
    degrees = np.asarray(graph.sum(axis=1)).ravel()
    inverse = np.divide(
        1.0,
        degrees,
        out=np.zeros_like(degrees, dtype=np.float64),
        where=degrees > 0,
    )
    return (sparse.diags(inverse) @ graph).tocsr()


def moran_geary_scores(
    matrix: np.ndarray,
    graph: sparse.csr_matrix,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute global Moran I and Geary C column-wise without smoothing."""

    values = np.asarray(matrix, dtype=np.float64)
    if values.ndim == 1:
        values = values[:, None]
    if values.ndim != 2 or values.shape[0] != graph.shape[0]:
        raise ValueError("Expression matrix rows must equal graph nodes")
    n_spots, n_features = values.shape
    if n_spots < 2:
        missing = np.full(n_features, np.nan, dtype=np.float64)
        return missing.copy(), missing.copy()

    weights = row_normalize_graph(graph)
    s0 = float(weights.sum())
    if s0 <= 0:
        missing = np.full(n_features, np.nan, dtype=np.float64)
        return missing.copy(), missing.copy()

    centered = values - values.mean(axis=0, keepdims=True)
    denominator = np.square(centered).sum(axis=0)
    weighted_centered = weights @ centered
    spatial_cross_product = (centered * weighted_centered).sum(axis=0)
    moran = np.divide(
        (n_spots / s0) * spatial_cross_product,
        denominator,
        out=np.full(n_features, np.nan, dtype=np.float64),
        where=denominator > 0,
    )

    row_sums = np.asarray(weights.sum(axis=1)).ravel()
    column_sums = np.asarray(weights.sum(axis=0)).ravel()
    weighted_values = weights @ values
    squared_difference_sum = (
        (np.square(values) * row_sums[:, None]).sum(axis=0)
        + (np.square(values) * column_sums[:, None]).sum(axis=0)
        - 2.0 * (values * weighted_values).sum(axis=0)
    )
    geary = np.divide(
        ((n_spots - 1.0) / (2.0 * s0)) * squared_difference_sum,
        denominator,
        out=np.full(n_features, np.nan, dtype=np.float64),
        where=denominator > 0,
    )
    return moran, geary


def permutation_pvalues(
    values: np.ndarray,
    graph: sparse.csr_matrix,
    *,
    n_permutations: int,
    seed: int,
) -> tuple[float, float, float, float]:
    """Return observed effects and one-sided empirical permutation p-values."""

    values = np.asarray(values, dtype=np.float64).ravel()
    if isinstance(n_permutations, bool) or int(n_permutations) < 1:
        raise ValueError("n_permutations must be a positive integer")
    observed_moran, observed_geary = moran_geary_scores(values[:, None], graph)
    moran_value = float(observed_moran[0])
    geary_value = float(observed_geary[0])
    if not np.isfinite(moran_value) or not np.isfinite(geary_value):
        return moran_value, geary_value, np.nan, np.nan

    generator = np.random.default_rng(int(seed))
    permutations = np.column_stack(
        [generator.permutation(values) for _ in range(int(n_permutations))]
    )
    null_moran, null_geary = moran_geary_scores(permutations, graph)
    moran_p = float(
        (1 + np.count_nonzero(null_moran >= moran_value))
        / (int(n_permutations) + 1)
    )
    geary_p = float(
        (1 + np.count_nonzero(null_geary <= geary_value))
        / (int(n_permutations) + 1)
    )
    return moran_value, geary_value, moran_p, geary_p


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """BH-adjust finite p-values and preserve missing values."""

    values = np.asarray(pvalues, dtype=np.float64)
    output = np.full(values.shape, np.nan, dtype=np.float64)
    finite = np.isfinite(values)
    if not finite.any():
        return output
    observed = values[finite]
    if ((observed < 0) | (observed > 1)).any():
        raise ValueError("P-values must be between zero and one")
    order = np.argsort(observed, kind="mergesort")
    ranked = observed[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    output[finite] = restored
    return output
