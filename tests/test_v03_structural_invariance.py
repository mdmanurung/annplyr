from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

import annplyr as ap
from annplyr._grouped import GroupedAnnData


@pytest.mark.filterwarnings("ignore:Observation names are not unique")
@pytest.mark.filterwarnings("ignore:Variable names are not unique")
def test_chained_grouped_fixture_preserves_structure_ownership_and_positional_order() -> None:
    obs = pd.DataFrame(
        {
            "group": pd.Categorical(["b", "a", "b", "a"], categories=["a", "b", "unused"]),
            "value": pd.array([1, 4, 3, 2], dtype="Int64"),
            "flag": pd.array([True, pd.NA, False, True], dtype="boolean"),
            "label": pd.array(["w", "x", "y", "z"], dtype="string"),
        },
        index=["dup", "dup", "dup", "dup"],
    )
    var = pd.DataFrame({"kind": pd.Categorical(["rna", "protein", "rna"])}, index=["gene", "gene", "gene"])
    values = np.arange(12, dtype=np.float32).reshape(4, 3)
    source = AnnData(X=sparse.csr_matrix(values), obs=obs, var=var)
    source.layers["counts"] = sparse.csr_matrix(values + 100)
    source.raw = source.copy()
    source.obsm["embedding"] = np.arange(8, dtype=np.float32).reshape(4, 2)
    source.varm["loading"] = np.arange(6, dtype=np.float32).reshape(3, 2)
    source.obsp["connectivity"] = sparse.csr_matrix(np.arange(16).reshape(4, 4))
    source.varp["similarity"] = sparse.csr_matrix(np.arange(9).reshape(3, 3))

    expected_positions = np.array([2, 0, 1, 3])
    expected_x = source.X[expected_positions, :].toarray()
    expected_layer = source.layers["counts"][expected_positions, :].toarray()
    expected_raw = source.raw.X[expected_positions, :].toarray()
    expected_obsm = source.obsm["embedding"][expected_positions, :].copy()
    expected_obsp = source.obsp["connectivity"][expected_positions, :][:, expected_positions].toarray()

    grouped = source.ap.group_by(obs="group")
    result = grouped.mutate(obs={"within_group": ap.row_number(), "double": ap.col("value") * 2}).arrange(
        obs=ap.desc("value")
    )

    assert isinstance(result, GroupedAnnData)
    assert result is not grouped
    assert result.group_vars() == ["group"]
    assert result.group_keys()["group"].tolist() == ["b", "a"]
    output = result.ungroup()
    assert output is not source
    assert output.obs_names.tolist() == ["dup"] * 4
    assert output.var_names.tolist() == ["gene"] * 3
    assert output.obs["value"].tolist() == [3, 1, 4, 2]
    assert output.obs["within_group"].tolist() == [2, 1, 1, 2]
    assert output.obs["double"].tolist() == [6, 2, 8, 4]
    assert isinstance(output.obs["group"].dtype, pd.CategoricalDtype)
    assert str(output.obs["value"].dtype) == "Int64"
    assert str(output.obs["flag"].dtype) == "boolean"
    assert str(output.obs["label"].dtype).startswith("string")

    assert sparse.isspmatrix_csr(output.X)
    assert sparse.isspmatrix_csr(output.layers["counts"])
    assert output.raw is not None and sparse.isspmatrix_csr(output.raw.X)
    assert np.array_equal(output.X.toarray(), expected_x)
    assert np.array_equal(output.layers["counts"].toarray(), expected_layer)
    assert np.array_equal(output.raw.X.toarray(), expected_raw)
    assert np.array_equal(output.obsm["embedding"], expected_obsm)
    assert np.array_equal(output.obsp["connectivity"].toarray(), expected_obsp)
    assert output.obsm["embedding"].shape == (4, 2)
    assert output.varm["loading"].shape == (3, 2)
    assert output.obsp["connectivity"].shape == (4, 4)
    assert output.varp["similarity"].shape == (3, 3)

    assert not np.shares_memory(output.X.data, source.X.data)
    assert not np.shares_memory(output.layers["counts"].data, source.layers["counts"].data)
    assert not np.shares_memory(output.raw.X.data, source.raw.X.data)
    assert not np.shares_memory(output.obsm["embedding"], source.obsm["embedding"])
    assert not np.shares_memory(output.varm["loading"], source.varm["loading"])
    assert not np.shares_memory(output.obsp["connectivity"].data, source.obsp["connectivity"].data)
    assert not np.shares_memory(output.varp["similarity"].data, source.varp["similarity"].data)

    output.obs.iloc[0, output.obs.columns.get_loc("value")] = 99
    output.X.data[0] = -1
    assert source.obs["value"].tolist() == [1, 4, 3, 2]
    assert source.X.toarray().tolist() == values.tolist()
