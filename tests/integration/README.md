# Representative downstream fixture

The integration fixture is a wholly synthetic 8-cell by 6-feature AnnData
object. Its explicit integer counts and metadata were written for annplyr; they
do not derive from a person, organism sample, published dataset, or external
database. The fixture is distributed under the repository's BSD-3-Clause
license.

`representative_fixture.py` is the authoritative source. It deliberately
contains categorical columns, nullable integer and Boolean columns, duplicate
observation names, CSR `X` and layer data, `raw`, dense `obsm`/`varm`, CSR
`obsp`, and CSC `varp`. Explicit values make regeneration deterministic across
Python and NumPy releases.

No generated H5AD or Zarr dataset is committed. To materialize both forms in a
disposable location, run:

```console
uv run --group integration python tests/integration/representative_fixture.py /tmp/annplyr-fixture
```

The integration tests regenerate into pytest's temporary directory and compare
values, pandas dtypes and categories, sparse formats, axis order, aligned
container shapes, independent ownership, positional identity, and grouping
state. They also reconstruct the fixture with `anndata.concat` and pass an
annplyr result through Scanpy normalization, log transformation, and PCA.
