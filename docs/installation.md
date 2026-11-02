# Installation

`annplyr` supports Python 3.12 and newer.

## Install from PyPI

Choose the package manager already used by your project:

```bash
pip install annplyr
```

or:

```bash
uv add annplyr
```

`annplyr` has no plotting framework or Scanpy dependency. It works with the
AnnData object already produced by your scverse workflow.

## Wrap your AnnData

Wrap an object once with `ap.tidy`, then call verbs on it:

```python
import anndata as ad
import annplyr as ap
from annplyr.data import obs

adata = ap.tidy(ad.read_h5ad("analysis.h5ad"))
filtered = adata.filter(obs.qc_pass)
```

`ap.tidy` does not copy or read the matrices. Importing `annplyr` also registers
the `adata.ap` accessor, so `adata.ap.filter(...)` works on a plain AnnData.


## Check the installation

```python
import annplyr

print(annplyr.__version__)
```

Continue with the {doc}`quickstart` for a complete, small cohort workflow.

## Install for development

Contributors can create the complete local environment from a clone:

```bash
git clone https://github.com/mdmanurung/annplyr
cd annplyr
uv sync --all-extras
```

The contributor workflow and verification commands are documented in
{doc}`contributing`.
