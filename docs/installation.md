# Installation

Install `annplyr` from PyPI:

```bash
pip install annplyr
```

For development from the GitHub repository:

```bash
git clone https://github.com/mdmanurung/annplyr
cd annplyr
uv sync --all-extras
```

Importing the package registers the AnnData accessor:

```python
import annplyr as ap

adata.ap
```

## Development Checks

The common local checks are:

```bash
pytest -q
ruff check src tests
ruff format src tests --check
python -m mypy src/annplyr
uvx hatch run docs:build
```

Use `UV_CACHE_DIR=/tmp/uv-cache` in restricted environments where the default
uv cache location is read-only.
