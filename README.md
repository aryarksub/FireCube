# PyroStack

Code used to generate the PyroStack wildfire trajectory dataset.

## Setup

Install the core datacube-generation environment with uv:

```bash
uv sync
```

Install the optional dependencies for the machine learning training code:

```bash
uv sync --extra ml
```

Run commands through uv so they use the project environment:

```bash
uv run python -c "import main"
```
