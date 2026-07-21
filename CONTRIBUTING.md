# Contributing to PackMem2

## Reporting bugs, asking for features

If you encounter a bug or would like to suggest a new feature, please create a [issue](https://github.com/patrickfuchs/PackMem2/issues). Please be very specific about the bug you encountered or your requirements. Provide as many examples as possible.

Do not create a pull request before being invited to do so by the development team

## Developping

Fork and/or clone the GitHub repo:

```bash
git clone git@github.com:patrickfuchs/PackMem2.git
```

We use [uv](https://docs.astral.sh/uv/) to manage all dependencies.

Always valide the current tests upon adding new features or fixing bugs:

```bash
# Fast tests
uv run pytest
# All tests, slow (2-3 minutes)
uv run pytest -m "not slowest"
```


