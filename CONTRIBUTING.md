# Contributing

This project values a small, verified feature set over a large, unverified one.

## Rules for scientific changes

1. **Pre-register before you measure.** A new benchmark or acceptance
   criterion is committed as a document in `docs/benchmarks/` *before* the
   code or results that it evaluates. A criterion is changed only through a
   dated amendment section explaining why. The original text is never
   rewritten.
2. **Report every run.** Benchmarks report all seeds and configurations that
   were run. Cherry-picking is not allowed.
3. **Independent references.** A new method is compared against an
   independent implementation (exact solution, reference solver, or an
   external open-source code). The method under test and its reference share
   only the problem specification.
4. **Cite at the point of use.** Equations, models and parameters carry a
   citation (with DOI) in the docstring or document where they appear.
   Material parameters require a primary source that has actually been read.
5. **Label claims.** State whether something is a reproduction of a published
   result, an independent implementation of a published method, a new
   contribution, or an unverified hypothesis.
6. **No copied copyleft code.** GPL/AGPL projects may be read for method
   understanding, but their code is not copied here.

## Workflow

- Small, focused pull requests with conventional-commit messages
  (`feat(fv): ...`, `docs(benchmark): ...`, `fix(pinn): ...`).
- `pytest` must pass locally. CI runs the same tests with and without torch.
- Open an issue for every known scientific or engineering limitation.

## Development setup

```bash
python -m venv .venv
.venv/bin/pip install torch --index-url https://download.pytorch.org/whl/cpu
.venv/bin/pip install -e ".[test]"
.venv/bin/pytest
```
