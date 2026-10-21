## Summary

<!-- What does this PR change, and why? Link the related issue(s). -->

## Verification evidence

<!-- Paste the commands you ran and their results. Include the benchmark manifest if results changed. -->

- Commands run:
- Seeds / runs reported (all of them, no cherry-picking):
- Benchmark manifest (if results regenerated):

## Checklist

- [ ] Tests pass locally (`pytest`).
- [ ] `ruff check .` passes.
- [ ] `ruff format --check .` passes.
- [ ] New acceptance criteria were pre-registered in `docs/benchmarks/` before results were produced (link the commit): <!-- commit hash -->
- [ ] Every seed and run used is reported, including failed or discarded ones.
- [ ] Every new equation or parameter is cited with a DOI at the point of use.
- [ ] Each claim is labelled: reproduction / independent implementation / new contribution / unverified hypothesis.
- [ ] No GPL or AGPL code was copied into this change.
- [ ] Results were regenerated with `benchmarks/*/run.py` and the manifest is included.
