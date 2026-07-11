# Task Breakdown 12 — #11 ODMD-Anywhere (signal-processing spinout)
Goal: standalone library + the head-to-head benchmark notebook that justifies its existence.

1. **Extraction refactor** — factor ODMD core out of the solver into `odmd/` with zero quantum/chemistry imports: `estimate(times, signal, n_modes) -> [(freq, decay, amplitude, error_bar)]`. Solver becomes client #1.
   ✓ Solver test suite green against the extracted package (proves the cut is clean). (L)
2. **Benchmark harness** — synthetic signals: known mode mixtures × noise levels × record lengths; competitors: classical DMD, ESPRIT, Prony (existing Python impls), matrix pencil.
   ✓ Harness runs all methods on identical data; metrics: frequency error, success rate vs SNR, error-bar calibration (ODMD's differentiator — do competitors even provide one?). (M)
3. **The verdict notebook (publish BEFORE the package)** — where does ODMD win, where does it lose. If it doesn't win anywhere meaningful, the spinout stops here — cheap kill, per the mini-spec's own risk note.
   ✓ Honest regime map; go/no-go decision written at the bottom. (M)
4. **(If go) One real-signal demo** — a classical MD velocity-autocorrelation function or public NMR FID → extracted spectrum vs reference.
   ✓ Real-data example in README. (M)
5. **(If go) Package + release** — API docs, PyPI, the benchmark notebook as the front page.
   ✓ `pip install` cold-start to first estimate in < 5 lines. (S)
