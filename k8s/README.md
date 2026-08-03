# Kubernetes deployment

A benchmark Job plus a weekly cross-check CronJob, over a shared data PVC.

```bash
docker build -t chem:0.1.0 .
kubectl apply -k k8s/overlays/dev     # namespace chem-dev  (CronJob suspended)
kubectl apply -k k8s/overlays/prod    # namespace chem-prod (CronJob active)
```

| Workload | Command | Cadence |
|---|---|---|
| `chem-benchmark-krylov` | `python benchmark_krylov.py` | on apply |
| `chem-cross-check` | `python cross_check.py` | Sun 03:00 |

`cross_check.py` is the scheduled canary because it requires four independent
methods to agree — a silent physics regression surfaces as disagreement rather
than a crash, which matches this repo's stated "falsifiable honesty" culture.

## What is deliberately not used

**`run_sprint_benchmarks.sh` is not wired up.** It imports
`orchestrate_hybrid_pipeline`, which `CLAUDE.md` quarantines as the old broken
core retained only as a regression fixture. It also writes `temp_runner.py` into
the working directory, which `readOnlyRootFilesystem` forbids outright. The
Jobs use the validated benchmarks named in CLAUDE.md instead.

DMRG/block2 is not installed in the default image either — block2 initialises
its own OpenMP runtime and segfaults if it loads into a process that already
imported pyscf or qiskit-aer. It belongs in a separate image, not this one.

## Threads must equal the CPU limit

This is the one setting most likely to bite. PySCF and the BLAS underneath size
their thread pools from `OMP_NUM_THREADS`, **not** from the cgroup. Set it above
the pod's CPU limit and the kernel throttles a process that thinks it has more
cores than it can use — slower than running single-threaded.

That is why the thread counts live in the `chem-env` ConfigMap and each overlay
sets them alongside its CPU value:

| Overlay | cpu | `OMP_NUM_THREADS` |
|---|---|---|
| base | 4 | 4 |
| dev | 1 | 1 |
| prod | 8 | 8 |

**Change one, change the other.** The image itself defaults to 1 so a plain
`docker run` stays predictable instead of grabbing every core on the host.

## Volumes

- `/app/data` → PVC. All benchmarks write relative `data/…` paths (e.g.
  `data/krylov_benchmark.csv`). Sized 20Gi/100Gi mostly for SCF checkpoints:
  `benchmark_nbn.py` caches a ground state to `data/nbn_scf.chk` to make the
  spin scan resumable, and those reach GBs.
- `/app/.dmrg_tmp` → emptyDir. Large, purely intermediate, worthless after the run.
- `/tmp` → emptyDir. `HOME` and `MPLCONFIGDIR` both point here; `/home/app` is
  on the read-only root filesystem, so anything writing a `$HOME` dotfile
  (PySCF's `~/.pyscf_conf.py`, matplotlib's cache) fails without it.

## Dev CronJob is suspended

A weekly multi-hour regression check has no value in a dev namespace and would
quietly burn a CPU-week a year. Trigger it manually instead:

```bash
kubectl create job --from=cronjob/chem-cross-check-dev manual-run -n chem-dev
```
