# Reporting and resource monitoring

`run_with_resource_monitor.py` is a side-band wrapper for long standalone
modules. It records a process-tree time series, command log, and JSON summary
containing wall time, whole-machine-normalized CPU, RSS/VMS, observed I/O,
project size, and filesystem free space. Thresholds warn but do not terminate
scientific work; thread/process limits still belong to the wrapped command.
The JSON summary converts paths inside the project to relative paths and
redacts other absolute paths. Child-process text written to the command log
cannot be sanitized automatically, so review that log before sharing it.

Example:

```bash
python \
  workflow/scripts/reporting/run_with_resource_monitor.py \
  --cwd . \
  --project-root . \
  --series-output logs/resources/my_step.series.tsv \
  --summary-output logs/resources/my_step.summary.json \
  --command-log logs/resources/my_step.command.log \
  --cpu-warn-percent 40 \
  --project-warn-gib 10 \
  --project-critical-gib 20 \
  -- env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  python path/to/module.py
```

CPU is reported as a percentage of total logical-machine capacity, not as the
usual per-core `ps` percentage. Concurrent wrappers each describe only their
own descendant tree, so cohort-wide concurrency must also be controlled by the
caller.

After a run, `summarize_resource_logs.py` consolidates monitor JSON files into
one TSV plus a small JSON overview. It preserves the distinction between
per-process-tree peaks, sampled project size, and non-additive concurrent wall
times.

```bash
python \
  workflow/scripts/reporting/summarize_resource_logs.py \
  --resource-dir logs/resources \
  --resource-summary results/module_a/resource_summary.json \
  --resource-summary results/module_b/resource_summary.json \
  --table-output results/reporting/resource_run_summary.tsv \
  --summary-output results/reporting/resource_run_summary.json
```

Discovery accepts both `*.summary.json` and `*.resources.json` payloads after
checking their schema; unrelated science JSON files are ignored. Generic names
such as `resource_summary.json` receive parent-directory context in `step_id`.
