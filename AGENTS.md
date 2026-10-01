# AGENTS.md

## Project Overview
A Nextflow pipeline to process single-cell RNA/ATAC-seq data and pseudobulk it.

There is a helper pipeline `fix_bad_beds.nf` with `fix_bad_beds.config` to repair some
incorrect files on the IGVF Portal. It is not called by `main.nf`. Ignore the helper pipeline
unless it is requested or affected by changes to shared modules, environments, or scripts.

### Environments
The main project environment is managed by pixi. There are subprojects:
- pseudobulk is a Python project managed by uv. It has software tools for performing the
  pseudobulking.
- igvf_portal is a Python project managed by pixi. It has software tools for interacting with the
  IGVF Portal: downloading, uploading, querying, etc.
- visualize_qc is a Python project managed by uv. It has software tools for making summary plots
  for pseudobulk QC.
- environments/CALL_PEAKS.yaml defines a conda environment used in the pipeline.

### Compute resources
The project is intended to be used locally on macOS or Linux, or on a SLURM cluster, specifically
`sherlock`.
- Use existing SSH authentication to access sherlock. If authentication fails, report the issue;
  do not create or modify credentials, keys, or SSH configuration.
- Notable slurm queues:
  - normal: the main submission queue
  - akundaje: a submission queue for the group that I am in
  - owners: a preemptible queue that will usually launch jobs more quickly
  - dev: a queue for testing short-duration low-resource jobs

## Coding Standards
- Follow Python best practices and use type annotations.
- Prefer type narrowing or runtime validation. Use `typing.cast` when the type is already
  established but the checker cannot infer it. Avoid blanket ignores.

## Verification
- Add focused tests for changed behavior when practical, using existing test infrastructure.
  Ask before introducing a new test framework or tests requiring external services or substantial
  resources.
- Run checks for the affected component during development. Run the full suite before completing
  code changes when feasible. Reviews and documentation-only changes do not require the full suite.
- If the full suite is not run, state why and list the checks that were run.
- Some check tasks run `ruff check --fix` and `ruff format`, which modify files. Inspect the
  resulting diff and preserve unrelated user changes.
- Introduce no new failures or warnings. Fix those caused by the change. Report unrelated existing
  failures separately, and do not describe incomplete or failing checks as passing.
- Do not weaken assertions or suppress diagnostics merely to pass checks. Update test expectations
  when the requested behavior changes, while preserving meaningful coverage.

### Nextflow verification
- Lint and preview do not fully verify channel behavior. For channel changes, use small synthetic
  fixtures covering empty input, single-item outputs, partial batches, and shared join keys where
  relevant.
- Isolate test launch directories, work directories, and resume caches from real runs.
- Mock Portal writes in tests; do not run a real upload as verification.

## Commands
The most important commands are all tasks in pixi.toml, which mostly pass arguments to bash scripts in `scripts/`.
This is the preferred design because pixi run [TASK] sets up the environment, but task semantics are much more
limited than a bash script.

| Check | Working directory | Command |
| --- | --- | --- |
| Nextflow syntax, both workflows | Repository root | `pixi run lint-nextflow` |
| Main workflow channel preview, no submitted tasks | Repository root | `pixi run lint-preview` |
| Shell scripts | Repository root | `pixi run lint-scripts` |
| Portal Python project | `igvf_portal/` | `pixi run checks` |
| Pseudobulk Python project, including tests | `pseudobulk/` | `pixi run poe checks` |
| QC visualization Python project | `visualize_qc/` | `pixi run poe checks` |
| Full suite | Repository root | `pixi run checks` |

The Python project checks include automatic fixes and formatting. Use the listed working
directories so the commands select the intended project.

## Boundaries
This section is guidance, not enforcement. Enforced allow/ask/deny rules live in
`.claude/settings.json`; personal overrides go in the untracked `.claude/settings.local.json`.
Keep the two in sync: a boundary added here that is pattern-matchable belongs there too.
- Allowed without prompting:
  - Creation of temporary files as long as they are cleaned up.
  - Read-only bash commands that will not alter or create files, and will not read non-temporary files outside the project folder.
  - Linting, type checking, formatting, or running python tests.
- Agent-submitted Sherlock jobs must request at most five minutes of walltime, two CPUs, and 8 GB
  of total memory. Set explicit limits; do not launch jobs requiring more resources under these
  instructions.
- Do not run resource-intensive jobs on the login node.
- Avoid adding new dependencies to environments.
- Read-only IGVF Portal queries and downloads are permitted. Do not submit files or create, update,
  or delete metadata on any Portal instance, including staging and sandbox.

Keep this root instruction file compact. Put detailed component-specific guidance in nested
`AGENTS.md` files when needed.
