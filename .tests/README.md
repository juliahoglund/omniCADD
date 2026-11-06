# Tests

This directory contains the test configurations and test data used to check that the workflow's job graph builds
under two data scenarios. The tests are **dry runs**: Snakemake builds and checks the DAG, it does not execute the
scripts or tools, so they catch rule, config and wiring errors but not runtime problems such as memory use or bugs
inside scripts.

## Contents

```
.tests/
├── config/
│   ├── config_standard.yaml   # Standard data tier test configuration
│   └── config_limited.yaml    # Limited data tier test configuration
├── resources/                 # Minimal chr18 test data (see resources/README.md)
├── annotation/, genome/, grantham/, snpEff/   # Small support files (used by config_production.yaml and snpeff_annotation.smk, not by the two test configs)
├── test_configs.sh            # Runs the dry run for both configs
└── test_debug.sh              # Debug helper for bash behaviour in CI (set -e with arithmetic)
```

## The two scenarios

Both use `sus_scrofa` and **chromosome 18 only**. Where they differ:

| Setting | Standard | Limited |
|---------|----------|---------|
| `data_tier` | `standard` | `limited` |
| Gene annotation | `gff` (existing annotation), no gene prediction | `augustus` (gene prediction on) |
| Variant annotator | VEP | SnpEff |
| Alignment entry | `44_mammals.epo`, full species list | `44_mammals.epo`, a short list of six species |
| `build_alignment` (minimap2) | not set | enabled |
| Preprocessing | none needed (`needs_*` false) | formatting, filtering and VCF processing on; chromosome-name standardization on |
| GERP tree | `resources/tree_44_mammals.nwk` | `resources/tree_custom.nwk` |
| phastCons / phyloP precision | `HIGH` | `MEDIUM` |
| Cross-validation | 10 folds, C in [0.1, 1.0, 10.0] | 5 folds, C = 1.0 |
| Final model C | 10.0 | 1.0 |

Notes:

- The `data_tier` value matters: `standard` and `limited` apply a preset from the `profiles:` section that
  **overrides** explicit settings of the same keys (for example the variant annotator and the ancestral-sequence
  source). See `docs/configuration.rst` and `docs/DATA_TIERS.md`.
- The limited scenario exercises the alignment-free mode, which is still under development.

## Running the tests

### One scenario, by hand

```bash
snakemake --dry-run test_dag --configfile .tests/config/config_standard.yaml
snakemake --dry-run test_dag --configfile .tests/config/config_limited.yaml
```

`test_dag` is a CI target defined in the `snakefile` that requests the main outputs of all eight modules.

### Both scenarios, with the script

```bash
bash .tests/test_configs.sh
```

For each scenario the script

1. **copies the test config over `config/config.yaml`** (and creates `config/config_refs.yaml` and
   `config/resources.yaml` from the `*.default.yaml` templates if they are missing),
2. runs `snakemake --dry-run test_dag --configfile config/config.yaml`,
3. prints PASSED or FAILED and writes the output to `/tmp/snakemake_<scenario>.log`.

> **Warning**: because `config/config.yaml` is overwritten, it is left holding the last scenario that ran. It is a
> tracked file, so restore it afterwards with `git checkout config/config.yaml`. Do not keep your own run configuration
> under that name.

## Continuous integration

`.github/workflows/snakemake-dry-run.yml` runs on pushes and pull requests to `main`, `master`, `dev` and
`restructure` that change a `.smk` file, the `snakefile`, `.tests/config/**` or `workflow/scripts/**`. It

1. checks that every Python script in `workflow/scripts` compiles (files of two lines or fewer are treated as placeholders),
2. copies the default `config_refs` and `resources` templates and runs `test_configs.sh`,
3. runs `snakemake --lint` with the standard test config and fails if more than two lint warnings are reported.

Other workflows check formatting and style: `snakefmt.yml`, `snakemake-lint.yml` and `pycodestyle.yml`.

## Adding a test scenario

1. Copy `config_standard.yaml` or `config_limited.yaml` and change the relevant sections (for example `data_tier`,
   `annotation`, `ancestral_sequence`). Keep all six entries in `containers:`.
2. Add a `test_config` call for it in `test_configs.sh`.
3. Document it in the table above.

## Test data

The two test configs read their data from `.tests/resources/`. The minimal chr18 data (a 600 bp synthetic
reference, a handful of variants, a five-block alignment, a tiny annotation) is described in [resources/README.md](resources/README.md). It is only meant to let the configs resolve their input
files; it is too small to produce meaningful scores.
