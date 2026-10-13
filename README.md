# <img src="images/hydragenetics.png" width=40 /> hydra-genetics/poppy_uppsala

#### Uppsala config and final output changes for Twist Myeloid GMS Poppy pipeline

![Lint](https://github.com/hydra-genetics/poppy_uppsala/actions/workflows/lint.yaml/badge.svg?branch=develop)
![Snakefmt](https://github.com/hydra-genetics/poppy_uppsala/actions/workflows/snakefmt.yaml/badge.svg?branch=develop)
![snakemake dry run](https://github.com/hydra-genetics/poppy_uppsala/actions/workflows/snakemake-dry-run.yaml/badge.svg?branch=develop)
![integration test](https://github.com/hydra-genetics/poppy_uppsala/actions/workflows/integration.yaml/badge.svg?branch=develop)

![pycodestyle](https://github.com/hydra-genetics/poppy_uppsala/actions/workflows/pycodestyle.yaml/badge.svg?branch=develop)
![pytest](https://github.com/hydra-genetics/poppy_uppsala/actions/workflows/pytest.yaml/badge.svg?branch=develop)

[![License: GPL-3](https://img.shields.io/badge/License-GPL3-yellow.svg)](https://opensource.org/licenses/gpl-3.0.html)

## :speech_balloon: Introduction

The module consists of alignment  ....

## :heavy_exclamation_mark: Dependencies

In order to use this module, the following dependencies are required:

[![hydra-genetics](https://img.shields.io/badge/hydragenetics-v1.3.0-blue)](https://github.com/hydra-genetics/)
[![pandas](https://img.shields.io/badge/pandas-1.3.1-blue)](https://pandas.pydata.org/)
[![python](https://img.shields.io/badge/python-3.8-blue)
[![snakemake](https://img.shields.io/badge/snakemake-6.8.0-blue)](https://snakemake.readthedocs.io/en/stable/)
[![singularity](https://img.shields.io/badge/singularity-3.0.0-blue)](https://sylabs.io/docs/)

## :school_satchel: Preparations

### Sample data

Input data should be added to [`samples.tsv`](https://github.com/hydra-genetics/poppy_uppsala/blob/develop/config/samples.tsv)
and [`units.tsv`](https://github.com/hydra-genetics/poppy_uppsala/blob/develop/config/units.tsv).
The following information need to be added to these files:

| Column Id | Description |
| --- | --- |
| **`samples.tsv`** |
| sample | unique sample/patient id, one per row |
| **`units.tsv`** |
| sample | same sample/patient id as in `samples.tsv` |
| type | data type identifier (one letter), can be one of **T**umor, **N**ormal, **R**NA |
| platform | type of sequencing platform, e.g. `NovaSeq` |
| machine | specific machine id, e.g. NovaSeq instruments have `@Axxxxx` |
| flowcell | identifer of flowcell used |
| lane | flowcell lane number |
| barcode | sequence library barcode/index, connect forward and reverse indices by `+`, e.g. `ATGC+ATGC` |
| fastq1/2 | absolute path to forward and reverse reads |
| adapter | adapter sequences to be trimmed, separated by comma |

## :white_check_mark: Testing

The checks below are the same as the ones run by the GitHub Actions. Run them from the pipeline environment (Python 3.9) with the test requirements installed:

```bash
$ pip install -r requirements.txt -r requirements.test.txt
```

### Snakemake dry run and lint

The workflow repository contains a small test dataset `.tests/integration`:

```bash
$ cd .tests/integration
$ snakemake -n -s ../../workflow/Snakefile \
    --configfiles config/config_poppy_v1.1.0.yaml ../../config/config_marvin.yaml config/config_test.yaml \
    --config POPPY_HOME=../.. POPPY_UU_HOME=../.. REFERENCE_DIRECTORY=. sequenceid="dry_run" poppy_version="v1.1.0" poppy_uu_version="v2.0.0"
$ snakemake --lint -n -s ../../workflow/Snakefile \
    --configfiles config/config_poppy_v1.1.0.yaml ../../config/config_marvin.yaml config/config_test.yaml \
    --config POPPY_HOME=../.. POPPY_UU_HOME=../.. REFERENCE_DIRECTORY=. sequenceid="lint" poppy_version="v1.1.0" poppy_uu_version="v2.0.0"
```
`config/config_poppy_v1.1.0.yaml` is the poppy GMS config, `../../config/config_marvin.yaml` the Uppsala config and `config/config_test.yaml` the test config. Config files given later overwrite overlapping variables in earlier ones.

### Python scripts

```bash
$ pytest workflow/scripts
$ pycodestyle --max-line-length=130 --statistics workflow/scripts
```

### Snakefmt

The GitHub Action runs snakefmt through [super-linter](https://github.com/super-linter/super-linter) (v8.7.0, which uses snakefmt 2.0.2). snakefmt 2 requires Python >= 3.11 and can therefore not be installed in the pipeline environment. Instead, run it with [uv](https://docs.astral.sh/uv/) (install once, e.g. `pip install uv`), which creates a small cached environment with the pinned version:

```bash
$ uvx --from snakefmt==2.0.2 --with shfmt-py snakefmt --config .github/linters/.snakefmt.toml --check workflow/
```
Remove `--check` to let snakefmt format the files. If super-linter is updated in `.github/workflows/snakefmt.yaml`, update the snakefmt version above to the one listed in super-linter's `dependencies/python/snakefmt.txt`.
## :rocket: Usage

To use this module in your workflow, follow the description in the
[snakemake docs](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules).
Add the module to your `Snakefile` like so:

```bash
module prealignment:
    snakefile:
        github(
            "poppy_uppsala",
            path="workflow/Snakefile",
            tag="1.0.0",
        )
    config:
        config


use rule * from poppy_uppsala as poppy_uppsala_*
```

### Output files

The following output files should be targeted via another rule:

| File | Description |
|---|---|
| `poppy_uppsala/PATH/FILE` | DESCRIPTION |

## :judge: Rule Graph
![rule_graph_reference](images/rulegraph.svg)