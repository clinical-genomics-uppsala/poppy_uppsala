# Offline package for poppy_uppsala

`build_conda.sh` packages poppy_uppsala, Poppy GMS, a conda environment, the hydra-genetics modules and, optionally, the apptainer images and the reference files, so the pipeline can be run on a cluster without internet access (e.g. Miarka).

The script is run on a machine with internet access (and, for the reference files, access to the marvin paths in the validation file). The result is a few archives that are copied to the cluster without internet access.

## What the script does

1. Creates a conda environment with the pipeline requirements and packs it (`env.tar.gz`).
2. Clones poppy_uppsala and Poppy GMS (the chosen versions), `snakemake-wrappers` and the hydra-genetics modules.
3. Builds the apptainer images, if `-a build` is supplied.
4. Writes offline versions of the configs:
   - container paths: `docker://...` is replaced by `{{APPTAINER_CACHE}}/<image>.sif`
   - reference paths: the marvin paths in the poppy_uppsala config are replaced by `{{DESIGN_REF_FOLDER}}/...`. The pairs of old and new folders are created by `extract_reference_paths.py` from the config and the reference validation file.
5. Downloads the reference and design files listed in the validation file, if `-R true` is supplied.
6. Packs everything into tar files.

## Requirements

- `conda` with `conda-pack` installed in the base environment (`conda install -n base -c conda-forge conda-pack`)
- `git`
- apptainer/singularity, only with `-a build`
- Access to the files in the reference validation file (`file:/` urls on marvin), only with `-R true`

## Build

The script must be run from the directory where the package should be created (`miarka_package` in the example below). The example uses the default values, use `--help` to list all options.

```bash
module load miniconda3
wget https://raw.githubusercontent.com/clinical-genomics-uppsala/poppy_uppsala/<version>/build/build_conda.sh

bash build_conda.sh
```

Build everything, including images and reference files, from a specific version:

```bash
../poppy_uppsala/build/build_conda.sh \
    --poppy-uu-version v2.1.1 --poppy-gms-version v1.1.0 \
    --apptainer-images build --download-references true
```

Options that are not given use the defaults listed in `--help`. If no reference validation file is given, `build/references.hg38.md5sums.yaml` of the cloned poppy_uppsala version is used. The validation file is read in every run, since it is used to create the reference path replacements, also when `-R` is `false`.

## Output of the build

All files are created in the directory the script is run from.

| Output | Created | Content |
|---|---|---|
| `repos-<date>.tar.gz` | always | the pipelines, modules, conda environment and the offline configs, see below |
| `design_and_ref_files-<date>.tar.gz` | `-R true` | reference and design files, stored at the `path:` given for each file in the validation file |
| `apptainer_cache/` | `-a build` | apptainer images (`.sif`), one per container in the configs. This folder is not part of the tar files and must be copied separately |
| `reference_download.log` | `-R true` | output of the reference download |

The folders that are packed (`repos-<date>`, `design_and_ref_files-<date>`) are removed after they have been packed.

Content of `repos-<date>.tar.gz`:

```
repos-<date>/
├── env.tar.gz                          packed conda environment (use this one, see below)
├── venv_poppy_uppsala-<version>/       the conda environment before packing, can be removed
├── poppy_uppsala-<version>/
│   └── config/config_offline.yaml      offline version of config_marvin.yaml
├── poppy-<version>/
│   └── config/config_poppy_offline.yaml  offline version of config.yaml from Poppy GMS
├── snakemake-wrappers/
└── hydra-genetics/<module>/            alignment, annotation, biomarker, cnv_sv, compression, filtering,
                                        fusions, misc, mitochondrial, parabricks, prealignment, qc, reports, snv_indels
```

The two offline configs are used instead of `config.yaml` (Poppy GMS) and `config_marvin.yaml` (poppy_uppsala). They contain the placeholders that must be given when the pipeline is started, see [Run the pipeline](#run-the-pipeline).

### Messages from the build

The script prints a line with the versions and repositories that are packaged, and a line with the reference validation file and python version that are used. After that, the output of `conda`, `git`, `pip`, `hydra-genetics` and `tar -v` follows (one line per packed file).

| Message | Meaning |
|---|---|
| `APPT_IMAGES must be skip or build, got: ...` / `REF_DOWNLOAD must be true or false, got: ...` | invalid option value, nothing has been created |
| `WARNING: <config> <key> <path> is inside a list, reference-path-update will not replace it` | printed by `extract_reference_paths.py`. `hydra-genetics` does not replace paths inside yaml lists, so the path is replaced by the `sed` command in the script instead. If the warning is about a path other than `merge_cnv_json/annotations`, that path will still point to marvin in the offline config |
| `Cannot create the path pairs:` followed by a list | `extract_reference_paths.py` found an absolute path in the config that is not in the validation file, or a file that has a different name in the config and the validation file. The script stops. Add the file to the validation file, or fix the config |
| `Downloading references.` | the reference download starts (`-R true`) |
| `UPDATED: <n>` / `NOT UPDATED: <n>` | printed by `hydra-genetics references download`: files downloaded, and files that were skipped because the checksum is unchanged |
| `Reference download had failures, see reference_download.log` | one or more files could not be downloaded. The script stops before packing. `hydra-genetics` itself reports failures only in the log |

`hydra-genetics prepare-environment create-singularity-files` does not stop if an image cannot be built. When building images (`-a build`), check the output for `singularity build` errors and that `apptainer_cache/` contains one `.sif` per container.

## Install on the cluster

Copy the archives (and `apptainer_cache/`) to the cluster and extract them at the desired location:

```bash
tar -xzvf repos-<date>.tar.gz
tar -xzvf design_and_ref_files-<date>.tar.gz
```

Unpack the conda environment into a new folder and fix its paths (conda-pack):

```bash
mkdir repos-<date>/venv
tar -xzf repos-<date>/env.tar.gz -C repos-<date>/venv
source repos-<date>/venv/bin/activate
conda-unpack
```

The apptainer images are expected in the `apptainer_cache/` folder that is given as `APPTAINER_CACHE` below. The miarka profile (`profile/miarka_drmaa/config.yaml`) points to a fixed apptainer cache (`singularity-prefix`) and to `snakemake-wrappers` in the repo folder (`wrapper-prefix`, set by the build script). Check that these paths match where the files are placed.

## Run the pipeline

Activate the environment and start snakemake from the run directory. Variables used in the command:

| Variable | Content | Notes |
|---|---|---|
| `poppy_gms_path` | `repos-<date>/poppy-<version>` |
| `poppy_uu_path` | `repos-<date>/poppy_uppsala-<version>` |
| `reference_directory` | the run-specific reference files, `${design_ref_path}/NovaSeq` or `${design_ref_path}/NextSeq`, replaces `{{REFERENCE_DIRECTORY}}` |
| `poppy_version`, `poppy_uu_version` | versions of Poppy GMS and poppy_uppsala |
| `sequenceid` | id of the sequencing run |
| `hydra_modules_path` | `repos-<date>/hydra-genetics`, the cloned modules. Given as `hydra_local_path`, otherwise the modules are fetched from GitHub | Specific for offline mode |
| `apptainer_cache` | folder with the `.sif` images |  Specific for offline mode |
| `design_ref_path` | extracted `design_and_ref_files-<date>`, replaces `{{DESIGN_REF_FOLDER}}` | Specific for offline mode |

```bash
source ${repos_path}/venv/bin/activate

snakemake --profile ${poppy_uu_path}/profile/miarka_drmaa -s ${poppy_uu_path}/workflow/Snakefile \
        --configfiles ${poppy_gms_path}/config/config_poppy_offline.yaml ${poppy_uu_path}/config/config_offline.yaml \
        --config poppy_version=${poppy_version} \
            poppy_uu_version=${poppy_uu_version} \
            POPPY_HOME=${poppy_gms_path} \
            POPPY_UU_HOME=${poppy_uu_path} \
            REFERENCE_DIRECTORY=${reference_directory} \
            sequenceid=${sequenceid} \
            hydra_local_path=${hydra_modules_path} \
            APPTAINER_CACHE=${apptainer_cache} \
            DESIGN_REF_FOLDER=${design_ref_path}
```

Notes:

- Without internet access, the pipeline starts with a short wait while it tries to reach GitHub, and then uses `POPPY_HOME` for Poppy GMS.
- The machine-specific files are stored in `NovaSeq/reference_files/` or `NextSeq/reference_files/` in the design and reference folder. `REFERENCE_DIRECTORY` should therefore point to `NovaSeq` or `NextSeq` in that folder.
- Any other config files used for a run (e.g. samples and units) are added to `--configfiles` as usual.

## Known limitations

- The hydra-genetics modules and `snakemake-wrappers` are cloned from their default branches. They are not pinned to the versions used by the pipeline.
- The apptainer images are not included in the tar files, copy `apptainer_cache/` separately. With `-a skip` no images are built, but the offline configs still point to `{{APPTAINER_CACHE}}`.
- The profile and the example paths are specific to Miarka.
