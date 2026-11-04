#!/usr/bin/env bash
# Package poppy_uppsala (and everything it needs) so it can be copied to a cluster without internet access.
#
# Steps:
#   1. create a conda env with the pipeline requirements and pack it
#   2. clone poppy_uppsala, Poppy GMS, snakemake-wrappers and the hydra-genetics modules
#   3. (optional) build the apptainer images
#   4. write offline versions of the configs (container and reference paths replaced)
#   5. (optional) download the reference and design files
#   6. tar everything
#
# Run with --help for the options.

set -eo pipefail


# ------------------------------------------------------------------------------
# Help text
# ------------------------------------------------------------------------------
usage() {
    cat <<USAGE
Usage: $(basename "$0") [options]

Package the poppy_uppsala pipeline with its conda env, Poppy GMS, hydra-genetics modules,
snakemake-wrappers and optionally apptainer images and reference files so it can be copied to the cluster.
Run this from the directory where the package should be created.

Outputs:
  <repo-folder>.tar.gz               pipelines, conda environment and modules
  <design-ref-folder>.tar.gz         reference and design files (only with -R true)
  apptainer_cache/                   apptainer images (.sif), only with -a build

Options:
  -n, --repo-folder NAME             folder for the pipelines, modules and conda environment
                                     (default: repos-YYYYMMDD)
  -d, --design-ref-folder NAME       folder for the reference and design files, used with -R true
                                     (default: design_and_ref_files-YYYYMMDD)
  --poppy-uu-repo URL                git repository of poppy_uppsala
                                     (default: https://github.com/clinical-genomics-uppsala/poppy_uppsala.git)
  --poppy-uu-version VERSION         poppy_uppsala tag or branch to package
                                     (default: v2.2.0)
  --poppy-gms-repo URL               git repository of Poppy GMS
                                     (default: https://github.com/genomic-medicine-sweden/poppy.git)
  --poppy-gms-version VERSION        Poppy GMS tag or branch to package
                                     (default: v1.1.0)
  --python-version VERSION           python version of the conda environment
                                     (default: 3.9)
  -R, --download-references true|false
                                     download the reference files listed in the validation file
                                     (default: false)
  -r, --reference-validation-file FILE
                                     validation file (.yaml or .json) listing the reference files to download,
                                     only used with -R true
                                     (default: build/references.hg38.md5sums.yaml of the poppy_uppsala version)
  -a, --apptainer-images build|skip
                                     build: download the apptainer images into apptainer_cache/
                                     skip: do not download; the configs still point to apptainer_cache/,
                                     so the images must already be there
                                     (default: skip)
  -h, --help                         show this help
USAGE
}


# ------------------------------------------------------------------------------
# Options and defaults
# ------------------------------------------------------------------------------
# Read the options. Anything not given is set in the defaults block below.
while [[ $# -gt 0 ]]; do
    case "$1" in
        -n|--repo-folder) REPO_FOLDER="$2"; shift 2 ;;
        -d|--design-ref-folder) DESIGN_REF_FOLDER="$2"; shift 2 ;;
        --poppy-uu-repo) POPPY_UU_REPO="$2"; shift 2 ;;
        -u|--poppy-uu-version) POPPY_UU_VERSION="$2"; shift 2 ;;
        --poppy-gms-repo) POPPY_GMS_REPO="$2"; shift 2 ;;
        -g|--poppy-gms-version) POPPY_GMS_VERSION="$2"; shift 2 ;;
        -p|--python-version) PYTHON_VERSION="$2"; shift 2 ;;
        -R|--download-references) REF_DOWNLOAD="$2"; shift 2 ;;
        -r|--reference-validation-file) REFERENCE_VALIDATION_FILE="$2"; shift 2 ;;
        -a|--apptainer-images) APPT_IMAGES="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

# All default values are collected here. They are set after the options are read
# so that a default can use the final value of another option.
REPO_FOLDER="${REPO_FOLDER:-repos-$(date +%Y%m%d)}"
DESIGN_REF_FOLDER="${DESIGN_REF_FOLDER:-design_and_ref_files-$(date +%Y%m%d)}"
POPPY_UU_REPO="${POPPY_UU_REPO:-https://github.com/clinical-genomics-uppsala/poppy_uppsala.git}"
POPPY_UU_VERSION="${POPPY_UU_VERSION:-v2.2.0}"
POPPY_GMS_REPO="${POPPY_GMS_REPO:-https://github.com/genomic-medicine-sweden/poppy.git}"
POPPY_GMS_VERSION="${POPPY_GMS_VERSION:-v1.1.0}"
PYTHON_VERSION="${PYTHON_VERSION:-3.9}"
REF_DOWNLOAD="${REF_DOWNLOAD:-false}"
REFERENCE_VALIDATION_FILE="${REFERENCE_VALIDATION_FILE:-./${REPO_FOLDER}/poppy_uppsala-${POPPY_UU_VERSION}/build/references.hg38.md5sums.yaml}"
APPT_IMAGES="${APPT_IMAGES:-skip}"

# Check the options that only accept fixed values
case "${REF_DOWNLOAD,,}" in
    true|false) ;;
    *) echo "REF_DOWNLOAD must be true or false, got: ${REF_DOWNLOAD}" >&2; exit 1 ;;
esac

case "${APPT_IMAGES,,}" in
    skip|build) ;;
    *) echo "APPT_IMAGES must be skip or build, got: ${APPT_IMAGES}" >&2; exit 1 ;;
esac


# ------------------------------------------------------------------------------
# Fixed settings (not options): folder layout, file names and the module list
# ------------------------------------------------------------------------------
# All paths are relative to the directory the script is run from.
ENV_DIR="./${REPO_FOLDER}/venv_poppy_uppsala-${POPPY_UU_VERSION}"   # conda env
UU_DIR="./${REPO_FOLDER}/poppy_uppsala-${POPPY_UU_VERSION}"         # clone of poppy_uppsala
GMS_DIR="./${REPO_FOLDER}/poppy-${POPPY_GMS_VERSION}"               # clone of Poppy GMS
APPTAINER_CACHE_DIR="apptainer_cache"                               # where the .sif images are built

# Config files that are read ...
GMS_CONFIG="${GMS_DIR}/config/config.yaml"
UU_CONFIG="${UU_DIR}/config/config_marvin.yaml"
# ... and the offline versions that are created. These are used instead of the originals.
GMS_OFFLINE_CONFIG="${GMS_DIR}/config/config_poppy_offline.yaml"
UU_OFFLINE_CONFIG="${UU_DIR}/config/config_offline.yaml"
UU_TMP_CONFIG="${UU_DIR}/config/config_uu_offline.apptainer.yaml"   # intermediate, removed after use

# Profile that contains the REPO_FOLDER placeholder
#TODO: make it optional which profile to use incase of several offline profiles
PROFILE_CONFIG="${UU_DIR}/profile/miarka_drmaa/config.yaml"

# hydra-genetics modules the pipeline imports.
# NOTE: these are cloned from the default branch, i.e. not pinned to the versions used by the pipeline.
HYDRA_MODULES=(
    alignment annotation cnv_sv filtering prealignment
    qc reports snv_indels
)


# ------------------------------------------------------------------------------
# Step 1: conda environment
# ------------------------------------------------------------------------------
echo "Packaging Poppy UU ${POPPY_UU_VERSION} from ${POPPY_UU_REPO} and Poppy GMS ${POPPY_GMS_VERSION} from ${POPPY_GMS_REPO}"
echo "Reference validation file used: ${REFERENCE_VALIDATION_FILE}. conda env python ${PYTHON_VERSION}"

# Clones poppy_uppsala (its requirements.txt is needed), creates the env, installs the requirements
# and packs the env into env.tar.gz inside the repo folder. The env stays activated for the later steps.
mkdir -p "${REPO_FOLDER}"

eval "$(conda shell.bash hook)"
conda create --prefix "${ENV_DIR}" python="${PYTHON_VERSION}" -y
conda activate "${ENV_DIR}"
conda install -c conda-forge pip -y

git clone --branch "${POPPY_UU_VERSION}" "${POPPY_UU_REPO}" "${UU_DIR}"

"${ENV_DIR}/bin/pip3" install --no-cache-dir -I --upgrade "pip<25" wheel
"${ENV_DIR}/bin/pip3" install --no-cache-dir -I -r "${UU_DIR}/requirements.txt"

# needs conda-pack installed in the base conda
conda pack --prefix "${ENV_DIR}" -o "./${REPO_FOLDER}/env.tar.gz"


# ------------------------------------------------------------------------------
# Step 2: clone the remaining repositories
# ------------------------------------------------------------------------------
# Poppy GMS: only the chosen tag/branch
git clone --single-branch --branch "${POPPY_GMS_VERSION}" "${POPPY_GMS_REPO}" "${GMS_DIR}"

git clone https://github.com/snakemake/snakemake-wrappers.git "./${REPO_FOLDER}/snakemake-wrappers"

mkdir -p "./${REPO_FOLDER}/hydra-genetics"
for module in "${HYDRA_MODULES[@]}"; do
    git clone "https://github.com/hydra-genetics/${module}.git" "./${REPO_FOLDER}/hydra-genetics/${module}"
done

# The miarka profile has a placeholder for the folder the package is unpacked into on the cluster.
sed -i -E "s/REPO_FOLDER/${REPO_FOLDER}/g" "${PROFILE_CONFIG}"


# ------------------------------------------------------------------------------
# Step 3: apptainer images
# ------------------------------------------------------------------------------
# Builds one .sif per container in the configs into apptainer_cache/ (only with -a build).
# hydra-genetics ignores failed builds, so check the output for errors.
if [ "${APPT_IMAGES,,}" == "build" ]; then
    hydra-genetics prepare-environment create-singularity-files -c "${GMS_CONFIG}" -o "${APPTAINER_CACHE_DIR}"
    hydra-genetics prepare-environment create-singularity-files -c "${UU_CONFIG}" -o "${APPTAINER_CACHE_DIR}"
fi


# ------------------------------------------------------------------------------
# Step 4: offline configs
# ------------------------------------------------------------------------------
# Replaces the docker:// containers with {{APPTAINER_CACHE}}/<image>.sif.
hydra-genetics prepare-environment container-path-update \
    -c "${GMS_CONFIG}" -p '{{APPTAINER_CACHE}}' -n "${GMS_OFFLINE_CONFIG}"
# Intermediate file, the reference paths are updated in the next step
hydra-genetics prepare-environment container-path-update \
    -c "${UU_CONFIG}" -p '{{APPTAINER_CACHE}}' -n "${UU_TMP_CONFIG}"

# Replaces the marvin reference paths in the poppy_uppsala config with {{DESIGN_REF_FOLDER}}/...,
# the folder with the files downloaded in step 5.
# extract_reference_paths.py compares the paths in the config with the url:s in the validation file
# and prints one "--reference-path old:new" per folder. It stops with an error if a path is missing in the validation file.
# The output is used unquoted below so it is split into separate arguments, i.e. it expects
# folder paths without spaces or special characters. Assigned on its own line so that set -e stops on errors.
REFERENCE_FLAGS="$(python3 "${UU_DIR}/build/extract_reference_paths.py" \
    -v "${REFERENCE_VALIDATION_FILE}" \
    -c "${UU_TMP_CONFIG}" \
    -o flags)"

hydra-genetics prepare-environment reference-path-update \
    -c "${UU_TMP_CONFIG}" \
    -n "${UU_OFFLINE_CONFIG}" \
    ${REFERENCE_FLAGS}

# reference-path-update does not replace values inside lists, so merge_cnv_json/annotations is updated here.
# The value is quoted since a plain scalar starting with {{ is not valid yaml. 
# Add more if WARNING message from extract_reference_paths.py (step 4) gives other files to that should be included, then rerun the build script.
sed -E -i 's#^( *- )/projects/wp2/nobackup/Twist_Myeloid/Bed_files/HG38/(.*)$#\1"{{DESIGN_REF_FOLDER}}/\2"#' \
    "${UU_OFFLINE_CONFIG}"

rm -f "${UU_TMP_CONFIG}"


# ------------------------------------------------------------------------------
# Step 5: reference and design files
# ------------------------------------------------------------------------------
# Downloads the files in the validation file into the design/reference folder and packs it (only with -R true).
if [ "${REF_DOWNLOAD,,}" == "true" ]; then
    echo "Downloading references."
    # hydra-genetics only logs failed downloads and exits 0, so the log is checked.
    # The log is written next to the tar files, i.e. outside the folders that are packed.
    hydra-genetics --debug references download -o "${DESIGN_REF_FOLDER}" -v "${REFERENCE_VALIDATION_FILE}" 2>&1 \
        | tee reference_download.log
    if grep -q "Failed to retrieve" reference_download.log; then
        echo "Reference download had failures, see reference_download.log" >&2
        exit 1
    fi

    tar -zcvf "${DESIGN_REF_FOLDER}.tar.gz" "${DESIGN_REF_FOLDER}"
    rm -rf "${DESIGN_REF_FOLDER}"
fi


# ------------------------------------------------------------------------------
# Step 6: package
# ------------------------------------------------------------------------------
# Packs the repo folder (repos, modules, conda env and the offline configs) and removes the folder.
tar -zcvf "${REPO_FOLDER}.tar.gz" "${REPO_FOLDER}"
rm -rf "${REPO_FOLDER}"


# ------------------------------------------------------------------------------
# Done
# ------------------------------------------------------------------------------
conda deactivate
