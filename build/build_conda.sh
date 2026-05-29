#!/usr/bin/env bash
set -e

eval "$(conda shell.bash hook)"

# ---- repos -----
PIPELINE_NAME="poppy_uppsala"
POPPY_UU_REPO="https://github.com/clinical-genomics-uppsala/poppy_uppsala.git"
POPPY_UU_VERSION="v2.0.0"
# poppy_uppsala is build on the top of poppy GMS
POPPY_GMS_VERSION="v3.0.0"
POPPY_GMS_REPO="https://github.com/genomic-medicine-sweden/poppy.git"
CONFIG_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/poppy_uppsala_config.git"
CONFIG_VERSION="develop"
PYTHON_VERSION="3.9"  # 3.11 funkar ocksa

# Clone git of poppy_uppsala to configure conda environment
git clone --branch ${POPPY_UU_VERSION} ${POPPY_UU_REPO}
cd ${PIPELINE_NAME}

# Create and activate conda environment in the current directory, then install pipeline requirements
conda create --prefix ./${PIPELINE_NAME}_${POPPY_UU_VERSION}_env python=${PYTHON_VERSION} -y
conda activate ./${PIPELINE_NAME}_${POPPY_UU_VERSION}_env
conda install -c conda-forge pip -y

if [ -d ${PIPELINE_NAME}_${POPPY_UU_VERSION} ];
then
    rm -fr ${PIPELINE_NAME}_${POPPY_UU_VERSION}
fi

# The directory ${PIPELINE_NAME}_${POPPY_UU_VERSION} is created for the files that are to be packaged and transferred
# elsewhere:
# - the pipeline code for poppy_uppsala as well as the "base" code from Poppy GMS
# - the conda environment used to run the pipeline
# - the snakemake-wrappers
# - the hydra-genetics modules
mkdir -p ${PIPELINE_NAME}_${POPPY_UU_VERSION}

# Clone git of poppy_uppsala and package the env
git clone --branch ${POPPY_UU_VERSION} ${POPPY_UU_REPO} ${PIPELINE_NAME}_${POPPY_UU_VERSION}/${PIPELINE_NAME}
./${PIPELINE_NAME}_${POPPY_UU_VERSION}_env/bin/pip3 install --no-cache-dir -I -r ${PIPELINE_NAME}_${POPPY_UU_VERSION}/${PIPELINE_NAME}/requirements.txt
conda pack --prefix ./${PIPELINE_NAME}_${POPPY_UU_VERSION}_env -o ${PIPELINE_NAME}_${POPPY_UU_VERSION}/env.tar.gz

# Clone the relevant branch only of Poppy GMS
git clone --single-branch --branch ${POPPY_GMS_VERSION} ${POPPY_GMS_REPO} ${PIPELINE_NAME}_${POPPY_UU_VERSION}/poppy


# Clone wrappers
git clone https://github.com/snakemake/snakemake-wrappers.git ${PIPELINE_NAME}_${POPPY_UU_VERSION}/snakemake-wrappers

# Clone hydra modules
mkdir -p ${PIPELINE_NAME}_${POPPY_UU_VERSION}/hydra-genetics
git clone https://github.com/hydra-genetics/alignment.git ${PIPELINE_NAME}_${POPPY_UU_VERSION}/hydra-genetics/alignment
git clone https://github.com/hydra-genetics/annotation.git ${PIPELINE_NAME}_${POPPY_UU_VERSION}/hydra-genetics/annotation
git clone https://github.com/hydra-genetics/cnv_sv.git ${PIPELINE_NAME}_${POPPY_UU_VERSION}/hydra-genetics/cnv_sv
git clone https://github.com/hydra-genetics/filtering.git ${PIPELINE_NAME}_${POPPY_UU_VERSION}/hydra-genetics/filtering
git clone https://github.com/hydra-genetics/fusions.git ${PIPELINE_NAME}_${POPPY_UU_VERSION}/hydra-genetics/fusions
git clone https://github.com/hydra-genetics/prealignment.git ${PIPELINE_NAME}_${POPPY_UU_VERSION}/hydra-genetics/prealignment
git clone https://github.com/hydra-genetics/qc.git ${PIPELINE_NAME}_${POPPY_UU_VERSION}/hydra-genetics/qc
git clone https://github.com/hydra-genetics/reports.git ${PIPELINE_NAME}_${POPPY_UU_VERSION}/hydra-genetics/reports
git clone https://github.com/hydra-genetics/snv_indels.git ${PIPELINE_NAME}_${POPPY_UU_VERSION}/hydra-genetics/snv_indels
git clone https://github.com/hydra-genetics/references.git ${PIPELINE_NAME}_${POPPY_UU_VERSION}/hydra-genetics/references # not used every time only when building pon

#  Replace POPPY_UU_VERSION with the chosen version of the pipeline for snakemake wrapper and hydragenetics paths
sed -i -E "s/POPPY_UU_VERSION/${POPPY_UU_VERSION}/g" ./${PIPELINE_NAME}_${POPPY_UU_VERSION}/${PIPELINE_NAME}/profile/miarka/config.yaml
sed -i -E "s/POPPY_UU_VERSION/${POPPY_UU_VERSION}/g" ./${PIPELINE_NAME}_${POPPY_UU_VERSION}/${PIPELINE_NAME}/config/config_offline.yaml # borde kanske vara separat med bara miarka paths?

# # Download containers and update container path in configs
if [ ${APPT_CACHE_STATUS} == "build" ];
then
    hydra-genetics prepare-environment create-singularity-files -c ${PIPELINE_NAME}_${POPPY_UU_VERSION}/poppy/config/config_static.yaml -o apptainer_cache
    cp ${PIPELINE_NAME}_${POPPY_UU_VERSION}/poppy/config/config_static.yaml ${PIPELINE_NAME}_${POPPY_UU_VERSION}/poppy/config/config_static.yaml.copy
    hydra-genetics prepare-environment container-path-update -c ${PIPELINE_NAME}_${POPPY_UU_VERSION}/poppy/config/config_static.yaml.copy -n ${PIPELINE_NAME}_${POPPY_UU_VERSION}/poppy/config/config_static.yaml -p ${PATH_TO_apptainer_cache}
    # config_marvin from poppy_uppsala also need to be updated
    hydra-genetics prepare-environment create-singularity-files -c ${PIPELINE_NAME}_${POPPY_UU_VERSION}/${PIPELINE_NAME}/config/config_marvin.yaml -o apptainer_cache
    cp ${PIPELINE_NAME}_${POPPY_UU_VERSION}/${PIPELINE_NAME}/config/config_marvin.yaml ${PIPELINE_NAME}_${POPPY_UU_VERSION}/${PIPELINE_NAME}/config/config_marvin.yaml.copy
    hydra-genetics prepare-environment container-path-update -c ${PIPELINE_NAME}_${POPPY_UU_VERSION}/${PIPELINE_NAME}/config/config_marvin.yaml.copy -n ${PIPELINE_NAME}_${POPPY_UU_VERSION}/${PIPELINE_NAME}/config/config_marvin.yaml -p ${PATH_TO_apptainer_cache}
fi

# If apptainer already available remote and no new build is needed (specified by user). Only update path to cache in config.
if [ ${APPT_CACHE_STATUS} == "update" ];
then
    cp ${PIPELINE_NAME}_${POPPY_UU_VERSION}/poppy/config/config_static.yaml ${PIPELINE_NAME}_${POPPY_UU_VERSION}/poppy/config/config_static.yaml.copy
    hydra-genetics prepare-environment container-path-update -c ${PIPELINE_NAME}_${POPPY_UU_VERSION}/poppy/config/config_static.yaml.copy -n ${PIPELINE_NAME}_${POPPY_UU_VERSION}/poppy/config/config_static.yaml -p ${PATH_TO_apptainer_cache}
    cp ${PIPELINE_NAME}_${POPPY_UU_VERSION}/${PIPELINE_NAME}/config/config_marvin.yaml ${PIPELINE_NAME}_${POPPY_UU_VERSION}/${PIPELINE_NAME}/config/config_marvin.yaml.copy
    hydra-genetics prepare-environment container-path-update -c ${PIPELINE_NAME}_${POPPY_UU_VERSION}/${PIPELINE_NAME}/config/config_marvin.yaml.copy -n ${PIPELINE_NAME}_${POPPY_UU_VERSION}/${PIPELINE_NAME}/config/config_marvin.yaml -p ${PATH_TO_apptainer_cache}
fi

# # Pack all cloned repositories
tar -zcvf ${PIPELINE_NAME}_${POPPY_UU_VERSION}.tar.gz ${PIPELINE_NAME}_${POPPY_UU_VERSION}


# Download references if given on command line
if [ "$#" != 0 ];
then
    for reference_config in "$@"
    do
        hydra-genetics --debug references download -o design_and_ref_files -v $reference_config
    done
    # Compress data
    tar -czvf design_and_ref_files.tar.gz design_and_ref_files
fi


conda deactivate

if [ -d ${PIPELINE_NAME}_${POPPY_UU_VERSION}_env ];
then
    rm -fr ${PIPELINE_NAME}_${POPPY_UU_VERSION}_env
fi

if [ -d ${PIPELINE_NAME}_${POPPY_UU_VERSION} ];
then
    rm -fr ${PIPELINE_NAME}_${POPPY_UU_VERSION}
fi
