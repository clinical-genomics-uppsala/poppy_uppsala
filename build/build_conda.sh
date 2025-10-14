#!/usr/bin/env bash
set -e

eval "$(conda shell.bash hook)"

# Clone git of poppy_uppsala to configure conda environment
git clone --branch ${TAG_OR_BRANCH} ${PIPELINE_GITHUB_REPO}
cd ${PIPELINE_NAME}

# Create and activate conda envrionment in the current directory, then install pipeline requirements
conda create --prefix ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env python=${PYTHON_VERSION} -y
conda activate ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env
conda install -c conda-forge pip -y

if [ -d ${PIPELINE_NAME}_${TAG_OR_BRANCH} ];
then
    rm -fr ${PIPELINE_NAME}_${TAG_OR_BRANCH}
fi

# The directory ${PIPELINE_NAME}_${TAG_OR_BRANCH} is created the files to be packaged and transferred elsewhere:
# - the pipeline code for poppy_uppsala as well as the "base" code from Poppy GMS
# - the conda environment used to run the pipeline
# - the snakemake-wrappers
# - the hydra-genetics modules
mkdir -p ${PIPELINE_NAME}_${TAG_OR_BRANCH}

# Clone git of poppy_uppsala
git clone --branch ${TAG_OR_BRANCH} ${PIPELINE_GITHUB_REPO} ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}
./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env/bin/pip3 install -r ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/requirements.txt
conda pack --prefix ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env -o ${PIPELINE_NAME}_${TAG_OR_BRANCH}/env.tar.gz

# Clone git of Poppy GMS
git clone https://github.com/genomic-medicine-sweden/poppy.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/poppy

# Clone snakemake-wrappers and hydra-genetics
mkdir -p ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics

git clone https://github.com/snakemake/snakemake-wrappers.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/snakemake-wrappers

git clone https://github.com/hydra-genetics/alignment.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/alignment
git clone https://github.com/hydra-genetics/annotation.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/annotation
git clone https://github.com/hydra-genetics/biomarker.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/biomarker
git clone https://github.com/hydra-genetics/cnv_sv.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/cnv_sv
git clone https://github.com/hydra-genetics/compression.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/compression
git clone https://github.com/hydra-genetics/filtering.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/filtering
git clone https://github.com/hydra-genetics/fusions.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/fusions
git clone https://github.com/hydra-genetics/misc.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/misc
git clone https://github.com/hydra-genetics/mitochondrial ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/mitochondrial
git clone https://github.com/hydra-genetics/parabricks ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/parabricks
git clone https://github.com/hydra-genetics/prealignment.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/prealignment
git clone https://github.com/hydra-genetics/qc.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/qc
git clone https://github.com/hydra-genetics/reports.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/reports
git clone https://github.com/hydra-genetics/snv_indels.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/snv_indels
git clone https://github.com/hydra-genetics/references.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/references

## Download the config files from the config repo
git clone --branch ${CONFIG_VERSION} ${CONFIG_GITHUB_REPO} ./poppy_uppsala_config
## copy resources.yaml files to the pipline config directory
cp poppy_uppsala_config/config/*.yaml ./${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/config/
cp -r poppy_uppsala_config/references ./${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/config/
cp -r poppy_uppsala_config/profiles/* ./${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/profiles/
## replace TAG_OR_BRANCH in profiles by the chosen version of the pipeline
sed -i -E "s/TAG_OR_BRANCH/${TAG_OR_BRANCH}/g" ./${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/profiles/miarka/config.yaml

# Pack all cloned repositories
tar -zcvf ${PIPELINE_NAME}_${TAG_OR_BRANCH}.tar.gz ${PIPELINE_NAME}_${TAG_OR_BRANCH}

# Download containers
conda activate ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env
hydra-genetics prepare-environment create-singularity-files -c config/config.yaml -o apptainer_cache

# Download references
for reference_config in "$@"
do
    hydra-genetics --debug references download -o design_and_ref_files -v $reference_config
done

conda deactivate

if [ -d ${PIPELINE_NAME}_${TAG_OR_BRANCH}_env ];
then
    rm -fr ${PIPELINE_NAME}_${TAG_OR_BRANCH}_env
fi

if [ -d ${PIPELINE_NAME}_${TAG_OR_BRANCH} ];
then
    rm -fr ${PIPELINE_NAME}_${TAG_OR_BRANCH}
fi

# Compress data
tar -czvf design_and_ref_files.tar.gz design_and_ref_files
