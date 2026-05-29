Build for packaging,
the script must be run from the directory where the pipeline will be packaged (`miarka_package` in the example below).


```bash
module load miniconda3
git clone https://github.com/clinical-genomics-uppsala/poppy_uppsala.git
mkdir miarka_package && cd miarka_package

POPPY_GMS_VERSION="v3.0.0" POPPY_GMS_REPO="https://github.com/genomic-medicine-sweden/poppy.git" \
POPPY_UU_VERSION="v2.0.0" PIPELINE_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/poppy_uppsala.git" \
PIPELINE_NAME="poppy_uppsala" PYTHON_VERSION="3.9" \
PATH_TO_apptainer_cache="path/to/.sif/on/remote" \
APPT_CACHE_STATUS="build"

bash build_conda.sh poppy_uppsala_config/config/miarka/references/references.GRCh38.md5sums.yaml poppy_uppsala_config/config/miarka/references/references.PoN.md5sums.yaml
```
!!! note
If you run on the compute node:
export SINGULARITY_CACHEDIR=/projects/wp4/nobackup/singularity_cache/{user}

Copy the following files and folders to the cluster (eg Miarka):
* poppy_uppsala_{version}.tar.gz
* design_and_ref_files.tar.gz
* apptainer_cache


Next step: Extract the tar files and the env file inside the poppy_uppsala_{version} folder:
```bash
cd poppy_uppsala_{version}
mkdir poppy_uppsala/venv/
tar -xzvf env.tar.gz -C poppy_uppsala/venv/
```