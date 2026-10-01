Build for packaging,
the script must be run from the directory where the pipeline will be packaged (`miarka_package` in the example below).

```bash
module load miniconda3
git clone https://github.com/clinical-genomics-uppsala/poppy_uppsala.git
mkdir miarka_package && cd miarka_package

POPPY_GMS_VERSION="v1.0.0" POPPY_GMS_REPO="https://github.com/genomic-medicine-sweden/poppy.git" \
TAG_OR_BRANCH="v0.4.1" CONFIG_VERSION="develop" PIPELINE_NAME="poppy_uppsala" PYTHON_VERSION="3.9" \
PIPELINE_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/poppy_uppsala.git" \
CONFIG_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/poppy_uppsala_config.git" \
bash ../poppy_uppsala/build/build_conda.sh poppy_uppsala_config/config/miarka/references/references.hg38.md5sums.yaml
```
