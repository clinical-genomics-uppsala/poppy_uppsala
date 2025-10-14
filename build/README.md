Build for packaging:

```bash
TAG_OR_BRANCH="miarka" CONFIG_VERSION="develop" PIPELINE_NAME="poppy_uppsala" PYTHON_VERSION="3.9" \
PIPELINE_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/poppy_uppsala.git" \
CONFIG_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/poppy_uppsala_config.git" \
bash build_conda.sh poppy_uppsala/config/miarka/references/references.hg38.md5sums.yaml
```
