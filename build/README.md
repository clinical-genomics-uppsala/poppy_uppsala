Build for packaging,
the script must be run from the directory where the pipeline will be packaged (`miarka_package` in the example below).

```bash
git clone https://github.com/clinical-genomics-uppsala/poppy_uppsala.git
mkdir miarka_package && cd miarka_package
bash ../poppy_uppsala/build/build_conda.sh poppy_uppsala/config/miarka/references/references.hg38.md5sums.yaml
```
