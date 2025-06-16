# Setup and configurations

Remember you need both 

Overload some config parameters: 
- Mosdepth coverage in exon regions only,
- Home folder of the analysis

Use a bash script to start the analysis.

## Requirements
**Recommended hardware**

- CPU: >10 cores per sample
- Memory: 6GB per core
- Storage: >75GB per sample

**Note**: Running the pipeline with less resources may work, but has not been tested.

**Software**

- [python](https://www.python.org/), version 3.9 or newer
- [pip3](https://pypi.org/project/pip/)
- [virtuelenv](https://docs.python.org/3/library/venv.html)
- [singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html)

**Nice to have**

- DRMAA compatible scheduler
- valid [SSH key to Github](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)

## Installation
A list of releases of the GMS Poppy pipeline can be found at: [Releases](https://github.com/genomic-medicine-sweden/poppy/releases).
A list of releases of the poppy-uppsala pipeline can be found at: [Releases](https://github.com/clinical-genomics-uppsala/poppy_uppsala/releases).

### Clone both the GMS Poppy repo and the poppy-uppsala git repo
We recommend that the poppy-uppsala repository is cloned to your working directory, on the same level as the GMS Poppy repository. 
```bash
# Set up a working directory path
WORKING_DIRECTORY="/path_working_to_directory"
```

Fetch pipeline
Choose the release you need for both GMS Poppy and poppy-uppsala, for instance `0.2.0` and `v0.2.1`:
```bash
# Set versions
VERSION="v0.2.0"
VERSION_UU="v0.2.1"

# Clone selected version, use SSH URL if you have configured a local SSH key to Github (preferred)
git clone --branch ${VERSION} https://github.com/genomic-medicine-sweden/poppy.git ${WORKING_DIRECTORY}
git clone --branch ${VERSION_UU} https://github.com/clinical-genomics-uppsala/poppy_uppsala.git ${WORKING_DIRECTORY}
```

### Create python environment
To run the poppy-uppsala pipeline, a python virtual environments is needed.
```bash
# Enter working directory
cd ${WORKING_DIRECTORY}

# Create a new virtual environment
python3 -m venv ${WORKING_DIRECTORY}/venv-poppy-uu
```

### Install pipeline requirements
Activate the virtual environment and install pipeline requirements specified in `requirements.txt`.
```bash
# Enter working directory
cd ${WORKING_DIRECTORY}

# Activate python environment
source venv-poppy-uu/bin/activate

# Install requirements
pip install -r requirements.txt
```

### Setup required data and config

**Download data**
```bash
# make sure hydra-genetics is available
# make sure that TMPDIR points to a location with a lot of storage, it
# will be required to fetch reference data
# export TMPDIR=/PATH_TO_STORAGE
hydra-genetics --debug --verbose references download -o design_and_ref_files  -v config/references/references.hg19.yaml -v config/references/design_files.hg19.yaml -v config/references/nextseq.hg19.pon.yaml
```

**Update config**
```yaml 
# file config/config.data.hg19.yaml
# change rows:
PROJECT_DESIGN_DATA: "PATH_TO/design_and_ref_files" # parent folder for GMS560 design, ex GMS560/design
PROJECT_PON_DATA: "PATH_TO/design_and_ref_files" # artifact/background/PoN, ex GMS560/PoN
PROJECT_REF_DATA: "PATH_TO/design_and_ref_files" # parent folder for ref_data, ex ref_data/hg19
```

poppy-uppsala overwrites some configurations of poppy.

NB: tricky versioning for MultiQC.

## Input sample files
The pipeline uses sample input files (`samples.tsv` and `units.tsv`) with information regarding sample information, sequencing meta information as well as the location of the fastq-files. Specification for the input files can be found at [Twist Solid schemas](https://github.com/genomic-medicine-sweden/Twist_Solid/blob/develop/workflow/schemas/). Using the python virtual environment created above it is possible to generate these files automatically using [hydra-genetics create-input-files](https://hydra-genetics.readthedocs.io/en/latest/run_pipeline/create_sample_files/):
```bash
hydra-genetics create-input-files -d path/to/fastq-files/
```

## Known variants and custom filtering
Specific table with known variants, artifacts (machine-specific), pindel regions to limit the computational cost,
