# Running poppy-uppsala

## Requirements

...

## Installation

### Clone the repositories 

We recommend that the repositories are cloned to your working directory. 

* Poppy:

```bash
git clone git@github.com:genomic-medicine-sweden/poppy.git
```

* Poppy Uppsala:

```bash
git clone git@github.com:clinical-genomics-uppsala/poppy_uppsala.git
```

### Create the python virtual environment (with the reauirements)

...

### Edit the configurations

The values in `config/NAME_OF_CONFIG` must match the resources you have, 
for instance the paths to BED files and to the reference genome hg38 fasta-file.

## Create the input sample files

...

## Run command
Using the activated python virtual environment created above, this is a basic command to run the pipeline in dry-run mode:

```bash
snakemake -n --profile profiles/NAME_OF_PROFILE -s workflow/Snakefile --configfiles NAME_OF_CONFIG_POPPY NAME_OF_CONFIG_POPPY_UPPSALA --config POPPY_HOME=. POPPY_UU_HOME=. sequenceid="dry_run" poppy_version="v0.2.0" poppy_uu_version="0.1.2"
```  

The configuration file for poppy **must** be specified first after the option `--configfiles`,
before the configuration file to use for poppy-uppsala as the latter file overrides some values in the first file.

<br />
The are many additional [snakemake running options](https://snakemake.readthedocs.io/en/stable/executing/cli.html#) some of which is listed below. However, options that are always used should be put in the [profile](https://hydra-genetics.readthedocs.io/en/latest/run_pipeline/profile/).

* --notemp - Saves all intermediate files. Good for development and testing different options.
* --until <rule> - Runs only rules dependent on the specified rule.

To see an example of valid command to run the pipeline, see the last line of `.github/workflows/snakemake-dry-run.yaml`
in the repository on GitHub.

<br />
**Note:** Remember to have singularity and drmaa available on the system where the pipeline will be run.
<br />