# Running poppy-uppsala

Once you are done with the setup and the configuration of the pipeline, you can proceed with running it.

## Run command
Using the activated python virtual environment created above, this is a basic command to run the pipeline in dry-run mode:

```bash
snakemake -n --profile profiles/NAME_OF_PROFILE -s workflow/Snakefile --configfiles NAME_OF_CONFIG_POPPY NAME_OF_CONFIG_POPPY_UPPSALA --config POPPY_HOME=. POPPY_UU_HOME=. sequenceid="dry_run" poppy_version="v0.2.0" poppy_uu_version="0.1.2"
```  

**Important note**: The configuration file for poppy **must** be specified first after the option `--configfiles`,
before the configuration file to use for poppy-uppsala as the latter file overrides some values in the first file.

### Config files and dictionary

That said, it means that if some keys are identical in the config files for Poppy GMS and poppy-uppsala, 
the values in the poppy-uppsala config are the ones that are actually used during the execution.

*Example:*

```yaml
# Poppy GMS
pindel2vcf:
  container: "docker://hydragenetics/pindel:0.2.5b9"
  extra: "-e 10 -mc 10 -is 5 -he 0.01 -G" # num reads support
  refname: "hg19"
  refdate: "2009"
```
and
```yaml
# poppy-uppsala
pindel2vcf:
  extra: "-e 10 -mc 10 -is 5 -he 0.01 -G" # num reads support
  refname: "GCA_000001405.15_GRCh38_no_alt"
  refdate: "2022"
```
which is aggregated as follows in the combined pipeline that is actually run:
```python3
# poppy GMS + poppy-uppsala
config_dict["pindel2vcf"] = {
  "container": "docker://hydragenetics/pindel:0.2.5b9",
  "extra": "-e 10 -mc 10 -is 5 -he 0.01 -G", 
  "refname": "GCA_000001405.15_GRCh38_no_alt",
  "refdate": "2022"
}
```


<br />
The are many additional [snakemake running options](https://snakemake.readthedocs.io/en/stable/executing/cli.html#) some of which is listed below. However, options that are always used should be put in the [profile](https://hydra-genetics.readthedocs.io/en/latest/run_pipeline/profile/).

* --notemp - Saves all intermediate files. Good for development and testing different options.
* --until <rule> - Runs only rules dependent on the specified rule.

To see an example of valid command to run the pipeline, see the last line of `.github/workflows/snakemake-dry-run.yaml`
in the repository on GitHub.

<br />
**Note:** Remember to have singularity and drmaa available on the system where the pipeline will be run.
<br />

## Example of custom start script used at Clinical Genomics Uppsala

See the [bash script](https://github.com/clinical-genomics-uppsala/pipeline_start_scripts/blob/develop/marvin/start_wp2_tm.sh) 
that we use on our HPC cluster at Clinical Genomics Uppsala.
You may use something similar to start your pipeline

### Fetch the source code to be used
The script clones the branch or tag of the Poppy GMS/Uppsala repositories that are chosen by the user
and then install the required virtual environments.

### Create the input files for hydra-genetics
The script determines what kind of sequencing machine was used to generate the data. 
The type of the machine is required to pass the correct options to the command `hydra-genetics create-input-files` 
that creates `samples.tsv` and `input.tsv`. 
For instance, one must pass the option `--default-barcode NNNNNNNN+NNNNNNNN` for NextSeq machines,
but not for NovaSeq machines.
Among other things, the input files let hydra-genetics know what wildcards to use for the samples in the analysis.

### Pick the correct configuration file
Different sequencing machines lead to different configurations to be used e.g. different artifact files, at CGU we 
chose the solution to have different configuration files for each machine type. 
After parsing the machine type, the script copies the relevant configuration file in the repository that was cloned.
At CGU, we have a configuration file for NextSeq machines (`config/config_uppsala_nextseq.yaml`) 
and one for NovaSeq machines (`config/config_uppsala_novaseq.yaml`).
The configuration file for Poppy GMS is the same regardless of the machine as we can override the keys if needed in 
in configurations of poppy-uppsala.

### Run the pipeline with Snakemake
Last, the scripts starts the analysis of the sequence data with Snakemake.
You may generate a rule graph to see the dependencies between the rules and the order in which they are executed.
Below is the command we would use in our CGU script:

```bash
...
snakemake --profile ${poppy_uu_path}poppy_uppsala/profile -s ${poppy_uu_path}poppy_uppsala/workflow/Snakefile \
        --configfiles config.yaml config_uppsala.yaml --config poppy_version=${poppy_version} poppy_uu_version=${poppy_uu_version} POPPY_HOME=${poppy_path}poppy \
        POPPY_UU_HOME=${poppy_uu_path}poppy_uppsala sequenceid=${sequenceid} --forceall --rulegraph | dot -Tpdf > rulegraph.pdf
```

Note that this command will create a log file in the `.snakemake/log/` directory.
This file contains only the printout "Building DAG of jobs...".

## Run time
The run time of the analysis depends on:
- the number of samples,
- the sequencing depth,
- the computational resources that are available.

At CGU, the analysis usually takes a couple of hours for 16 samples with 700-800x coverage on a NextSeq machine.