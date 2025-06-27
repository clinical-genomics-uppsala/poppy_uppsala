# Running poppy-uppsala

Once you are done with the setup and the configuration of the pipeline, you can proceed with running it.

## Run command
Using the activated python virtual environment created above, this is a basic command to run the pipeline in dry-run mode:

```bash
snakemake -n --profile profiles/NAME_OF_PROFILE -s workflow/Snakefile --configfiles NAME_OF_CONFIG_POPPY NAME_OF_CONFIG_POPPY_UPPSALA --config POPPY_HOME=. POPPY_UU_HOME=. sequenceid="dry_run" poppy_version="v0.2.0" poppy_uu_version="0.1.2"
```  

**Important note**: The configuration file for poppy **must** be specified first after the option `--configfiles`,
before the configuration file to use for poppy-uppsala as the latter file overrides some values in the first file.

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