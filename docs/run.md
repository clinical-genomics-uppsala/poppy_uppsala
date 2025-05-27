# Running poppy-uppsala

Execute poppy, then poppy uppsala.

## Run command
Using the activated python virtual environment created above, this is a basic command for running the pipeline:
```bash
snakemake --profile profiles/NAME_OF_PROFILE -s workflow/Snakefile
```  
<br />
The are many additional [snakemake running options](https://snakemake.readthedocs.io/en/stable/executing/cli.html#) some of which is listed below. However, options that are always used should be put in the [profile](https://hydra-genetics.readthedocs.io/en/latest/run_pipeline/profile/).

* --notemp - Saves all intermediate files. Good for development and testing different options.
* --until <rule> - Runs only rules dependent on the specified rule.

<br />
**Note:** Remember to have singularity and drmaa available on the system where the pipeline will be run.
<br />