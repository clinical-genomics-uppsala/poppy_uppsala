__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2024, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"

import itertools
import numpy as np
import pandas as pd
import pathlib
import re
from snakemake.utils import validate
from snakemake.utils import min_version
import yaml
from datetime import datetime

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *
from hydra_genetics import min_version as hydra_min_version
from hydra_genetics.utils.software_versions import add_version_files_to_multiqc
from hydra_genetics.utils.software_versions import export_pipeline_version_as_file
from hydra_genetics.utils.software_versions import get_pipeline_version
from hydra_genetics.utils.software_versions import touch_pipeline_version_file_name
from hydra_genetics.utils.misc import replace_dict_variables
from hydra_genetics.utils.misc import export_config_as_file

min_version("6.8.0")

## Version logging for MultiQC
date_string = datetime.now().strftime("%Y%m%d")

# Create pipeline version file, populate since only one onstart *(loaded last)
pipeline_version = get_pipeline_version(workflow, pipeline_name="poppy_uppsala")
version_files = touch_pipeline_version_file_name(pipeline_version, date_string=date_string, directory="versions/software")
add_version_files_to_multiqc(config, version_files)
export_pipeline_version_as_file(pipeline_version, date_string=date_string, directory="versions/software")

### Set and validate config file

if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/--configfiles, by command line or a profile!")
config = replace_dict_variables(config)

try:
    validate(config, schema="../schemas/config.schema.yaml")
except WorkflowError as we:
    # Probably a validation error, but the original exception in lost in
    # snakemake. Pull out the most relevant information instead of a potentially
    # *very* long error message.
    if not we.args[0].lower().startswith("error validating config file"):
        raise
    error_msg = "\n".join(we.args[0].splitlines()[:2])
    parent_rule_ = we.args[0].splitlines()[3].split()[-1]
    if parent_rule_ == "schema:":
        sys.exit(error_msg)
    else:
        schema_hiearachy = parent_rule_.split()[-1]
        schema_section = ".".join(re.findall(r"\['([^']+)'\]", schema_hiearachy)[1::2])
        sys.exit(f"{error_msg} in {schema_section}")

### Read and validate resources file

config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], comment="#").set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = (
    pandas.read_table(config["units"], dtype=str)
    .set_index(["sample", "type", "flowcell", "lane", "barcode"], drop=False)
    .sort_index()
)

validate(units, schema="../schemas/units.schema.yaml")

### Read and validate output file

with open(config["output"]) as output:
    if config["output"].endswith("json"):
        output_spec = json.load(output)
    elif config["output"].endswith("yaml") or config["output"].endswith("yml"):
        output_spec = yaml.safe_load(output.read())

validate(output_spec, schema="../schemas/output_files.schema.yaml")

first_sample = config.get("first_sample", "")
last_sample = config.get("last_sample", "")
### Set wildcard constraints
if first_sample != "" and last_sample != "":

    wildcard_constraints:
        sample="|".join(samples.index[int(first_sample) : int(last_sample)]),
        type="N|T|R",

else:

    wildcard_constraints:
        sample="|".join(samples.index),
        type="N|T|R",


def compile_output_file_list(wildcards):
    outdir = pathlib.Path(output_spec["directory"])
    output_files = []

    callers = config["bcbio_variation_recall_ensemble"]["callers"]
    wc_df = pd.DataFrame(np.repeat(units.values, len(callers), axis=0))
    wc_df.columns = units.columns
    caller_gen = itertools.cycle(callers)
    wc_df = wc_df.assign(caller=[next(caller_gen) for i in range(wc_df.shape[0])])
    wc_df = wc_df.assign(sequenceid=[config["sequenceid"] for i in range(wc_df.shape[0])])

    for f in output_spec["files"]:
        outputpaths = set(expand(f["output"], zip, **wc_df.to_dict("list")))
        if len(outputpaths) == 0:
            # Using expand with zip on a pattern without any wildcards results
            # in an empty list. Then just add the output filename as it is.
            outputpaths = [f["output"]]
        for op in outputpaths:
            output_files.append(outdir / Path(op))

    return output_files


def generate_copy_rules(output_spec):
    output_directory = pathlib.Path(output_spec["directory"])
    rulestrings = []
    rule_code = ""

    for f in output_spec["files"]:
        if f["input"] is None:
            continue

        rule_name = "copy_{}".format("_".join(re.sub(r"[\"'-.,]", "", f["name"].strip().lower()).split()))
        input_file = pathlib.Path(f["input"])
        output_file = output_directory / pathlib.Path(f["output"])

        mem_mb = config.get("copy", {}).get("mem_mb", config["default_resources"]["mem_mb"])
        mem_per_cpu = config.get("copy", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"])
        partition = config.get("copy", {}).get("partition", config["default_resources"]["partition"])
        threads = config.get("copy", {}).get("threads", config["default_resources"]["threads"])
        time = config.get("copy", {}).get("time", config["default_resources"]["time"])
        copy_container = config.get("_copy", {}).get("container", config["default_container"])

        rule_code += f'@workflow.rule(name="{rule_name}")\n'
        rule_code += f'@workflow.input("{input_file}")\n'
        if f["output"].endswith("/"):
            rule_code += f'@workflow.output(directory("{output_file}"))\n'
        else:
            rule_code += f'@workflow.output("{output_file}")\n'

        rule_code += f'@workflow.log("logs/{rule_name}_{output_file.name}.log")\n'
        rule_code += f'@workflow.container("{copy_container}")\n'
        rule_code += f'@workflow.resources(time = "{time}", threads = {threads}, mem_mb = {mem_mb}, mem_per_cpu = {mem_per_cpu}, partition = "{partition}")\n'
        rule_code += '@workflow.shellcmd("cp -r {input} {output}")\n\n'
        rule_code += "@workflow.run\n"
        rule_code += (
            f"def __rule_{rule_name}(input, output, params, wildcards, threads, resources, log, version, rule, "
            "conda_env, container_img, singularity_args, use_singularity, env_modules, bench_record, jobid, is_shell, "
            "bench_iteration, cleanup_scripts, shadow_dir, edit_notebook, conda_base_path, basedir, runtime_sourcecache_path, "
            "__is_snakemake_rule_func=True):\n"
            '\tshell ( "(cp -r {input[0]} {output[0]}) &> {log}" , bench_record=bench_record, bench_iteration=bench_iteration)\n\n'
        )

    exec(compile(rule_code, "copy_result_files", "exec"), workflow.globals)


generate_copy_rules(output_spec)
