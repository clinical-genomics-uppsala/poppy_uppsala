__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2024, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"

from datetime import datetime
from hydra_genetics.utils.software_versions import get_pipeline_version
from hydra_genetics.utils.software_versions import touch_pipeline_version_file_name


date_string = datetime.now().strftime("%Y%m%d")
pipeline_version = get_pipeline_version(workflow, pipeline_name="poppy")
poppy_yaml = touch_pipeline_version_file_name(pipeline_version, date_string=date_string, directory="versions/software")
pipeline_uu_version = get_pipeline_version(workflow, pipeline_name="poppy_uppsala")
poppy_uu_yaml = touch_pipeline_version_file_name(pipeline_uu_version, date_string=date_string, directory="versions/software")


rule version_update_poppy:
    input:
        yaml=poppy_yaml,
        yaml_uu=poppy_uu_yaml,
    output:
        touch_file=temp(touch("versions/update_poppy.temp")),
        yaml=f"versions/software_{date_string}/poppy_mqc_versions.yaml",
        yaml_uu=f"versions/software_{date_string}/poppy_uppsala_mqc_versions.yaml",
    params:
        poppy_version=config["poppy_version"],
        poppy_uu_version=config["poppy_uu_version"]
    log:
        "versions/update_poppy.temp.log",
    benchmark:
        repeat("versions/update_poppy.temp.benchmark.tsv", config.get("version_update_poppy", {}).get("benchmark_repeats", 1))
    threads: config.get("version_update_poppy", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("version_update_poppy", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("version_update_poppy", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("version_update_poppy", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("version_update_poppy", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("version_update_poppy", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("version_update_poppy", {}).get("container", config["default_container"])
    localrule: True
    message:
        "{rule}: update poppy_version for multiqc since it now has poppy_uppsala version"
    shell:
        """
        echo 'poppy: {params.poppy_version}'>{output.yaml}
        echo 'poppy_uppsala: {params.poppy_uu_version}'>{output.yaml_uu}
        rm -f {input.yaml} {input.yaml_uu}
        """
