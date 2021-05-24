configfile: "configuration/snakemake.yaml"
container: "docker://condaforge/miniforge:4.10.1-0"

BENCHS = config["bench_dir"]
ENVS = config["env_dir"]
LOGS = config["log_dir"]

rule all:
        input:

rule request:
        params:
                gene_ids = config["request"]["gene_ids"],
        log:
                LOGS + "request.log"
        benchmark:
                BENCHS + "request.txt"
        conda:
                ENVS + "request.yml"
        script:
                "scripts/request.py"
