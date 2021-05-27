configfile: "configuration/snakemake.yaml"
# container: "docker://condaforge/miniforge:4.10.1-0"

BENCHS = config["bench_dir"]
ENVS = config["env_dir"]
LOGS = config["log_dir"]
RESULTS = config["results_dir"]

rule all:
        input:
                RESULTS + "process/sorted_isoforms.xlsx",

rule request:
        params:
                gene_ids = config["request"]["gene_ids"],
        output:
                data = expand(RESULTS + "request/{gene}_message.tsv", gene=config["request"]["gene_ids"]),
        log:
                LOGS + "request.log"
        benchmark:
                BENCHS + "request.txt"
        threads:
                6
        conda:
                ENVS + "request.yml"
        script:
                "scripts/request.py"

rule process:
        input:
                data = expand(RESULTS + "request/{gene}_message.tsv", gene=config["request"]["gene_ids"]),
        output:
                data = RESULTS + "process/sorted_isoforms.xlsx",
        log:
                LOGS + "process.log"
        benchmark:
                BENCHS + "process.txt"
        threads:
                6
        conda:
                ENVS + "process.yml"
        script:
                "scripts/process.py"
