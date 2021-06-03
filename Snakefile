configfile: "configuration/snakemake.yaml"
container: "docker://condaforge/mambaforge:4.10.1-0"

BENCHS = config["bench_dir"]
DATA = config["data_dir"]
ENVS = config["env_dir"]
LOGS = config["log_dir"]
RESULTS = config["results_dir"]

rule all:
        input:
                RESULTS + "process/sorted_isoforms.xlsx",

rule ids:
        params:
                url = config["gencode_url"],
        output:
                lut = DATA + "gencode.v26.annotation",
        log:
                LOGS + "ids.log"
        benchmark:
                BENCHS + "ids.txt"
        shell:
                "wget -qO- {params.url} | "
                "gzip -dc | "
                "grep -w gene | "
                "cut -f9 | "
                "cut -d ' ' -f2,6 | "
                "tr -d ';\"' "
                "> {output.lut} "
                "2> {log}"

rule mane:
        params:
                url = config["mane_url"],
        output:
                data = DATA + "MANE.csv",
        log:
                LOGS + "mane.log"
        benchmark:
                BENCHS + "mane.txt"
        conda:
                ENVS + "mane.yml"
        script:
                "scripts/mane.py"

rule request:
        input:
                lut = DATA + "gencode.v26.annotation",
        params:
                gene_ids = config["gene_ids"],
                region = config["region"],
        output:
                data = expand(RESULTS + "request/{gene}_message.csv", gene=config["gene_ids"]),
        log:
                LOGS + "request.log"
        benchmark:
                BENCHS + "request.txt"
        threads:
                workflow.cores
        conda:
                ENVS + "request.yml"
        script:
                "scripts/request.py"

rule biomart:
        input:
                data = expand(RESULTS + "request/{gene}_message.csv", gene=config["gene_ids"]),
        output:
                data = expand(RESULTS + "biomart/{gene}_message.csv", gene=config["gene_ids"]),
        log:
                LOGS + "biomart.log"
        benchmark:
                BENCHS + "biomart.txt"
        threads:
                workflow.cores
        conda:
                ENVS + "biomart.yml"
        script:
                "scripts/biomart.py"

rule process:
        input:
                gtex = expand(RESULTS + "request/{gene}_message.csv", gene=config["gene_ids"]),
                bm = expand(RESULTS + "biomart/{gene}_message.csv", gene=config["gene_ids"]),
                mane = DATA + "MANE.csv",
        output:
                data = RESULTS + "process/sorted_isoforms.xlsx",
        log:
                LOGS + "process.log"
        benchmark:
                BENCHS + "process.txt"
        threads:
                workflow.cores
        conda:
                ENVS + "process.yml"
        script:
                "scripts/process.py"
