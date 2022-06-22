def get_vidjil_fragments(wildcards):
    R1 = samples.loc[wildcards.patient, wildcards.condition]["R1"]
    R2 = samples.loc[wildcards.patient, wildcards.condition]["R2"]

    return [R1, R2]


rule concat_vidjil_reads:
    input:
        get_vidjil_fragments,
    output:
        concat_reads=temp(f"{OUTDIR}/vidjil/{{patient}}/merge_{{condition}}.fastq.gz"),
    benchmark:
        f"{OUTDIR}/benchmarks/concat_vidjil_reads/{{patient}}_{{condition}}_concat.txt"
    shell:
        "cat {input} > {output}"


def get_vidjil_input(wildcards):
    R1 = samples.loc[wildcards.patient, wildcards.condition]["R1"]
    R2 = samples.loc[wildcards.patient, wildcards.condition]["R2"]

    if pd.isna(R2):
        return R1
    else:
        return f"{OUTDIR}/vidjil/{{patient}}/merge_{{condition}}.fastq.gz"
        

rule vidjil:
    input:
        get_vidjil_input,
    output:
        tvidjil=f"{OUTDIR}/vidjil/{{patient}}/{{patient}}_{{condition}}_vidjil.tsv",
        vvidjil=f"{OUTDIR}/vidjil/{{patient}}/{{patient}}_{{condition}}_vidjil.vidjil",
    threads: get_resource("vidjil", "threads")
    resources:
        mem=get_resource("vidjil", "mem"),
        walltime=get_resource("vidjil", "walltime"),
    params:
        vidjil_name=f"{{patient}}_{{condition}}_vidjil",
        vidjil_dir=directory(f"{OUTDIR}/vidjil/{{patient}}"),
        germline="-g /mnt/tools/vidjil-algo-2021.04/germline/homo-sapiens.g:IGH,IGH+,IGK,IGK+",
        log_file=f"{OUTDIR}/log/vidjil/{{patient}}_{{condition}}_vidjil.log",
        log_error=f"{OUTDIR}/log/vidjil/{{patient}}_{{condition}}_vidjil.error",
    log:
        f"{LOGDIR}/vidjil/{{patient}}/{{condition}}/{{patient}}_{{condition}}_vidjil.log",
    benchmark:
        f"{OUTDIR}/benchmarks/vidjil/{{patient}}_{{condition}}_vidjil.txt"
    shell:
        "mkdir -p {OUTDIR}/log/vidjil; vidjil -o {params.vidjil_dir} -b {params.vidjil_name} {params.germline} {input} > {params.log_file} 2> {params.log_error}"


def get_vidjil_tsv(wildcards):
    conditions = samples.loc[wildcards.patient]["Condition"]
    tsv_list = []

    for c in conditions:
        tsv_list += [
            f"{OUTDIR}/vidjil/{wildcards.patient}/{wildcards.patient}_{c}_vidjil.tsv"
        ]
    return tsv_list


rule filter_vidjil:
    input:
        vidjil_output=get_vidjil_tsv,
    output:
        vidjil_filtered=f"{OUTDIR}/vidjil/{{patient}}/{{patient}}_vidjil_filtered.csv",
    threads: get_resource("default", "threads")
    resources:
        mem=get_resource("default", "mem"),
        walltime=get_resource("default", "walltime"),
    conda:
        f"{ABS_PATH}/envs/python_pkgs.yaml"
    log:
        f"{LOGDIR}/vidjil/{{patient}}/{{patient}}_filter_vidjil.log",
    benchmark:
        f"{OUTDIR}/benchmarks/filter_vidjil/{{patient}}_filter_vidjil.txt"
    script:
        f"{ABS_PATH}/scripts/filter_vidjil.py"


rule join_vidjil_duplicates:
    input:
        filtered=rules.filter_vidjil.output.vidjil_filtered,
    output:
        joined=f"{OUTDIR}/vidjil/{{patient}}/{{patient}}_vidjil_joined.csv",
    threads: get_resource("default", "threads")
    resources:
        mem=get_resource("default", "mem"),
        walltime=get_resource("default", "walltime"),
    conda:
        f"{ABS_PATH}/envs/python_pkgs.yaml"
    log:
        f"{LOGDIR}/vidjil/{{patient}}/{{patient}}_join_vidjil.log",
    benchmark:
        f"{OUTDIR}/benchmarks/join_vidjil_duplicates/{{patient}}_join_vidjil.txt"
    script:
        f"{ABS_PATH}/scripts/join_duplicates.py"


rule harmonize_VDJ_vidjil:
    input:
        joined=rules.join_vidjil_duplicates.output.joined,
    output:
        harmonized=f"{OUTDIR}/vidjil/{{patient}}/{{patient}}_vidjil_harmonized.csv",
    threads: get_resource("default", "threads")
    resources:
        mem=get_resource("default", "mem"),
        walltime=get_resource("default", "walltime"),
    conda:
        f"{ABS_PATH}/envs/python_pkgs.yaml"
    log:
        f"{LOGDIR}/vidjil/{{patient}}/{{patient}}_harmonize_vidjil.log",
    benchmark:
        f"{OUTDIR}/benchmarks/harmonize_VDJ_vidjil/{{patient}}_harmonize_vidjil.txt"
    script:
        f"{ABS_PATH}/scripts/harmonize_VDJ.py"
