def get_mixcr_input(wildcards):
    R1 = samples.loc[wildcards.patient, wildcards.condition]["R1"]
    R2 = samples.loc[wildcards.patient, wildcards.condition]["R2"]

    if pd.isna(R2):
        return R1
    else:
        return [R1, R2]

rule mixcr_align:
    input:
        get_mixcr_input
    output:
        contigs = f"{OUTDIR}/mixcr/{{patient}}/tmp/{{condition}}/{{patient}}_{{condition}}_mixcr.contigs.clns"
    threads:
        get_resource('mixcr', 'threads')
    resources:
        mem = get_resource('mixcr', 'mem'),
        walltime = get_resource('mixcr', 'walltime')
    priority: 1
    params:
        receptor = 'BCR',
        report = f"{OUTDIR}/mixcr/{{patient}}/reports/{{patient}}_{{condition}}_mixcr.report",
        report_path = f"{OUTDIR}/mixcr/{{patient}}/reports",
        mixcr_outname = f"{OUTDIR}/mixcr/{{patient}}/tmp/{{condition}}/{{patient}}_{{condition}}_mixcr",
        extra = '--contig-assembly --only-productive --align "-OallowPartialAlignments=true"',
        log_file = f"{OUTDIR}/log/mixcr/{{patient}}_{{condition}}_mixcr.log",
        log_error = f"{OUTDIR}/log/mixcr/{{patient}}_{{condition}}_mixcr.error"
    log:
        f"{LOGDIR}/mixcr_align/{{patient}}/{{condition}}/{{patient}}_{{condition}}_mixcr_align.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mixcr_align/{{patient}}_{{condition}}_mixcr_align.txt"
    shell: 'mkdir -p {params.report_path} {OUTDIR}/log/mixcr; mixcr -Xmx32g analyze shotgun -s HomoSapiens --starting-material dna --report {params.report} {params.extra} --receptor-type {params.receptor} -t {threads} {input} {params.mixcr_outname} > {params.log_file} 2> {params.log_error}'
    

rule mixcr_export:
    input:
        contigs = rules.mixcr_align.output.contigs
    output:
        fullclones = f"{OUTDIR}/mixcr/{{patient}}/fullClones/{{patient}}_{{condition}}_full_clones.txt"
    threads:
        get_resource('mixcr', 'threads')
    resources:
        mem = get_resource('mixcr', 'mem'),
        walltime = get_resource('mixcr', 'walltime')
    params:
        fullclones_path = f"{OUTDIR}/mixcr/{{patient}}/fullClones/"
    log:
        f"{LOGDIR}/mixcr_export/{{patient}}/{{condition}}/{{patient}}_{{condition}}_mixcr_export.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mixcr_export/{{patient}}_{{condition}}_mixcr_export.txt"
    shell: 'mkdir -p {params.fullclones_path}; mixcr exportClones -c IG -p fullImputed {input.contigs} {output.fullclones}'
    

def get_mixcr_fullclones(wildcards):
    conditions = samples.loc[wildcards.patient]["Condition"]
    fullclones_list = []

    for c in conditions:
        fullclones_list += [f"{OUTDIR}/mixcr/{wildcards.patient}/fullClones/{wildcards.patient}_{c}_full_clones.txt"]

    return fullclones_list


rule filter_mixcr:
    input:
        mixcr_output = get_mixcr_fullclones
    output:
        mixcr_filtered =f"{OUTDIR}/mixcr/{{patient}}/{{patient}}_mixcr_filtered.csv"
    threads:
        get_resource('default', 'threads')
    resources:
        mem = get_resource('default', 'mem'),
        walltime = get_resource('default', 'walltime')
    conda:
        f"{ABS_PATH}/envs/python_pkgs.yaml"
    log:
        f"{LOGDIR}/filter_mixcr/{{patient}}/{{patient}}_filter_mixcr.log"
    benchmark:
        f"{OUTDIR}/benchmarks/filter_mixcr/{{patient}}_filter_mixcr.txt"
    script:
        f"{ABS_PATH}/scripts/filter_mixcr.py"


rule join_mixcr_duplicates:
    input:
        filtered = rules.filter_mixcr.output.mixcr_filtered
    output:
        joined = f"{OUTDIR}/mixcr/{{patient}}/{{patient}}_mixcr_joined.csv"
    threads:
        get_resource('default', 'threads')
    resources:
        mem = get_resource('default', 'mem'),
        walltime = get_resource('default', 'walltime')
    conda:
        f"{ABS_PATH}/envs/python_pkgs.yaml"
    log:
        f"{LOGDIR}/join_mixcr_duplicates/{{patient}}/{{patient}}_join_mixcr.log"
    benchmark:
        f"{OUTDIR}/benchmarks/join_mixcr_duplicates/{{patient}}_join_mixcr.txt"
    script:
        f"{ABS_PATH}/scripts/join_duplicates.py"


rule harmonize_VDJ_mixcr:
    input:
        joined = rules.join_mixcr_duplicates.output.joined
    output:
        harmonized = f"{OUTDIR}/mixcr/{{patient}}/{{patient}}_mixcr_harmonized.csv"
    threads:
        get_resource('default', 'threads')
    resources:
        mem = get_resource('default', 'mem'),
        walltime = get_resource('default', 'walltime')
    conda:
        f"{ABS_PATH}/envs/python_pkgs.yaml"
    log:
        f"{LOGDIR}/harmonize_VDJ_mixcr/{{patient}}/{{patient}}_harmonize_mixcr.log"
    benchmark:
        f"{OUTDIR}/benchmarks/harmonize_VDJ_mixcr/{{patient}}_harmonize_mixcr.txt"
    script:
        f"{ABS_PATH}/scripts/harmonize_VDJ.py"
