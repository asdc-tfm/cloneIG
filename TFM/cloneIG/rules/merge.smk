rule merge_tools:
    input:
        mixcr_df = f"{OUTDIR}/mixcr/{{patient}}/{{patient}}_mixcr_harmonized.csv",
        vidjil_df = f"{OUTDIR}/vidjil/{{patient}}/{{patient}}_vidjil_harmonized.csv"
    output:
        merge_df = f"{OUTDIR}/merge/{{patient}}/{{patient}}_merge.csv"
    threads:
        get_resource('default', 'threads')
    resources:
        mem = get_resource('default', 'mem'),
        walltime = get_resource('default', 'walltime')
    params:
        warnings_file = f"{OUTDIR}/merge/samples_warning.txt"
    conda:
        f"{ABS_PATH}/envs/python_pkgs.yaml"
    log:
        f"{LOGDIR}/merge/{{patient}}/{{patient}}_merged.log"
    benchmark:
        f"{OUTDIR}/benchmarks/merge_tools/{{patient}}_merge_tools.txt"
    script:
        f"{ABS_PATH}/scripts/merge_tools.py"

rule get_fishplot:
    input:
        merge_df = rules.merge_tools.output.merge_df
    output:
        mixcr_fishplot = f"{OUTDIR}/merge/{{patient}}/{{patient}}_mixcr.html",
        vidjil_fishplot = f"{OUTDIR}/merge/{{patient}}/{{patient}}_vidjil.html"
    threads:
        get_resource('default', 'threads')
    resources:
        mem = get_resource('default', 'mem'),
        walltime = get_resource('default', 'walltime')
    conda:
        f"{ABS_PATH}/envs/fishplot.yaml"
    log:
        f"{LOGDIR}/merge/{{patient}}/{{patient}}_fishplot.log"
    benchmark:
        f"{OUTDIR}/benchmarks/get_fishplot/{{patient}}_get_fishplot.txt"
    script:
        f"{ABS_PATH}/scripts/get_fishplot.R"
