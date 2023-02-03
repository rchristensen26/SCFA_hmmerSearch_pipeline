rule combineAllCSV:
    input:
        expand("workflow/out/summary/{pathway}/csv_summary/compiled_{gene}_hits.csv",gene=GENES, pathway=PATHWAYS)
    output:
        "workflow/out/summary/{pathway}/compiled_hits_{pathway}.csv"
    params:
        input_dir="workflow/out/summary/{pathway}/csv_summary"
    shell:
        """
        for f in {params.input_dir}/*.csv ; do cat $f ; done > {output}
        """

rule compileHitInfo:
    input:
        "workflow/out/summary/{pathway}/compiled_hits_{pathway}.csv"
    output:
        "workflow/out/summary/{pathway}/maxHitScoreDF_{pathway}.csv"
    shell:
        """
        python3 workflow/scripts/compileHitInfo.py {input} {output}
        """