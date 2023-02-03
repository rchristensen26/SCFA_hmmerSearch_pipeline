#Translate to NT to AA with PRODIGAL
rule runProdigal:
    input:
        join(config["strainGenomeDir"], "{strain}.fa")
    output:
        proteinSeqs=join(config["proteinSeqDir"], "{strain}_prodigal.faa")
    conda:
        config["prodigalEnv"]
    shell:
        """
        prodigal -i {input} -a {output.proteinSeqs}
        """

# run HMMER search on all genes
rule runHMMER:
    input:
        join(config["proteinSeqDir"], "{strain}_prodigal.faa")
    output:
        hmmOut="workflow/out/hmmer_search/{pathway}/{gene}/hmmer_output/{strain}_{gene}.hmm.out",
        domOut="workflow/out/hmmer_search/{pathway}/{gene}/hmmer_output/{strain}_{gene}.domtblout",
        msa="workflow/out/hmmer_search/{pathway}/{gene}/hmmer_output/{strain}_{gene}.sto"
    params:
        hmm_profile="config/profileHMMs/{pathway}/{gene}.HMM"
    conda:
        config["hmmerEnv"]
    shell:
        """
        hmmsearch -o {output.hmmOut} --domtblout {output.domOut} -A {output.msa} {params.hmm_profile} {input} 
        """

rule parseHMMER:
    input:
        "workflow/out/hmmer_search/{pathway}/{gene}/hmmer_output/{strain}_{gene}.domtblout"
    output:
        "workflow/out/hmmer_search/{pathway}/{gene}/csv_summary/{strain}_{gene}_hits.csv"
    shell:
        """
        python3 workflow/scripts/parse_hmmer_domtable.py {input} {output}
        """

rule combineCSV:
    input:
        expand("workflow/out/hmmer_search/{pathway}/{gene}/csv_summary/{strain}_{gene}_hits.csv",strain=STRAINS, gene=GENES, pathway=PATHWAYS)
    output:
        "workflow/out/summary/{pathway}/csv_summary/compiled_{gene}_hits.csv"
    params:
        input_dir="workflow/out/hmmer_search/{pathway}/{gene}/csv_summary"
    shell:
        """
        for f in {params.input_dir}/*_hits.csv ; do cat $f ; done > {output}
        """

# convert the HMMER output MSA from stockholm format (.sto) to fasta (.faa)
rule convertMSA_toFASTA:
    input:
        "workflow/out/hmmer_search/{pathway}/{gene}/hmmer_output/{strain}_{gene}.sto"
    output:
        "workflow/out/hmmer_search/{pathway}/{gene}/faa_summary/{strain}_{gene}_hits.faa"
    shell:
        """
        python3 workflow/scripts/convertFASTA.py {input} {output}
        """

rule combineMSA:
    input:
        expand("workflow/out/hmmer_search/{pathway}/{gene}/faa_summary/{strain}_{gene}_hits.faa",strain=STRAINS, gene=GENES, pathway=PATHWAYS)
    output:
        "workflow/out/summary/{pathway}/faa_summary/compiled_{gene}_hits.faa"
    params:
        input_dir="workflow/out/hmmer_search/{pathway}/{gene}/faa_summary"
    shell:
        """
        for f in {params.input_dir}/*.faa ; do cat $f ; done > {output}
        """
