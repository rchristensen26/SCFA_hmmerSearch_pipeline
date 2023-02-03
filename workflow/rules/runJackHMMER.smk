#Translate to NT to AA with PRODIGAL
rule runProdigal:
    input:
        join(config["ntGenomeDir"], "{strain}.fa")
    output:
        proteinSeqs=join(config["proteinSeqDir"], "{strain}_prodigal.faa")
    conda:
        config["prodigalEnv"]
    shell:
        """
        prodigal -i {input} -a {output.proteinSeqs}
        """

# run jackHMMER search on all genes
rule runJackHMMER:
    input:
        join(config["proteinSeqDir"], "{strain}_prodigal.faa")
    output:
        hmmOut="workflow/out/{gene}/jackhmmer_output/{strain}_{gene}.hmm.out",
        domOut="workflow/out/{gene}/jackhmmer_output/{strain}_{gene}.domtblout",
        msa="workflow/out/{gene}/jackhmmer_output/{strain}_{gene}.sto"
    params:
        seq_file="config/referenceGeneFiles/{gene}.faa"
    conda:
        config["hmmerEnv"]
    shell:
        """
        jackhmmer -o {output.hmmOut} --domtblout {output.domOut} -A {output.msa} {params.seq_file} {input} 
        """

rule parseHMMER:
    input:
        "workflow/out/{gene}/jackhmmer_output/{strain}_{gene}.domtblout"
    output:
        "workflow/out/{gene}/jackhmmer_csv_summary/{strain}_{gene}_hits.csv"
    shell:
        """
        python3 workflow/scripts/parse_hmmer_domtable.py {input} {output}
        """

rule combineCSV:
    input:
        expand("workflow/out/{gene}/jackhmmer_csv_summary/{strain}_{gene}_hits.csv",strain=STRAINS, gene=GENES)
    output:
        "workflow/out/summary_all/jackhmmer_csv_summary/compiled_{gene}_hits.csv"
    params:
        input_dir="workflow/out/{gene}/jackhmmer_csv_summary"
    shell:
        """
        for f in {params.input_dir}/*_hits.csv ; do cat $f ; done > {output}
        """

# convert the HMMER output MSA from stockholm format (.sto) to fasta (.faa)
rule convertMSA_toFASTA:
    input:
        "workflow/out/{gene}/jackhmmer_output/{strain}_{gene}.sto"
    output:
        "workflow/out/{gene}/jackhmmer_faa_summary/{strain}_{gene}_hits.faa"
    shell:
        """
        python3 workflow/scripts/convertFASTA.py {input} {output}
        """

rule combineMSA:
    input:
        expand("workflow/out/{gene}/jackhmmer_faa_summary/{strain}_{gene}_hits.faa",strain=STRAINS, gene=GENES)
    output:
        "workflow/out/summary_all/jackhmmer_faa_summary/compiled_{gene}_hits.faa"
    params:
        input_dir="workflow/out/{gene}/jackhmmer_faa_summary"
    shell:
        """
        for f in {params.input_dir}/*.faa ; do cat $f ; done > {output}
        """