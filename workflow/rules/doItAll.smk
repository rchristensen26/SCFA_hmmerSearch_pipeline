# DOWNLOAD ALL STRAIN GENOMES

rule retrieveStrainGenomes:
    output:
        join(config["strainGenomeDir"], "{strain}.fa")
    conda:
        "edirect"
    shell:
        """
        esearch -db assembly -query "{wildcards.strain} [ORGN] AND representative [PROP]" \
        | elink -target nuccore -name assembly_nuccore_refseq \
        | efetch -format fasta > {output}
        """

# BUILD PROFILE HMMS

rule makeMSA:
    input:
        "config/referenceGeneFiles/{pathway}/{gene}.faa"
    output:
        "config/referenceGeneMSA/{pathway}/{gene}_msa.faa"
    shell:
        """
        /usr/local/bin/mafft {input} > {output}
        """

rule buildProfileHMM:
    input:
        "config/referenceGeneMSA/{pathway}/{gene}_msa.faa"
    output:
        "config/profileHMMs/{pathway}/{gene}.HMM"
    conda:
        config["hmmerEnv"]
    shell:
        """
        hmmbuild {output} {input}
        """

# RUN HMMER SEARCH
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
        prodigal=join(config["proteinSeqDir"], "{strain}_prodigal.faa"),
        hmm_profile="config/profileHMMs/{pathway}/{gene}.HMM"
    output:
        hmmOut=(join(config["hmmOutDir"],"{pathway}/{gene}/hmmer_output/{strain}_{gene}.hmm.out")),
        domOut=(join(config["hmmOutDir"],"{pathway}/{gene}/hmmer_output/{strain}_{gene}.domtblout")),
        msa=(join(config["hmmOutDir"],"{pathway}/{gene}/hmmer_output/{strain}_{gene}.sto"))
    # params:
    #     hmm_profile="config/profileHMMs/{pathway}/{gene}.HMM"
    conda:
        config["hmmerEnv"]
    shell:
        """
        hmmsearch -o {output.hmmOut} --domtblout {output.domOut} -A {output.msa} {input.hmm_profile} {input.prodigal} 
        """

# run jackHMMER search on all genes
# rule runJackHMMER:
#     input:
#         prodigal=join(config["proteinSeqDir"], "{strain}_prodigal.faa"),
#         seq_file="config/referenceGeneFiles/P3/Roseburia-inulivorans-pduCDE.faa"
#     output:
#         hmmOut="workflow/out/hmmer_search/{pathway}/{jackhammer_gene}/hmmer_output/{strain}_{jackhammer_gene}.hmm.out",
#         domOut="workflow/out/hmmer_search/{pathway}/{jackhammer_gene}/hmmer_output/{strain}_{jackhammer_gene}.domtblout",
#         msa="workflow/out/hmmer_search/{pathway}/{jackhammer_gene}/hmmer_output/{strain}_{jackhammer_gene}.sto"
#     params:
#         seq_file="config/referenceGeneFiles/{pathway}/{jackhammer_gene}.faa"
#     conda:
#         config["hmmerEnv"]
#     shell:
#         """
#         jackhmmer -o {output.hmmOut} --domtblout {output.domOut} -A {output.msa} {input.seq_file} {input.prodigal}
#         """

rule parseHMMER:
    input:
        (join(config["hmmOutDir"],"{pathway}/{gene}/hmmer_output/{strain}_{gene}.domtblout"))
    output:
        (join(config["hmmOutDir"],"{pathway}/{gene}/csv_summary/{strain}_{gene}_hits.csv"))
    shell:
        """
        python3 workflow/scripts/parse_hmmer_domtable.py {input} {output}
        """

rule combineCSV:
    # input:
    #     expand("workflow/out/hmmer_search/{pathway}/{gene}/csv_summary/{strain}_{gene}_hits.csv",strain=STRAINS, gene=GENES, pathway=PATHWAYS)
    output:
        "workflow/out/summary/{pathway}/csv_summary/compiled_{gene}_hits.csv"
    params:
        input_dir=(join(config["hmmOutDir"],"{pathway}/{gene}/csv_summary"))
    shell:
        """
        for f in {params.input_dir}/*_hits.csv ; do cat $f ; done > {output}
        """

# convert the HMMER output MSA from stockholm format (.sto) to fasta (.faa)
rule convertMSA_toFASTA:
    input:
        (join(config["hmmOutDir"],"{pathway}/{gene}/hmmer_output/{strain}_{gene}.sto"))
    output:
        (join(config["hmmOutDir"],"{pathway}/{gene}/faa_summary/{strain}_{gene}_hits.faa"))
    shell:
        """
        python3 workflow/scripts/convertFASTA.py {input} {output}
        """

rule combineMSA:
    # input:
        # expand("workflow/out/hmmer_search/{pathway}/{gene}/faa_summary/{strain}_{gene}_hits.faa",strain=STRAINS, gene=GENES, pathway=PATHWAYS)
    output:
        "workflow/out/summary/{pathway}/faa_summary/compiled_{gene}_hits.faa"
    params:
        input_dir=(join(config["hmmOutDir"],"{pathway}/{gene}/faa_summary"))
    shell:
        """
        for f in {params.input_dir}/*.faa ; do cat $f ; done > {output}
        """

# COMPILE HITS

rule combineAllCSV:
    # input:
    #     expand("workflow/out/summary/{pathway}/csv_summary/compiled_{gene}_hits.csv",gene=GENES, pathway=PATHWAYS),
    #     expand("workflow/out/summary/{pathway}/csv_summary/compiled_{gene}_hits.csv", gene=JACKHAMMER_GENES, pathway=PATHWAYS),
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