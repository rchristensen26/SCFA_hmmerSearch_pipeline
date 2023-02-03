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
