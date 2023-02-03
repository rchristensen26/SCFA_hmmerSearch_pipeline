rule retrieveStrainGenomes:
    output:
        "config/strain_genomes/{strain}.fa"
    conda:
        "edirect"
    shell:
        """
        esearch -db assembly -query "{wildcards.strain} [ORGN] AND representative [PROP]" \
        | elink -target nuccore -name assembly_nuccore_refseq \
        | efetch -format fasta > {output}
        """
